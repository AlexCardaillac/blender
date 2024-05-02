/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

#include "util/array.h"

/* Scattering phase functions */
typedef enum PhaseFunction {
  PHASE_HENYEY_GREENSTEIN = 0,
  PHASE_FOURNIER_FORAND = 1,
} PhaseFunction;

CCL_NAMESPACE_BEGIN

/* VOLUME EXTINCTION */

ccl_device void volume_extinction_setup(ccl_private ShaderData *sd, Spectrum weight)
{
  if (sd->flag & SD_EXTINCTION) {
    sd->closure_transparent_extinction += weight;
  }
  else {
    sd->flag |= SD_EXTINCTION;
    sd->closure_transparent_extinction = weight;
  }
}

/* HENYEY-GREENSTEIN CLOSURE */

typedef struct HenyeyGreensteinVolume {
  SHADER_CLOSURE_BASE;

  float g;
  float IoR;
  float B;
  int16_t phase;
} HenyeyGreensteinVolume;

static_assert(sizeof(ShaderClosure) >= sizeof(HenyeyGreensteinVolume),
              "HenyeyGreensteinVolume is too large!");

#define ANGLE_EPSILON 1e-6f

/* Given cosine between rays, return probability density that a photon bounces
 * to that direction. The g parameter controls how different it is from the
 * uniform sphere. g=0 uniform diffuse-like, g=1 close to sharp single ray. */
ccl_device float single_peaked_henyey_greenstein(float cos_theta, float g)
{
  return ((1.0f - g * g) / safe_powf(1.0f + g * g - 2.0f * g * cos_theta, 1.5f)) *
         (M_1_PI_F * 0.25f);
};

ccl_device float ff_sigma(float n, float theta)
{
  float u = 4.0f * powf(sinf(theta / 2.0f), 2.0f);
  return u / (3.0f * powf(n - 1.0f, 2.0f));
};

/* Given cosine between rays, return probability density that a photon bounces
 * to that direction. The n parameter is the particle index of refraction and
 * controls how much of the light is refracted. B is the particle backscatter
 * fraction, B = b_b / b. */
ccl_device float single_peaked_fournier_forand(float cos_theta, float n, float B)
{
  float theta = acosf(cos_theta);
  if (theta < ANGLE_EPSILON)
  {
    theta = ANGLE_EPSILON;
  }
  float s90 = ff_sigma(n, M_PI_2_F);
  float s180 = ff_sigma(n, M_PI_F);
  float s_theta = ff_sigma(n, theta);
  float slope = 2.0f * (logf(2.0f * B * (s90 - 1.0f) + 1.0f) /
                logf(s90)) + 3.0f;
  float v = (3.0f - slope) / 2.0f;
  float pow_stheta_v = powf(s_theta, v);
  float pow_s180_v = powf(s180, v);
  float pf = 1.0f / (4.0f * M_PI_F * powf(1.0f - s_theta, 2.0f) * pow_stheta_v);
  pf *= (v * (1.0f - s_theta) - (1.0f - pow_stheta_v) + (s_theta * (1.0f -
        pow_stheta_v) - v * (1.0f - s_theta)) * powf(sinf(theta / 2.0f),-2.0f));
  pf += ((1.0f - pow_s180_v) / (16.0f * M_PI_F * (s180 - 1.0f) * pow_s180_v)) *
        (3.0f * powf(cosf(theta), 2.0f) - 1.0f);
  return pf;
};


#define CDF_RESOLUTION 1800

typedef float CDFTable[CDF_RESOLUTION][2];

ccl_device CDFTable &get_fournier_forand_cdf_table()
{
  static CDFTable ff_cdf = {0.0};
  return ff_cdf;
}

ccl_device void create_fournier_forand_cdf_table(ccl_private HenyeyGreensteinVolume *volume)
{
  static float IoR = -1.0f;
  static float B = -1.0f;

  if (volume->B == B && volume->IoR == IoR)
  {
    return;
  }
  B = volume->B;
  IoR = volume->IoR;
  CDFTable &ff_cdf = get_fournier_forand_cdf_table();
  ff_cdf[0][0] = 0.0;
  ff_cdf[0][1] = 1.0;
}

ccl_device int volume_henyey_greenstein_setup(ccl_private HenyeyGreensteinVolume *volume)
{
  volume->type = CLOSURE_VOLUME_HENYEY_GREENSTEIN_ID;

  if (volume->phase == PHASE_FOURNIER_FORAND) {
    // printf("--------------\nPHASE_FOURNIER_FORAND\n--------------\n");
      // *pdf = single_peaked_fournier_forand(cos_theta, svc->IoR, svc->B);
  }
  else {
    // printf("--------------\nPHASE_HENYEY_GREENSTEIN\n--------------\n");

    /* clamp anisotropy to avoid delta function */
    volume->g = signf(volume->g) * min(fabsf(volume->g), 1.0f - 1e-3f);
  }

  return SD_SCATTER;
}

ccl_device Spectrum volume_henyey_greenstein_eval_phase(ccl_private const ShaderVolumeClosure *svc,
                                                        const float3 wi,
                                                        float3 wo,
                                                        ccl_private float *pdf)
{
  float g = svc->g;

  /* note that wi points towards the viewer */
  if (fabsf(g) < 1e-3f) {
    *pdf = M_1_PI_F * 0.25f;
  }
  else {
    float cos_theta = dot(-wi, wo);
    if (svc->phase == PHASE_FOURNIER_FORAND) {
      *pdf = single_peaked_fournier_forand(cos_theta, svc->IoR, svc->B);
    }
    else {
      *pdf = single_peaked_henyey_greenstein(cos_theta, g);
    }
  }

  return make_spectrum(*pdf);
}

ccl_device float3 henyey_greenstrein_sample(float3 D, float g, float2 rand, ccl_private float *pdf)
{
  /* match pdf for small g */
  float cos_theta;
  bool isotropic = fabsf(g) < 1e-3f;

  if (isotropic) {
    cos_theta = (1.0f - 2.0f * rand.x);
    if (pdf) {
      *pdf = M_1_PI_F * 0.25f;
    }
  }
  else {
    float k = (1.0f - g * g) / (1.0f - g + 2.0f * g * rand.x);
    cos_theta = (1.0f + g * g - k * k) / (2.0f * g);
    if (pdf) {
      // if (svc->phase == PHASE_FOURNIER_FORAND) {
      // *pdf = single_peaked_fournier_forand(cos_theta, svc->IoR, svc->B);
      // }
      // else {
        *pdf = single_peaked_henyey_greenstein(cos_theta, g);
      // }
    }
  }

  float sin_theta = sin_from_cos(cos_theta);
  float phi = M_2PI_F * rand.y;
  float3 dir = make_float3(sin_theta * cosf(phi), sin_theta * sinf(phi), cos_theta);

  float3 T, B;
  make_orthonormals(D, &T, &B);
  dir = dir.x * T + dir.y * B + dir.z * D;

  return dir;
}

ccl_device int volume_henyey_greenstein_sample(ccl_private const ShaderVolumeClosure *svc,
                                               float3 wi,
                                               float2 rand,
                                               ccl_private Spectrum *eval,
                                               ccl_private float3 *wo,
                                               ccl_private float *pdf)
{
  float g = svc->g;

  /* note that wi points towards the viewer and so is used negated */
  *wo = henyey_greenstrein_sample(-wi, g, rand, pdf);
  *eval = make_spectrum(*pdf); /* perfect importance sampling */

  return LABEL_VOLUME_SCATTER;
}

/* VOLUME CLOSURE */

ccl_device Spectrum volume_phase_eval(ccl_private const ShaderData *sd,
                                      ccl_private const ShaderVolumeClosure *svc,
                                      float3 wo,
                                      ccl_private float *pdf)
{
  return volume_henyey_greenstein_eval_phase(svc, sd->wi, wo, pdf);
}

ccl_device int volume_phase_sample(ccl_private const ShaderData *sd,
                                   ccl_private const ShaderVolumeClosure *svc,
                                   float2 rand,
                                   ccl_private Spectrum *eval,
                                   ccl_private float3 *wo,
                                   ccl_private float *pdf)
{
  return volume_henyey_greenstein_sample(svc, sd->wi, rand, eval, wo, pdf);
}

/* Volume sampling utilities. */

/* todo: this value could be tweaked or turned into a probability to avoid
 * unnecessary work in volumes and subsurface scattering. */
#define VOLUME_THROUGHPUT_EPSILON 1e-6f

ccl_device Spectrum volume_color_transmittance(Spectrum sigma, float t)
{
  return exp(-sigma * t);
}

ccl_device float volume_channel_get(Spectrum value, int channel)
{
  return GET_SPECTRUM_CHANNEL(value, channel);
}

ccl_device int volume_sample_channel(Spectrum albedo,
                                     Spectrum throughput,
                                     float rand,
                                     ccl_private Spectrum *pdf)
{
  /* Sample color channel proportional to throughput and single scattering
   * albedo, to significantly reduce noise with many bounce, following:
   *
   * "Practical and Controllable Subsurface Scattering for Production Path
   *  Tracing". Matt Jen-Yuan Chiang, Peter Kutz, Brent Burley. SIGGRAPH 2016. */
  Spectrum weights = fabs(throughput * albedo);
  float sum_weights = reduce_add(weights);

  if (sum_weights > 0.0f) {
    *pdf = weights / sum_weights;
  }
  else {
    *pdf = make_spectrum(1.0f / SPECTRUM_CHANNELS);
  }

  float pdf_sum = 0.0f;
  FOREACH_SPECTRUM_CHANNEL (i) {
    pdf_sum += GET_SPECTRUM_CHANNEL(*pdf, i);
    if (rand < pdf_sum) {
      return i;
    }
  }
  return SPECTRUM_CHANNELS - 1;
}

CCL_NAMESPACE_END
