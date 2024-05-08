/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

#include "kernel/closure/volume_util.h"
#include "util/array.h"

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

/* VOLUME SCATTERING CLOSURE */

typedef struct ScatteringVolume {
  SHADER_CLOSURE_BASE;

  float g;
  float IoR;
  float B;
  int16_t phase;
} ScatteringVolume;

static_assert(sizeof(ShaderClosure) >= sizeof(ScatteringVolume), "ScatteringVolume is too large!");

/* FOURNIER-FORAND PHASE */

/* Given cosine between rays, return probability density that a photon bounces
 * to that direction. The n parameter is the particle index of refraction and
 * controls how much of the light is refracted. B is the particle backscatter
 * fraction, B = b_b / b. */
ccl_device float single_peaked_fournier_forand(float cos_theta, float n, float B)
{
  float theta = acosf(cos_theta);
  if (theta < ANGLE_EPSILON) {
    theta = ANGLE_EPSILON;
  }
  float s90 = fournier_forand_sigma(n, sinf(M_PI_2_F / 2.0f));
  float s180 = fournier_forand_sigma(n, sinf(M_PI_F / 2.0f));
  float sin_theta_2 = sinf(theta / 2.0f);
  float s_theta = fournier_forand_sigma(n, sin_theta_2);
  float slope = 2.0f * (logf(2.0f * B * (s90 - 1.0f) + 1.0f) / logf(s90)) + 3.0f;
  float v = (3.0f - slope) / 2.0f;
  float pow_stheta_v = powf(s_theta, v);
  float pow_s180_v = powf(s180, v);
  float pf = 1.0f / (4.0f * M_PI_F * sqr(1.0f - s_theta) * pow_stheta_v);
  pf *= (v * (1.0f - s_theta) - (1.0f - pow_stheta_v) +
         (s_theta * (1.0f - pow_stheta_v) - v * (1.0f - s_theta)) * (1.0f / sqr(sin_theta_2)));
  pf += ((1.0f - pow_s180_v) / (16.0f * M_PI_F * (s180 - 1.0f) * pow_s180_v)) *
        (3.0f * sqr(cosf(theta)) - 1.0f);
  return pf;
};

// TO DO
ccl_device float3 fournier_forand_sample(float3 D, float g, float2 rand, ccl_private float *pdf)
{
  return make_float3(0.0f, 0.0f, 0.0f);
}

/* HENYEY-GREENSTEIN PHASE */

/* Given cosine between rays, return probability density that a photon bounces
 * to that direction. The g parameter controls how different it is from the
 * uniform sphere. g=0 uniform diffuse-like, g=1 close to sharp single ray. */
ccl_device float single_peaked_henyey_greenstein(float cos_theta, float g)
{
  return ((1.0f - g * g) / safe_powf(1.0f + g * g - 2.0f * g * cos_theta, 1.5f)) *
         (M_1_PI_F * 0.25f);
};

ccl_device float3 henyey_greenstein_sample(float3 D, float g, float2 rand, ccl_private float *pdf)
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
      *pdf = single_peaked_henyey_greenstein(cos_theta, g);
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

/* VOLUME CLOSURE */

ccl_device Spectrum volume_phase_eval(ccl_private const ShaderData *sd,
                                      ccl_private const ShaderVolumeClosure *svc,
                                      float3 wo,
                                      ccl_private float *pdf)
{
  /* note that wi points towards the viewer */
  float cos_theta;

  switch (svc->phase) {
    case CLOSURE_VOLUME_FOURNIER_FORAND_ID:
      cos_theta = dot(-sd->wi, wo);
      *pdf = single_peaked_fournier_forand(cos_theta, svc->IoR, svc->B);
      break;
    default:  // CLOSURE_VOLUME_HENYEY_GREENSTEIN_ID
      float g = svc->g;

      if (fabsf(g) < 1e-3f) {
        *pdf = M_1_PI_F * 0.25f;
      }
      else {
        cos_theta = dot(-sd->wi, wo);
        *pdf = single_peaked_henyey_greenstein(cos_theta, g);
      }
      break;
  }
  return make_spectrum(*pdf);
}

ccl_device int volume_phase_sample(ccl_private const ShaderData *sd,
                                   ccl_private const ShaderVolumeClosure *svc,
                                   float2 rand,
                                   ccl_private Spectrum *eval,
                                   ccl_private float3 *wo,
                                   ccl_private float *pdf)
{
  /* note that wi points towards the viewer and so is used negated */
  switch (svc->phase) {
    case CLOSURE_VOLUME_FOURNIER_FORAND_ID:
      *wo = fournier_forand_sample(-sd->wi, svc->g, rand, pdf);
      break;
    default:  // CLOSURE_VOLUME_HENYEY_GREENSTEIN_ID
      *wo = henyey_greenstein_sample(-sd->wi, svc->g, rand, pdf);
      break;
  }
  *eval = make_spectrum(*pdf); /* perfect importance sampling */
  return LABEL_VOLUME_SCATTER;
}

ccl_device int volume_phase_setup(ccl_private ScatteringVolume *volume)
{
  volume->type = CLOSURE_VOLUME_HENYEY_GREENSTEIN_ID;

  switch (volume->phase) {
    case CLOSURE_VOLUME_FOURNIER_FORAND_ID:
      /* clamp backscatter fraction to avoid delta function */
      volume->B = min(fabsf(volume->B), 0.5f - 1e-3f);
      break;
    default:  // CLOSURE_VOLUME_HENYEY_GREENSTEIN_ID
      /* clamp anisotropy to avoid delta function */
      volume->g = signf(volume->g) * min(fabsf(volume->g), 1.0f - 1e-3f);
      break;
  }

  return SD_SCATTER;
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
