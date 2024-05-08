/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

CCL_NAMESPACE_BEGIN

#define DEG2RADF(_deg) ((_deg) * (float)(M_PI / 180.0))
#define ANGLE_EPSILON 1e-6f
#define CDF_RESOLUTION 4000
#define STEP_MIN 1e-6f

typedef float CDFTable[CDF_RESOLUTION][2];

ccl_device CDFTable &get_fournier_forand_cdf_table()
{
  static CDFTable ff_cdf = {0.0};
  return ff_cdf;
}

ccl_device float logistic_fn(float theta, float B)
{
  return 1.0f / (1.0f + expf(-0.04f * (0.5f - B) * (theta - 90.0f)));
}

ccl_device float logistic_norm_fn(float theta, float B)
{
  float y = logistic_fn(theta, B) - logistic_fn(0.0f, B);
  y /= logistic_fn(180.0f, B) - logistic_fn(0.0f, B);
  return y;
}

ccl_device float fournier_forand_sigma(float n, float sin_theta_2)
{
  float u = 4.0f * sqr(sin_theta_2 * sin_theta_2);
  return u / (3.0f * sqr(n - 1.0f));
};

ccl_device float fournier_forand_cdf(float theta, float n, float B)
{
  if (theta == 0.0f) {
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
  float cdf = 1.0f / ((1.0f - s_theta) * pow_stheta_v);
  cdf *= ((1.0f - powf(s_theta, v + 1.0f)) - (1.0f - pow_stheta_v) * sqr(sin_theta_2));
  cdf += (1.0f / 8.0f) * ((1.0f - pow_s180_v) / ((s180 - 1.0f) * pow_s180_v)) * cosf(theta) *
         sqr(sinf(theta));
  return cdf;
}

ccl_device uint create_fournier_forand_cdf_table(float n, float b)
{
  static float IoR = -1.0f;
  static float B = -1.0f;

  if (b == B && n == IoR) {
    return 0;
  }
  B = b;
  IoR = n;
  CDFTable &ff_cdf = get_fournier_forand_cdf_table();

  float l0 = logistic_fn(0.0f, B);
  float norm = logistic_fn(180.0f, B) - l0;

  float step = STEP_MIN;
  float theta = step;
  ff_cdf[0][1] = 0.0f;
  ff_cdf[0][1] = fournier_forand_cdf(0.0f, IoR, B);
  uint i = 1;
  while (theta < 180.0f && i < CDF_RESOLUTION) {
    ff_cdf[i][0] = DEG2RADF(theta);
    ff_cdf[i][1] = fournier_forand_cdf(DEG2RADF(theta), IoR, B);
    step = (logistic_fn(theta, B) - l0) / norm;
    if (step < STEP_MIN) {
      step = STEP_MIN;
    }
    theta += step;
    i += 1;
  }
  if (i < CDF_RESOLUTION) {
    ff_cdf[i][0] = theta;
    ff_cdf[i][1] = 0.0f;
  }
  return i;
}

CCL_NAMESPACE_END
