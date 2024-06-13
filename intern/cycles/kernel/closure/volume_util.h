/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

CCL_NAMESPACE_BEGIN

#define RAD2DEGF(_rad) ((_rad) * (float)(180.0 / M_PI))
#define DEG2RADF(_deg) ((_deg) * (float)(M_PI / 180.0))
#define ANGLE_EPSILON 1e-6f
#define ANGLE_TOL 1.0f / 180.0f
#define CDF_TOL 0.01f

ccl_device float fournier_forand_sigma(float n, float sin_theta_2)
{
  float u = 4.0f * sqr(sin_theta_2);
  return u / (3.0f * sqr(n - 1.0f));
};

ccl_device float fournier_forand_cdf(float theta, float n, float B)
{
  if (theta <= 0.0f) {
    return 0.0f;
  }
  else if (theta >= M_PI_F) {
    return 1.0f;
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

ccl_device float interpolate_linear(float ax, float ay, float bx, float by, float x)
{
  return ay + (by - ay) * ((x - ax) / (bx - ax));
}

ccl_device float find_fournier_forand_angle(float rand, float B, float n)
{
  int it = 0;
  float l_low, l_up, m = 0.0f;
  float theta = 0.8726646259971648f;  // 50 degrees
  float fm = fournier_forand_cdf(theta, n, B);
  float err = fm - rand;
  if (err < 0.0f) {
    l_low = theta;
    l_up = M_PI_F;
  }
  else if (err > 0.0f) {
    l_low = 0.0f;
    l_up = theta;
  }
  else {
    return theta;
  }

  while (it < 100 && (fabsf(l_low - l_up) > ANGLE_TOL || fabsf(err) > CDF_TOL)) {
    m = (l_low + l_up) / 2.0f;
    fm = fournier_forand_cdf(m, n, B);
    err = fm - rand;
    it += 1;

    if (signf(fournier_forand_cdf(l_low, n, B) - rand) == signf(err)) {
      l_low = m;
    }
    else if (signf(fournier_forand_cdf(l_up, n, B) - rand) == signf(err)) {
      l_up = m;
    }
  }
  m = (l_low + l_up) / 2.0f;
  fm = fournier_forand_cdf(m, n, B);
  err = fm - rand;
  if (err < 0.0f) {
    return interpolate_linear(fm, m, fournier_forand_cdf(l_up, n, B), l_up, rand);
  }
  else if (err > 0.0f) {
    return interpolate_linear(fournier_forand_cdf(l_low, n, B), l_low, fm, m, rand);
  }
  return m;
}

CCL_NAMESPACE_END
