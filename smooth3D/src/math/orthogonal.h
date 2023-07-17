/*
 * orthogonal.h
 *
 *  Created on: 6 juil. 2015
 *      Author: weilljc
 */

#ifndef MATH_ORTHOGONAL_H_
#define MATH_ORTHOGONAL_H_

#include "Real3.h"

namespace Smooth3D {

inline Real orthogonal2D(Real3 v1, Real3 v2, Real3 v3, Real target) {
  Real3 u = v2 - v1, v = v3 - v1;
  Real cos = math::scaMul(u, v) / std::sqrt(u.abs2() * v.abs2());
  return (cos - target) * (cos - target);
}

inline Real3 gradOrtho2DV1(Real3 v1, Real3 v2, Real3 v3, Real target) {
  Real3 u = v2 - v1, v = v3 - v1;
  Real dot = math::scaMul(u, v);
  Real prod = u.abs2() * v.abs2();
  Real racine = std::sqrt(prod);

  Real cos = dot / racine;
  Real den = u.abs2() * v.abs2();

  Real3 grad = 2. * cos
      * ((-u - v) / racine
	  - dot * ((-v * u.abs2() - u * v.abs2()) / (prod * racine)));
  return grad;
}

inline Real3 gradOrtho2DV2(Real3 v1, Real3 v2, Real3 v3) {
  Real3 u = v2 - v1, v = v3 - v1;
  Real dot = math::scaMul(u, v);
  Real num = (dot * dot);
  Real den = u.abs2() * v.abs2();
  Real3 grad_num = 2.0 * (v3 - v1) * dot;
  Real3 grad_den = 2.0 * (v2 - v1) * v.abs2();
  return (grad_num - grad_den * num / den) / den;
}

}

#endif /* MATH_ORTHOGONAL_H_ */
