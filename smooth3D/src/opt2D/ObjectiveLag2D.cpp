/*
 * ObjectiveLag2D.cpp

 *
 *  Created on: 13 nov. 2013
 *      Author: weilljc
 *
 */

#include "ObjectiveLag2D.h"

namespace Smooth3D {

Real ObjectiveLag2D::evalF(const Real2 & uv) {
  return m_fdf->evalF(m_lag_map->eval(uv));
}

void ObjectiveLag2D::evalDF(const Real2 & uv, Real2 & gradient) {
  Real3 real_grad, dmdu, dmdv;

  m_fdf->evalDF(m_lag_map->eval(uv), real_grad);
  m_lag_map->gradient(uv, dmdu, dmdv);
  gradient = Real2(math::scaMul(real_grad, dmdu),
      math::scaMul(real_grad, dmdv));
}

Real ObjectiveLag2D::evalFDF(const Real2 & uv, Real2 & gradient) {
  Real3 real_grad, dmdu, dmdv;

  Real f = m_fdf->evalFDF(m_lag_map->eval(uv), real_grad);
  m_lag_map->gradient(uv, dmdu, dmdv);
  gradient = Real2(math::scaMul(real_grad, dmdu),
      math::scaMul(real_grad, dmdv));
  return f;
}

}
