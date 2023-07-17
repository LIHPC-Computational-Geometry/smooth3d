/*
 * Objective3D.h
 *
 *  Created on: 13 nov. 2013
 *      Author: weilljc
 */

#ifndef OBJECTIVE3D_H_
#define OBJECTIVE3D_H_

#include "math/Optimize.h"

namespace Smooth3D {

class IObjective3D: public Optimize::Fonction<Real3> {

public:
  IObjective3D() {
  }
  virtual ~ IObjective3D() {
  }

  virtual Real
  evalF(const Real3 & vec) = 0;
  virtual void
  evalDF(const Real3 & vec, Real3 & gradient) = 0;
  virtual Real
  evalFDF(const Real3 & vec, Real3 & gradient) = 0;
};

} /* namespace Smooth3D */

#endif /* OBJECTIVE3D_H_ */
