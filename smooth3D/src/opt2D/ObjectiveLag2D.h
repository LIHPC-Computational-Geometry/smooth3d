/*
 * ObjectiveLag2D.h
 *
 *  Created on: 13 nov. 2013
 *      Author: weilljc
 */

#ifndef OBJECTIVELAG2D_H_
#define OBJECTIVELAG2D_H_

#include "IObjective3D.h"
#include "math/Optimize.h"
#include "math/LagrangeMapping.h"


namespace Smooth3D {

class ObjectiveLag2D :public Optimize::Fonction < Real2>{
private:

  LagrangeMapping2D * m_lag_map;
  IObjective3D * m_fdf;

public:
  ObjectiveLag2D(IObjective3D * fonction, LagrangeMapping2D * lag) : m_lag_map(lag), m_fdf(fonction)
  {}
  virtual ~ObjectiveLag2D()
  {}

  virtual Real
   evalF (const Real2 & vec);


   virtual void
   evalDF (const Real2 & vec, Real2 & gradient);


   virtual Real
   evalFDF (const Real2 & vec, Real2 & gradient);

};

} /* namespace Smooth3D */

#endif /* OBJECTIVELAG2D_H_ */
