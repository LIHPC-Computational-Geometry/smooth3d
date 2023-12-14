/*
 * ObjectiveOrthogonal2D.h
 *
 *  Created on: 6 juil. 2015
 *      Author: weilljc
 */

#ifndef OPT2D_OBJECTIVEORTHOGONAL2D_H_
#define OPT2D_OBJECTIVEORTHOGONAL2D_H_

#include "IObjective3D.h"
#include <gmds/ig/Node.h>
#include <gmds/ig/Face.h>

namespace Smooth3D {

class ObjectiveOrthogonal2D: public IObjective3D {
public:
  ObjectiveOrthogonal2D();
  ObjectiveOrthogonal2D(gmds::Node node);

  virtual ~ObjectiveOrthogonal2D();
  void setNode(gmds::Node& node);
  virtual Real evalF(const Real3 & p);

  virtual void evalDF(const Real3 & p, Real3 & grad);

  virtual Real evalFDF(const Real3 & p, Real3 & grad);

private:

  gmds::Node m_node;
};

} /* namespace Smooth3D */

#endif /* OPT2D_OBJECTIVEORTHOGONAL2D_H_ */
