/*
 * ObjectiveCondition2D.h
 *
 *  Created on: 13 nov. 2013
 *      Author: weilljc
 */

#ifndef OBJECTIVECONDITION2D_H_
#define OBJECTIVECONDITION2D_H_

#include "IObjective3D.h"
#include "gmds/ig/Node.h"
#include "gmds/ig/Face.h"

namespace Smooth3D {

class ObjectiveCondition2D: public IObjective3D {

public:

	ObjectiveCondition2D();
	ObjectiveCondition2D(gmds::Node node);
	virtual ~ObjectiveCondition2D();

	void set_node(gmds::Node& node);

	virtual Real evalF(const Real3 & p);

	virtual void evalDF(const Real3 & p, Real3 & grad);

	virtual Real evalFDF(const Real3 & p, Real3 & grad);

private:

	gmds::Node m_node;

};

} /* namespace Smooth3D */

#endif /* OBJECTIVECONDITION2D_H_ */
