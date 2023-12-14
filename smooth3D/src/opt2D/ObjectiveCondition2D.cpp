/*
 * ObjectiveCondition2D.cpp
 *
 *  Created on: 7 juil. 2014
 *      Author: ledouxf
 */


#include "ObjectiveCondition2D.h"

#include "math/condition.h"

namespace Smooth3D {

ObjectiveCondition2D::ObjectiveCondition2D(){
	m_node= gmds::Node();
}

ObjectiveCondition2D::ObjectiveCondition2D(gmds::Node node) : m_node(node){;}

ObjectiveCondition2D::~ObjectiveCondition2D() {}

void ObjectiveCondition2D::set_node(gmds::Node& node) {
  m_node = node;
}

Real ObjectiveCondition2D::evalF(const Real3 & p) {
  Real f_obj = 0.0;
  std::vector<gmds::Face> faces = m_node.get<gmds::Face>();

  for (int face_local_id = 0; face_local_id < faces.size(); face_local_id++) {
    gmds::Face  thisFace = faces[face_local_id];
    gmds::Node  adj1, adj2;
    thisFace.getAdjacentNodes(m_node, adj1, adj2);
    Real3 n1, n2, n12, n21;
    n1 = Real3(adj1.X(), adj1.Y(), adj1.Z());
    n2 = Real3(adj2.X(), adj2.Y(), adj2.Z());

    gmds::Node adj12, adj21, temp;
    thisFace.getAdjacentNodes(adj1, adj12, temp);
    if (temp.id() != m_node.id())
  	  adj12 = temp;
    thisFace.getAdjacentNodes(adj2, adj21, temp);
    if (temp.id() != m_node.id())
  	  adj21 = temp;

    n12 = Real3(adj12.X(), adj12.Y(), adj12.Z());
    n21 = Real3(adj21.X(), adj21.Y(), adj21.Z());

    f_obj += condition2D(p, n1, n2) + condition2D(n1, p, n12)
	  + condition2D(n2, p, n21);

  }
  return f_obj;
}

void ObjectiveCondition2D::evalDF(const Real3 & p, Real3 & grad) {
  Real f_obj = 0.0;
  grad = Real3::null();
  std::vector<gmds::Face > faces = m_node.get<gmds::Face>();

  for (int face_local_id = 0; face_local_id < faces.size(); face_local_id++) {
    gmds::Face thisFace = faces[face_local_id];
    gmds::Node adj1, adj2;
    thisFace.getAdjacentNodes(m_node, adj1, adj2);
    Real3 n1, n2, n12, n21;
    n1 = Real3(adj1.X(), adj1.Y(), adj1.Z());
    n2 = Real3(adj2.X(), adj2.Y(), adj2.Z());

    gmds::Node adj12, adj21, temp;
    thisFace.getAdjacentNodes(adj1, adj12, temp);
    if (temp.id() != m_node.id())
  	  adj12 = temp;
    thisFace.getAdjacentNodes(adj2, adj21, temp);
    if (temp.id() != m_node.id())
  	  adj21 = temp;

    n12 = Real3(adj12.X(), adj12.Y(), adj12.Z());
    n21 = Real3(adj21.X(), adj21.Y(), adj21.Z());

    grad += gradCond2DV1(p, n1, n2);
    grad += gradCond2DV2(n1, p, n12) + gradCond2DV2(n2, p, n21);
  }
}

Real ObjectiveCondition2D::evalFDF(const Real3 & p, Real3 & grad) {
  Real f_obj = 0.0;
  grad = Real3::null();
  std::vector<gmds::Face> faces = m_node.get<gmds::Face>();

  for (int face_local_id = 0; face_local_id < faces.size(); face_local_id++) {
    gmds::Face  thisFace = faces[face_local_id];
    gmds::Node adj1, adj2;
    thisFace.getAdjacentNodes(m_node, adj1, adj2);
    Real3 n1, n2, n12, n21;
    n1 = Real3(adj1.X(), adj1.Y(), adj1.Z());
    n2 = Real3(adj2.X(), adj2.Y(), adj2.Z());

    gmds::Node adj12, adj21, temp;
    thisFace.getAdjacentNodes(adj1, adj12, temp);
    if (temp.id() != m_node.id())
  	  adj12 = temp;
    thisFace.getAdjacentNodes(adj2, adj21, temp);
    if (temp.id() != m_node.id())
  	  adj21 = temp;

    n12 = Real3(adj12.X(), adj12.Y(), adj12.Z());
    n21 = Real3(adj21.X(), adj21.Y(), adj21.Z());

    f_obj += condition2D(p, n1, n2) + condition2D(n1, p, n12)
	  + condition2D(n2, p, n21);
    grad += gradCond2DV1(p, n1, n2);
    grad += gradCond2DV2(n1, p, n12) + gradCond2DV2(n2, p, n21);
  }
  return f_obj;
}
}

