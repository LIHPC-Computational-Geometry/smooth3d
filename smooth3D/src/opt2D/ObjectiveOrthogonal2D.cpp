/*
 * ObjectiveOrthogonal2D.cpp
 *
 *  Created on: 6 juil. 2015
 *      Author: weilljc
 */

#include "ObjectiveOrthogonal2D.h"
#include "math/orthogonal.h"
#include <cmath>

namespace Smooth3D {
inline double compute_cos_target(int number) {
  if (number == 4)
    return 0.0;
  else if (number == 3)
    return -0.5;
  else
    return std::cos(2 * M_PI / number);
}
ObjectiveOrthogonal2D::ObjectiveOrthogonal2D() {
  m_node = gmds::Node();
}

ObjectiveOrthogonal2D::ObjectiveOrthogonal2D(gmds::Node node) :
    m_node(node) {
}

ObjectiveOrthogonal2D::~ObjectiveOrthogonal2D() {
}

void ObjectiveOrthogonal2D::setNode(gmds::Node& node) {
  m_node = node;
}

Real ObjectiveOrthogonal2D::evalF(const Real3& p) {
  Real f_obj = 0.0;

  std::vector<gmds::Face> faces = m_node.get<gmds::Face>();

  if (faces.size() != 0) {
    double cos_target = compute_cos_target(faces.size());
    for (int face_local_id = 0; face_local_id < faces.size(); face_local_id++) {
      gmds::Face this_face = faces[face_local_id];
      gmds::Node adj1, adj2;
      this_face.getAdjacentNodes(m_node, adj1, adj2);
      Real3 n1, n2, n12, n21;
      n1 = Real3(adj1.X(), adj1.Y(), adj1.Z());
      n2 = Real3(adj2.X(), adj2.Y(), adj2.Z());

      f_obj += orthogonal2D(p, n1, n2, cos_target);

    }

    f_obj = f_obj / faces.size();
  }
  return f_obj;

}

void ObjectiveOrthogonal2D::evalDF(const Real3& p, Real3& grad) {
  grad = 0.0;
  std::vector<gmds::Face> faces = m_node.get<gmds::Face>();
  if (faces.size() != 0) {
    double cos_target = compute_cos_target(faces.size());

    for (int face_local_id = 0; face_local_id < faces.size(); face_local_id++) {
      gmds::Face this_face = faces[face_local_id];
      gmds::Node adj1, adj2;
      this_face.getAdjacentNodes(m_node, adj1, adj2);
      Real3 n1, n2, n12, n21;
      n1 = Real3(adj1.X(), adj1.Y(), adj1.Z());
      n2 = Real3(adj2.X(), adj2.Y(), adj2.Z());
      grad += gradOrtho2DV1(p, n1, n2, cos_target);
    }
    grad = grad * (1.0 / faces.size());
  }


}

Real ObjectiveOrthogonal2D::evalFDF(const Real3& p, Real3& grad) {
  Real f_obj = 0.0;
  grad = 0.0;
  std::vector<gmds::Face> faces = m_node.get<gmds::Face>();
  if (faces.size() != 0) {
    double cos_target = compute_cos_target(faces.size());
    for (int face_local_id = 0; face_local_id < faces.size(); face_local_id++) {
      gmds::Face this_face = faces[face_local_id];
      gmds::Node adj1, adj2;
      this_face.getAdjacentNodes(m_node, adj1, adj2);
      Real3 n1, n2, n12, n21;
      n1 = Real3(adj1.X(), adj1.Y(), adj1.Z());
      n2 = Real3(adj2.X(), adj2.Y(), adj2.Z());

      grad += gradOrtho2DV1(p, n1, n2, cos_target);

      f_obj += orthogonal2D(p, n1, n2, cos_target);
    }
    grad = grad * (1.0 / faces.size());
    f_obj /= faces.size();
  }
  return f_obj;
}

} /* namespace Smooth3D */
