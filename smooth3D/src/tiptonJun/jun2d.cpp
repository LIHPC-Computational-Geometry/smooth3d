/*
 * jun2d.cpp
 *
 *  Created on: 18 nov. 2015
 *      Author: weilljc
 */

#include "smooth3D/smooth.h"
#include "math/Real3.h"
#include <vector>

namespace Smooth3D {
static inline
void computeAijQuad(Real alpha, Real beta, const Real3 & u1, const Real3 & u2,
    const Real3 & u3, const Real3 & u4, Real & a11, Real & a12, Real & a21,
    Real & a22) {
  Real3 df_a = u2 + 2 * beta * u4;
  Real3 df_b = u3 + 2 * alpha * u4;

  a11 = math::scaMul(df_b, df_b);
  a22 = math::scaMul(df_a, df_a);
  a12 = a21 = -math::scaMul(df_a, df_b);
}

/*!
 * \brief Lissage 2D par Algo de Tipton sur une surface
 */
extern "C" int S3_Jun2D(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double *x, double *y, double *z, const double *weights, const double *relax,
    int_type n_iter) {

  std::vector<Real3> nouvelle_pos(nb_nodes);
  std::vector<Real> somme_poids(nb_nodes);
  std::vector<bool> non_valid(nb_nodes);
  const Real3 kNullVector = Real3(0., 0., 0.);

  for (int iter_n = 0; iter_n < n_iter; ++iter_n) {

    for (int_type vtx = 0; vtx < nb_nodes; ++vtx) {
      nouvelle_pos[vtx] = kNullVector;
      somme_poids[vtx] = 0.0;
      non_valid[vtx] = false;
    }
    int nb_face_processed = 0;

    int_type cell;
    const int_type * f_node_cell;
    for (f_node_cell = nodes_number, cell = 0; cell < nb_cells; ++cell) {

      s3_cell_3d_t my_cell;
      s3_edge_t my_edges[12];
      int nb_edges, edge;
      my_cell.nb_nodes = static_cast<s3_geom_3d_t>(nb_nodes_per_cell[cell]);
      for (int vtx = 0; vtx < nb_nodes_per_cell[cell]; ++vtx)
	my_cell.vtx[vtx] = *(f_node_cell++);

      switch (my_cell.nb_nodes) {
      case 3: {
	Real3 u[3];
	++nb_face_processed;
	u[0] = Real3(x[my_cell.vtx[0]], y[my_cell.vtx[0]], z[my_cell.vtx[0]]);
	u[1] = Real3(x[my_cell.vtx[1]], y[my_cell.vtx[1]], z[my_cell.vtx[1]]);

	u[2] = Real3(x[my_cell.vtx[2]], y[my_cell.vtx[2]], z[my_cell.vtx[2]]);

	Real3 u1 = u[1] - u[0];
	Real3 u2 = u[2] - u[0];

	Real a11, a12, a22;
	a11 = math::scaMul(u2, u2);
	a22 = math::scaMul(u1, u1);
	a12 = -math::scaMul(u1, u2);

	Real stiffness_matrix[3][3];
	stiffness_matrix[0][0] = a11 + 2.0 * a12 + a22;
	stiffness_matrix[0][1] = -a11 - a12;
	stiffness_matrix[0][2] = -a22 - a12;

	stiffness_matrix[1][0] = -a11 - a12;
	stiffness_matrix[1][1] = a11;
	stiffness_matrix[1][2] = a12;

	stiffness_matrix[2][0] = -a22 - a12;
	stiffness_matrix[2][1] = a12;
	stiffness_matrix[2][2] = a22;

	for (int i = 0; i < 3; ++i) {
	  int_type id1 = my_cell.vtx[i];
	  for (int j = i; j < 3; ++j) {
	    int_type id2 = my_cell.vtx[j];
	    if (i == j)
	      somme_poids[id1] -= stiffness_matrix[i][i];
	    else {
	      nouvelle_pos[id1] += stiffness_matrix[i][j]
		  * Real3(x[id2], y[id2], z[id2]);
	      nouvelle_pos[id2] += stiffness_matrix[j][i]
		  * Real3(x[id1], y[id1], z[id1]);
	    }
	  }
	}
      }
	break;
      case 4: {
	Real3 u[4];
	++nb_face_processed;
	for (int i = 0; i < 4; ++i)
	u[i] = Real3(x[my_cell.vtx[i]], y[my_cell.vtx[i]], z[my_cell.vtx[i]]);

	Real3 u1 = (u[0] + u[1] + u[2] + u[3]) * 0.5;
	Real3 u2 = (u[0] + u[1] - u[2] - u[3]) * 0.5;
	Real3 u3 = (-u[0] + u[1] + u[2] - u[3]) * 0.5;
	Real3 u4 = (-u[0] + u[1] - u[2] + u[3]) * 0.5;

	Real a11[5], a22[5], a12[5], a21[5];

	computeAijQuad(0.0, 0.0, u1, u2, u3, u4, a11[0], a12[0], a21[0],
	    a22[0]);

	computeAijQuad(0.5, -0.5, u1, u2, u3, u4, a11[1], a12[1], a21[1],
	    a22[1]);
	computeAijQuad(0.5, 0.5, u1, u2, u3, u4, a11[2], a12[2], a21[2],
	    a22[2]);
	computeAijQuad(-0.5, 0.5, u1, u2, u3, u4, a11[3], a12[3], a21[3],
	    a22[3]);
	computeAijQuad(-0.5, -0.5, u1, u2, u3, u4, a11[4], a12[4], a21[4],
	    a22[4]);

	Real stiffness_matrix[4][4];
	stiffness_matrix[0][0] = a11[0] - a12[0] + a22[0] + a11[1] + a22[1];
	stiffness_matrix[0][1] = -a22[0] - a22[1];
	stiffness_matrix[0][2] = a12[0];
	stiffness_matrix[0][3] = -a11[0] - a11[1];

	stiffness_matrix[1][0] = -a22[0] - a22[2];
	stiffness_matrix[1][1] = a11[0] + a12[0] + a22[0] + a11[2] + a22[2];
	stiffness_matrix[1][2] = -a11[0] - a11[2];
	stiffness_matrix[1][3] = -a12[0];

	stiffness_matrix[2][0] = a12[0];
	stiffness_matrix[2][1] = -a11[0] - a11[3];
	stiffness_matrix[2][2] = a11[0] - a12[0] + a22[0] + a11[3] + a22[3];
	stiffness_matrix[2][3] = -a22[0] - a22[3];

	stiffness_matrix[3][0] = -a11[0] - a11[4];
	stiffness_matrix[3][1] = -a12[0];
	stiffness_matrix[3][2] = -a22[0] - a22[4];
	stiffness_matrix[3][3] = a11[0] + a12[0] + a22[0] + a11[4] + a22[4];
	for (int i = 0; i < 4; ++i) {
	  int_type id1 = my_cell.vtx[i];
	  for (int j = i; j < 4; ++j) {
	    int_type id2 = my_cell.vtx[j];
	    if (i == j)
	      somme_poids[id1] -= stiffness_matrix[i][i];
	    else {
	      nouvelle_pos[id1] += stiffness_matrix[i][j]
		  * Real3(x[id2], y[id2], z[id2]);
	      nouvelle_pos[id2] += stiffness_matrix[j][i]
		  * Real3(x[id1], y[id1], z[id1]);
	    }
	  }
	}

      }

	break;

      default:
	for (int i = 0; i < my_cell.nb_nodes; ++i)
	  non_valid[my_cell.vtx[i]] = true;
	break;
      }
    }
    for (int_type vtx = 0; vtx < nb_nodes; ++vtx) {
      double rel = relax[vtx];
      int seen = 0;
      if (rel != 0.0 && non_valid[vtx] == false) {
	double m_x0, m_y0, m_z0, my_self;
	my_self = somme_poids[vtx];

	if (my_self != 0.0) {
	  Real3 calcul = nouvelle_pos[vtx] / my_self;

	  x[vtx] = (1 - rel) * x[vtx] + rel * calcul.x();
	  y[vtx] = (1 - rel) * y[vtx] + rel * calcul.y();
	  z[vtx] = (1 - rel) * z[vtx] + rel * calcul.z();
	}
      }
    }
  }
  return (1);
}
}

