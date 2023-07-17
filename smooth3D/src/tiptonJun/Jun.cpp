/*
 * Jun.c
 *
 *
 * methode de Tipton/Jun en 3D
 *
 *  Created on: 20 mars 2013
 *      Author: weilljc
 */

#include "smooth3D/smooth.h"
#include "math/Real3.h"
#include <vector>

namespace Smooth3D {

static const Real kTHexa[8][8] = { { 1, 1, 1, 1, 1, 1, 1, 1 }, { -1, 1, 1, -1, -1, 1,
    1, -1 }, { -1, -1, 1, 1, -1, -1, 1, 1 }, { -1, -1, -1, -1, 1, 1, 1, 1 }, {
    1, -1, 1, -1, 1, -1, 1, -1 }, { 1, 1, -1, -1, -1, -1, 1, 1 }, { 1, -1, -1,
    1, -1, 1, 1, -1 }, { -1, 1, -1, 1, 1, -1, 1, -1 } };

static inline void computeAijT3D8(Real alpha, Real beta, Real gamma,
    const Real3 U[8], Real a[3][3]) {
  Real3 df[3];

  df[0] = U[1] * 0.5 + U[4] * beta + U[6] * gamma + 2.0 * U[7] * beta * gamma;
  df[1] = U[2] * 0.5 + U[4] * alpha + U[5] * gamma + 2.0 * U[7] * alpha * gamma;
  df[2] = U[3] * 0.5 + U[5] * beta + U[6] * alpha + 2.0 * U[7] * beta * alpha;

  Real g[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = i; j < 3; ++j)
      g[i][j] = g[j][i] = math::scaMul(df[i], df[j]);

  a[0][0] = (g[1][1] * g[2][2] - g[2][1] * g[2][1]);
  a[1][1] = (g[0][0] * g[2][2] - g[2][0] * g[2][0]);
  a[2][2] = (g[0][0] * g[1][1] - g[1][0] * g[1][0]);
  a[0][1] = a[1][0] = (g[2][1] * g[2][0] - g[0][1] * g[2][2]);
  a[1][2] = a[2][1] = -(g[0][0] * g[2][1] - g[2][0] * g[0][1]);
  a[0][2] = a[2][0] = (g[0][1] * g[1][2] - g[1][1] * g[0][2]);
}

/*
 * Calcul du Jacobien et de la matrice aij pour une pyramide
 * La fonction de forme de reference est :
 * f(a,b,c) = = (1 - a - b - c + a*b + Min[a, b]*c)*
 *    n1 + (a - a*b - Min[a, b]*c)n2 + (a*b + Min[a, b]*c)n3 + (b - a b -
 *        Min[a, b]*c) n4  + c n5)
 *
 * Le changement de variables associe est
 * u1 = n1
 * u2 = -n1 + n2
 * u3 = -n1 + n4
 * u4 = -n1 + n5
 * u5 = n1 - n2 + n3 - n4
 *
 * On trouve alors
 * f(a,b,c) = u1 + a u2 + b u3 + c u4 + ( a b + c * Min[a,b]) u5
 *
 * Dans la fonction suivante a == alpha
 *                           b == beta
 *                           c == gamma
 */

static inline void computeAijPyramide(Real alpha, Real beta, Real gamma,
    const Real3 U[5], Real a[3][3]) {
  Real3 df[3];
  // Calcul des derivees partielles de F / alpha,beta,gamma
  if (alpha < beta) {
    df[0] = U[1] + (beta + gamma) * U[4];
    df[1] = U[2] + alpha * U[4];
    df[2] = U[3] + alpha * U[4];
  } else {
    df[0] = U[1] + beta * U[4];
    df[1] = U[2] + (alpha + gamma) * U[4];
    df[2] = U[3] + beta * U[4];
  }

  // Calcul de la matrice G et de son inverse
  Real g[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = i; j < 3; ++j)
      g[i][j] = g[j][i] = math::scaMul(df[i], df[j]);

  a[0][0] = (g[1][1] * g[2][2] - g[2][1] * g[2][1]);
  a[1][1] = (g[0][0] * g[2][2] - g[2][0] * g[2][0]);
  a[2][2] = (g[0][0] * g[1][1] - g[1][0] * g[1][0]);
  a[0][1] = a[1][0] = (g[2][1] * g[2][0] - g[0][1] * g[2][2]);
  a[1][2] = a[2][1] = -(g[0][0] * g[2][1] - g[2][0] * g[0][1]);
  a[0][2] = a[2][0] = (g[0][1] * g[1][2] - g[1][1] * g[0][2]);

}

/*
 * Calcul du Jacobien et de la matrice aij pour une prisme (pentaedre)
 * La fonction de forme de reference est :
 * f(a,b,c) = = c (a n1 + b n2 + (1 - b - a) n0)
 *                + (1-c) (a n4 + b n5 + (1 - a - b ) n3)
 *
 *
 * Le changement de variables associe est
 * u1 = n4
 * u2 = n5 - n4
 * u3 = n6 - n4
 * u4 = n1 - n4
 * u5 = -n1 + n2 + n4 - n5
 * u6 = -n1 + n3 + n4 - n6
 *
 * On trouve alors
 * f(a,b,c) = u1 + a u2 + b u3 + c u4 + a c u5 + bc u6
 *
 * Dans la fonction suivante a == alpha
 *                           b == beta
 *                           c == gamma
 */

static /*inline*/ void computeAijPenta(Real alpha, Real beta, Real gamma,
    const Real3 U[6], Real a[3][3]) {
  Real3 df[3];

  df[0] = U[1] + gamma * U[4];
  df[1] = U[2] + gamma * U[5];
  df[2] = U[3] + alpha * U[4] + beta * U[5];

  // Calcul de la matrice G et de son inverse
  Real g[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = i; j < 3; ++j)
      g[i][j] = g[j][i] = math::scaMul(df[i], df[j]);

  a[0][0] = (g[1][1] * g[2][2] - g[2][1] * g[2][1]);
  a[1][1] = (g[0][0] * g[2][2] - g[2][0] * g[2][0]);
  a[2][2] = (g[0][0] * g[1][1] - g[1][0] * g[1][0]);
  a[0][1] = a[1][0] = (g[2][1] * g[2][0] - g[0][1] * g[2][2]);
  a[1][2] = a[2][1] = -(g[0][0] * g[2][1] - g[2][0] * g[0][1]);
  a[0][2] = a[2][0] = (g[0][1] * g[1][2] - g[1][1] * g[0][2]);

}

/*!
 * \brief Lissage 3D par algo de Jun
 * Une iteration de Jun en 3D
 */

extern "C" int S3_Jun(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double *x, double *y, double *z, const double *weights, const double *relax,
    int_type n_iter) {
  const Real3 kNullVector = Real3(0., 0., 0.);
  std::vector<Real3> nouvelle_pos(nb_nodes);
  std::vector<Real> somme_poids(nb_nodes);
  int_type vtx, cell;
  const int_type *f_node_cell;
  for (int iter = 0; iter < n_iter; ++iter) {

    for (vtx = 0; vtx < nb_nodes; ++vtx) {
      nouvelle_pos[vtx] = kNullVector;
      somme_poids[vtx] = 0.0;
    }

    for (f_node_cell = nodes_number, cell = 0; cell < nb_cells; ++cell) {

      s3_cell_3d_t my_cell;
      s3_edge_t my_edges[12];
      int nb_edges, edge;
      my_cell.nb_nodes = static_cast<s3_geom_3d_t>(nb_nodes_per_cell[cell]);
      for (int vtx = 0; vtx < nb_nodes_per_cell[cell]; ++vtx)
	my_cell.vtx[vtx] = *(f_node_cell++);

      // La maille participe au lissage
      switch (my_cell.nb_nodes) {

      case 4: {

	Real3 u[4];

	/* U = T Nodes */

	u[0] = Real3(x[my_cell.vtx[0]], y[my_cell.vtx[0]], z[my_cell.vtx[0]]);
	u[1] = Real3(x[my_cell.vtx[1]], y[my_cell.vtx[1]], z[my_cell.vtx[1]])
	    - u[0];
	u[2] = Real3(x[my_cell.vtx[2]], y[my_cell.vtx[2]], z[my_cell.vtx[2]])
	    - u[0];
	u[3] = Real3(x[my_cell.vtx[3]], y[my_cell.vtx[3]], z[my_cell.vtx[3]])
	    - u[0];

	Real g[3][3];
	for (int i = 0; i < 3; ++i)
	  for (int j = i; j < 3; ++j)
	    g[i][j] = g[j][i] = math::scaMul(u[i + 1], u[j + 1]);

	Real a11, a12, a23, a22, a13, a33;
	a11 = (g[1][1] * g[2][2] - g[2][1] * g[2][1]);
	a22 = (g[0][0] * g[2][2] - g[2][0] * g[2][0]);
	a33 = (g[0][0] * g[1][1] - g[1][0] * g[1][0]);
	a12 = (g[2][1] * g[2][0] - g[0][1] * g[2][2]);
	a23 = -(g[0][0] * g[2][1] - g[2][0] * g[0][1]);
	a13 = (g[0][1] * g[1][2] - g[1][1] * g[0][2]);

	Real matrix[4][4];
	matrix[0][0] = a11 + a22 + a23 + 2 * a12 + 2 * a13 + 2 * a23;
	matrix[0][1] = -a11 - a12 - a13;
	matrix[0][2] = -a12 - a22 - a23;
	matrix[0][3] = -a13 - a23 - a33;

	matrix[1][0] = -a11 - a12 - a13;
	matrix[1][1] = a11;
	matrix[1][2] = a12;
	matrix[1][3] = a23;

	matrix[2][0] = -a12 - a22 - a23;
	matrix[2][1] = a12;
	matrix[2][2] = a22;
	matrix[2][3] = a23;

	matrix[3][0] = -a13 - a23 - a33;
	matrix[3][1] = a13;
	matrix[3][2] = a23;
	matrix[3][3] = a33;
	Real w = weights[cell] * (1.0 / 6.0);

	for (int i = 0; i < 4; ++i) {
	  int id1 = my_cell.vtx[i];
	  somme_poids[id1] -= w * matrix[i][i];
	  for (int j = i + 1; j < 4; ++j) {
	    int id2 = my_cell.vtx[j];
	    nouvelle_pos[id1] += w * matrix[i][j] * Real3(x[id2], y[id2], z[id2]);
	    nouvelle_pos[id2] += w * matrix[j][i] * Real3(x[id1], y[id1], z[id1]);
	  }
	}
      }
	break;
      case 5: {
	Real3 u[5];

	// Changement de variable
	u[0] = Real3(x[my_cell.vtx[0]], y[my_cell.vtx[0]], z[my_cell.vtx[0]]);
	u[1] = Real3(x[my_cell.vtx[1]], y[my_cell.vtx[1]], z[my_cell.vtx[1]])
	    - u[0];
	u[2] = Real3(x[my_cell.vtx[2]], y[my_cell.vtx[2]], z[my_cell.vtx[2]])
	    - u[0];
	u[3] = Real3(x[my_cell.vtx[3]], y[my_cell.vtx[3]], z[my_cell.vtx[3]])
	    - u[0];
	u[4] = u[1] + u[2] - u[3];

	// Calcul de la matrice de rigidite
	Real matrix[5][5];
	Real a[3][3];
	for (int i = 0; i < 5; ++i)
	  for (int j = 0; j < 5; ++j)
	    matrix[i][j] = 0.0;

	// Calcul des coefficient AIJ pour le point 1/4,1/2,1/4
	computeAijPyramide(0.25, 0.5, 0.25, u, a);
	Real a11, a12, a23, a22, a13, a33;
	a11 = a[0][0];
	a12 = a[0][1];
	a13 = a[0][2];
	a22 = a[1][1];
	a23 = a[1][2];
	a33 = a[2][2];

	matrix[0][0] = a11 / 16. + (3 * a12) / 8. + (3 * a13) / 8. + (9 * a22) / 16.
	    + (9 * a23) / 8. + (9 * a33) / 16.;
	matrix[0][1] = -a11 / 16. - a12 / 8. - a13 / 8. + (3 * a22) / 16.
	    + (3 * a23) / 8. + (3 * a33) / 16.;
	matrix[0][2] = (-3 * a11) / 16. - (5 * a12) / 8. - (5 * a13) / 8.
	    - (3 * a22) / 16. - (3 * a23) / 8. - (3 * a33) / 16.;
	matrix[0][3] = (3 * a11) / 16. + (3 * a12) / 8. + (5 * a13) / 8.
	    - (9 * a22) / 16. - (3 * a23) / 8. + (3 * a33) / 16.;
	matrix[0][4] = -a13 / 4. - (3 * a23) / 4. - (3 * a33) / 4.;

	matrix[1][0] = matrix[0][1];
	matrix[1][1] = a11 / 16. - a12 / 8. - a13 / 8. + a22 / 16. + a23 / 8.
	    + a33 / 16.;
	matrix[1][2] = (3 * a11) / 16. - a12 / 8. - a13 / 8. - a22 / 16. - a23 / 8.
	    - a33 / 16.;
	matrix[1][3] = (-3 * a11) / 16. + (3 * a12) / 8. + a13 / 8. - (3 * a22) / 16.
	    - a23 / 8. + a33 / 16.;
	matrix[1][4] = a13 / 4. - a23 / 4. - a33 / 4.;

	matrix[2][0] = matrix[0][2];
	matrix[2][1] = matrix[1][2];
	matrix[2][2] = (9 * a11) / 16. + (3 * a12) / 8. + (3 * a13) / 8. + a22 / 16.
	    + a23 / 8. + a33 / 16.;
	matrix[2][3] = (-9 * a11) / 16. + (3 * a12) / 8. - (3 * a13) / 8.
	    + (3 * a22) / 16. + a23 / 8. - a33 / 16.;
	matrix[2][4] = (3 * a13) / 4. + a23 / 4. + a33 / 4.;

	matrix[3][0] = matrix[0][3];
	matrix[3][1] = matrix[1][3];
	matrix[3][2] = matrix[2][3];
	matrix[3][3] = (9 * a11) / 16. - (9 * a12) / 8. + (3 * a13) / 8.
	    + (9 * a22) / 16. - (3 * a23) / 8. + a33 / 16.;
	matrix[3][4] = (-3 * a13) / 4. + (3 * a23) / 4. - a33 / 4.;

	matrix[4][0] = matrix[0][4];
	matrix[4][1] = matrix[1][4];
	matrix[4][2] = matrix[2][4];
	matrix[4][3] = matrix[3][4];
	matrix[4][4] = a33;

	//  Calcul des coefficient AIJ pour le point 1/2,1/4,1/4
	computeAijPyramide(0.5, 0.25, 0.25, u, a);
	a11 = a[0][0];
	a12 = a[0][1];
	a13 = a[0][2];
	a22 = a[1][1];
	a23 = a[1][2];
	a33 = a[2][2];

	matrix[0][0] += (9 * a11) / 16. + (3 * a12) / 8. + (9 * a13) / 8. + a22 / 16.
	    + (3 * a23) / 8. + (9 * a33) / 16.;
	matrix[0][1] += (-9 * a11) / 16. + (3 * a12) / 8. - (3 * a13) / 8.
	    + (3 * a22) / 16. + (5 * a23) / 8. + (3 * a33) / 16.;
	matrix[0][2] += (-3 * a11) / 16. - (5 * a12) / 8. - (3 * a13) / 8.
	    - (3 * a22) / 16. - (5 * a23) / 8. - (3 * a33) / 16.;
	matrix[0][3] += (3 * a11) / 16. - a12 / 8. + (3 * a13) / 8. - a22 / 16.
	    - a23 / 8. + (3 * a33) / 16.;
	matrix[0][4] += (-3 * a13) / 4. - a23 / 4. - (3 * a33) / 4.;

	matrix[1][0] += (-9 * a11) / 16. + (3 * a12) / 8. - (3 * a13) / 8.
	    + (3 * a22) / 16. + (5 * a23) / 8. + (3 * a33) / 16.;
	matrix[1][1] += (9 * a11) / 16. - (9 * a12) / 8. - (3 * a13) / 8.
	    + (9 * a22) / 16. + (3 * a23) / 8. + a33 / 16.;
	matrix[1][2] += (3 * a11) / 16. + (3 * a12) / 8. + a13 / 8. - (9 * a22) / 16.
	    - (3 * a23) / 8. - a33 / 16.;
	matrix[1][3] += (-3 * a11) / 16. + (3 * a12) / 8. - a13 / 8.
	    - (3 * a22) / 16. + a23 / 8. + a33 / 16.;
	matrix[1][4] += (3 * a13) / 4. - (3 * a23) / 4. - a33 / 4.;

	matrix[2][0] += (-3 * a11) / 16. - (5 * a12) / 8. - (3 * a13) / 8.
	    - (3 * a22) / 16. - (5 * a23) / 8. - (3 * a33) / 16.;
	matrix[2][1] += (3 * a11) / 16. + (3 * a12) / 8. + a13 / 8. - (9 * a22) / 16.
	    - (3 * a23) / 8. - a33 / 16.;
	matrix[2][2] += a11 / 16. + (3 * a12) / 8. + a13 / 8. + (9 * a22) / 16.
	    + (3 * a23) / 8. + a33 / 16.;
	matrix[2][3] += -a11 / 16. - a12 / 8. - a13 / 8. + (3 * a22) / 16. - a23 / 8.
	    - a33 / 16.;
	matrix[2][4] += a13 / 4. + (3 * a23) / 4. + a33 / 4.;

	matrix[3][0] += (3 * a11) / 16. - a12 / 8. + (3 * a13) / 8. - a22 / 16.
	    - a23 / 8. + (3 * a33) / 16.;
	matrix[3][1] += (-3 * a11) / 16. + (3 * a12) / 8. - a13 / 8.
	    - (3 * a22) / 16. + a23 / 8. + a33 / 16.;
	matrix[3][2] += -a11 / 16. - a12 / 8. - a13 / 8. + (3 * a22) / 16. - a23 / 8.
	    - a33 / 16.;
	matrix[3][3] += a11 / 16. - a12 / 8. + a13 / 8. + a22 / 16. - a23 / 8.
	    + a33 / 16.;
	matrix[3][4] += -a13 / 4. + a23 / 4. - a33 / 4;

	matrix[4][0] += (-3 * a13) / 4. - a23 / 4. - (3 * a33) / 4.;
	matrix[4][1] += (3 * a13) / 4. - (3 * a23) / 4. - a33 / 4.0;
	matrix[4][2] += a13 / 4. + (3 * a23) / 4. + a33 / 4.;
	matrix[4][3] += -a13 / 4. + a23 / 4. - a33 / 4;
	matrix[4][4] += a33;

	Real w = weights[cell] / 12.0; // 1/6 pour chaque point * 1/2 overall

	for (int i = 0; i < 5; ++i) {
	  int id1 = my_cell.vtx[i];
	  somme_poids[id1] -= w * matrix[i][i];
	  for (int j = i + 1; j < 5; ++j) {
	    int id2 = my_cell.vtx[j];
	    nouvelle_pos[id1] += w * matrix[i][j] * Real3(x[id2], y[id2], z[id2]);
	    nouvelle_pos[id2] += w * matrix[j][i] * Real3(x[id1], y[id1], z[id1]);
	  }
	}

      }
	break;

      case 6: {
	Real3 u[6];

	/* U = T Nodes */

	u[0] = Real3(x[my_cell.vtx[3]], y[my_cell.vtx[3]], z[my_cell.vtx[3]]);
	u[1] = Real3(x[my_cell.vtx[4]], y[my_cell.vtx[4]], z[my_cell.vtx[4]])
	    - u[0];
	u[2] = Real3(x[my_cell.vtx[5]], y[my_cell.vtx[5]], z[my_cell.vtx[5]])
	    - u[0];
	u[3] = Real3(x[my_cell.vtx[0]], y[my_cell.vtx[0]], z[my_cell.vtx[0]])
	    - u[0];
	u[4] = -u[0] - u[1]
	    + Real3(x[my_cell.vtx[1]], y[my_cell.vtx[1]], z[my_cell.vtx[1]]);
	u[5] = Real3(x[my_cell.vtx[2]], y[my_cell.vtx[2]], z[my_cell.vtx[2]])
	    - u[2] - u[0];

	Real matrix[6][6];
	Real a[3][3];

	for (int i = 0; i < 6; ++i)
	  for (int j = 0; j < 6; ++j)
	    matrix[i][j] = 0.0;

	// Calcul des coefficient AIJ pour le point 1/3,1/3,1/2
	computeAijPenta(1.0 / 3.0, 1.0 / 3.0, 0.5, u, a);

	Real a11, a12, a23, a22, a13, a33;
	a11 = a[0][0];
	a12 = a[0][1];
	a13 = a[0][2];
	a22 = a[1][1];
	a23 = a[1][2];
	a33 = a[2][2];

	matrix[0][0] = a11 * 0.25 + a12 * 0.5 - a13 / 3. + a22 * 0.25 - a23 / 3.
	    + a33 / 9.;
	matrix[0][1] = -a11 * 0.25 - a12 * 0.25 - a23 / 6. + a33 / 9.;
	matrix[0][2] = -a12 * 0.25 - a13 / 6. - a22 * 0.25 + a33 / 9.;
	matrix[0][3] = a11 * 0.25 + a12 * 0.5 + a22 * 0.25 - a33 / 9.;
	matrix[0][4] = -a11 * 0.25 - a12 * 0.25 + a13 / 3. + a23 / 6. - a33 / 9.;
	matrix[0][5] = -a12 * 0.25 + a13 / 6. - a22 * 0.25 + a23 / 3. - a33 / 9.;

	matrix[1][0] = matrix[0][1];
	matrix[1][1] = a11 * 0.25 + a13 / 3. + a33 / 9.;
	matrix[1][2] = a12 * 0.25 + a13 / 6. + a23 / 6. + a33 / 9.;
	matrix[1][3] = -a11 * 0.25 - a12 * 0.25 - a13 / 3. - a23 / 6. - a33 / 9.;
	matrix[1][4] = a11 * 0.25 - a33 / 9.;
	matrix[1][5] = a12 * 0.25 - a13 / 6. + a23 / 6. - a33 / 9.;

	matrix[2][0] = matrix[0][2];
	matrix[2][1] = matrix[1][2];
	matrix[2][2] = a22 * 0.25 + a23 / 3. + a33 / 9.;
	matrix[2][3] = -a12 * 0.25 - a13 / 6. - a22 * 0.25 - a23 / 3. - a33 / 9.;
	matrix[2][4] = a12 * 0.25 + a13 / 6. - a23 / 6. - a33 / 9.;
	matrix[2][5] = a22 * 0.25 - a33 / 9.;

	matrix[3][0] = matrix[0][3];
	matrix[3][1] = matrix[1][3];
	matrix[3][2] = matrix[2][3];
	matrix[3][3] = a11 * 0.25 + a12 * 0.5 + a13 / 3. + a22 * 0.25 + a23 / 3.
	    + a33 / 9.;
	matrix[3][4] = -a11 * 0.25 - a12 * 0.25 + a23 / 6. + a33 / 9.;
	matrix[3][5] = -a12 * 0.25 + a13 / 6. - a22 * 0.25 + a33 / 9.;

	matrix[4][0] = matrix[0][4];
	matrix[4][1] = matrix[1][4];
	matrix[4][2] = matrix[2][4];
	matrix[4][3] = matrix[3][4];
	matrix[4][4] = a11 * 0.25 - a13 / 3. + a33 / 9.;
	matrix[4][5] = a12 * 0.25 - a13 / 6. - a23 / 6. + a33 / 9.;

	matrix[5][0] = matrix[0][5];
	matrix[5][1] = matrix[1][5];
	matrix[5][2] = matrix[2][5];
	matrix[5][3] = matrix[3][5];
	matrix[5][4] = matrix[4][5];
	matrix[5][5] = a22 * 0.25 - a23 / 3. + a33 / 9.0;

	Real w = 0.5 * weights[cell];
	for (int i = 0; i < 6; ++i) {
	  int id1 = my_cell.vtx[i];
	  somme_poids[id1] -= w * matrix[i][i];
	  for (int j = i + 1; j < 6; ++j) {
	    int id2 = my_cell.vtx[j];
	    nouvelle_pos[id1] += w * matrix[i][j] * Real3(x[id2], y[id2], z[id2]);
	    nouvelle_pos[id2] += w * matrix[j][i] * Real3(x[id1], y[id1], z[id1]);
	  }
	}

      }

	break;

      case 8: {

	Real3 u[8];

	// Calcul du produit U = kTHexa Nodes
	for (int i = 0; i < 8; ++i) {
	  Real3 t = kNullVector;

	  for (int j = 0; j < 8; ++j) {

	    t += Real3(x[my_cell.vtx[j]], y[my_cell.vtx[j]], z[my_cell.vtx[j]])
		* kTHexa[i][j];
	  }
	  u[i] = t;
	}

	typedef Real mat[3][3];
	mat b[9];

	computeAijT3D8(0, 0, 0, u, b[0]);

	computeAijT3D8(-0.5, -0.5, -0.5, u, b[1]);
	computeAijT3D8(0.5, -0.5, -0.5, u, b[2]);
	computeAijT3D8(0.5, 0.5, -0.5, u, b[3]);
	computeAijT3D8(-0.5, 0.5, -0.5, u, b[4]);

	computeAijT3D8(-0.5, -0.5, 0.5, u, b[5]);
	computeAijT3D8(0.5, -0.5, 0.5, u, b[6]);
	computeAijT3D8(0.5, 0.5, 0.5, u, b[7]);
	computeAijT3D8(-0.5, 0.5, 0.5, u, b[8]);

	Real a[3][3];
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
	    a[i][j] = b[0][i][j];
	// Calcul de la matrice de rigidite
	Real matrix[8][8];

	for (int i = 0; i < 3; i++)
	  a[i][i] = b[0][i][i] + b[1][i][i];

	matrix[0][0] = a[0][0] + a[1][1] + a[2][2] + a[0][1] + a[0][2] + a[1][2];
	matrix[0][1] = -a[0][0];
	matrix[0][2] = -a[0][1];
	matrix[0][3] = -a[1][1];
	matrix[0][4] = -a[2][2];
	matrix[0][5] = -a[0][2];
	matrix[0][6] = 0;
	matrix[0][7] = -a[1][2];

	for (int i = 0; i < 3; i++)
	  a[i][i] = b[0][i][i] + b[2][i][i];

	matrix[1][0] = -a[0][0];
	matrix[1][1] = a[0][0] + a[1][1] + a[2][2] - a[0][1] - a[0][2] + a[1][2];
	matrix[1][2] = -a[1][1];
	matrix[1][3] = a[0][1];
	matrix[1][4] = a[0][2];
	matrix[1][5] = -a[2][2];
	matrix[1][6] = -a[1][2];
	matrix[1][7] = 0;

	for (int i = 0; i < 3; i++)
	  a[i][i] = b[0][i][i] + b[3][i][i];

	matrix[2][0] = -a[0][1];
	matrix[2][1] = -a[1][1];
	matrix[2][2] = a[0][0] + a[1][1] + a[2][2] + a[0][1] - a[0][2] - a[1][2];
	matrix[2][3] = -a[0][0];
	matrix[2][4] = 0;
	matrix[2][5] = a[1][2];
	matrix[2][6] = -a[2][2];
	matrix[2][7] = a[0][2];

	for (int i = 0; i < 3; i++)
	  a[i][i] = b[0][i][i] + b[4][i][i];

	matrix[3][0] = -a[1][1];
	matrix[3][1] = a[0][1];
	matrix[3][2] = -a[0][0];
	matrix[3][3] = a[0][0] + a[1][1] + a[2][2] - a[0][1] + a[0][2] - a[1][2];
	matrix[3][4] = a[1][2];
	matrix[3][5] = 0;
	matrix[3][6] = -a[0][2];
	matrix[3][7] = -a[2][2];

	for (int i = 0; i < 3; i++)
	  a[i][i] = b[0][i][i] + b[5][i][i];

	matrix[4][0] = -a[2][2];
	matrix[4][1] = a[0][2];
	matrix[4][2] = 0.0;
	matrix[4][3] = a[1][2];
	matrix[4][4] = a[0][0] + a[1][1] + a[2][2] + a[0][1] - a[0][2] - a[1][2];
	matrix[4][5] = -a[0][0];
	matrix[4][6] = -a[0][1];
	matrix[4][7] = -a[1][1];

	for (int i = 0; i < 3; i++)
	  a[i][i] = b[0][i][i] + b[6][i][i];

	matrix[5][0] = -a[0][2];
	matrix[5][1] = -a[2][2];
	matrix[5][2] = a[1][2];
	matrix[5][3] = 0;
	matrix[5][4] = -a[0][0];
	matrix[5][5] = a[0][0] + a[1][1] + a[2][2] - a[0][1] + a[0][2] - a[1][2];
	matrix[5][6] = -a[1][1];
	matrix[5][7] = a[0][1];

	for (int i = 0; i < 3; i++)
	  a[i][i] = b[0][i][i] + b[7][i][i];

	matrix[6][0] = 0;
	matrix[6][1] = -a[1][2];
	matrix[6][2] = -a[2][2];
	matrix[6][3] = -a[0][2];
	matrix[6][4] = -a[0][1];
	matrix[6][5] = -a[1][1];
	matrix[6][6] = a[0][0] + a[1][1] + a[2][2] + a[0][1] + a[0][2] + a[1][2];
	matrix[6][7] = -a[0][0];

	for (int i = 0; i < 3; i++)
	  a[i][i] = b[0][i][i] + b[8][i][i];

	matrix[7][0] = -a[1][2];
	matrix[7][1] = 0;
	matrix[7][2] = a[0][2];
	matrix[7][3] = -a[2][2];
	matrix[7][4] = -a[1][1];
	matrix[7][5] = a[0][1];
	matrix[7][6] = -a[0][0];
	matrix[7][7] = a[0][0] + a[1][1] + a[2][2] - a[0][1] - a[0][2] + a[1][2];

	Real w = weights[cell];
	//	  	  pinfo() << " poids dans Jun " << 0.25*w ;
	for (int i = 0; i < 8; ++i) {
	  int id1 = my_cell.vtx[i];
	  somme_poids[id1] -= w * matrix[i][i];
	  for (int j = i + 1; j < 8; ++j) {
	    int id2 = my_cell.vtx[j];
	    nouvelle_pos[id1] += w * matrix[i][j] * Real3(x[id2], y[id2], z[id2]);
	    nouvelle_pos[id2] += w * matrix[j][i] * Real3(x[id1], y[id1], z[id1]);
	  }
	}

      }
	break;
      default:
	break;
      }
    }
    for (vtx = 0; vtx < nb_nodes; ++vtx) {
      double rel = relax[vtx];
      int seen = 0;
      if (rel != 0.0) {
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


  } /* End for (iter) */

  return (1);
}


} /* End namespace */


