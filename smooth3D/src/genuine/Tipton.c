/*
 * Methode de Winslow Crowley avec Poids
 * 
 * Article de base "Grid Optimization by Equipotential Relaxation"
 *    R.E. Tipton . LLNL Jult 15, 92
 *
 *
 * $Id$
 * $Log$
 * Revision 1.5  2005/07/06 13:30:17  weilljc
 * Adaptation TERA 10 et RH4
 *
 * Revision 1.4  2004/11/08 08:02:12  weilljc
 * Tipton toutes sortes de mailles version Jun
 *
 * Revision 1.3  2001/01/26 09:11:29  weilljc
 * debug
 *
 * Revision 1.2  2000/11/23 09:33:31  weilljc
 * Modification de Tipton...
 *    Petites optimisations Memoire + Temps
 *
 * Revision 1.1  2000/11/16 11:47:49  weilljc
 * Version 0.0.3
 *
 *
 *    Ajout de S3_Tipton...
 *
 */

#include "smooth3D/smooth.h"
#include <math.h>

#ifndef M_SQRT2
 #define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
 #define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#endif

struct point_s {
  double x;
  double y;
  double z;
};

typedef struct point_s point_t;

static const double kTHexa[8][8] = { { 1, 1, 1, 1, 1, 1, 1, 1 }, { 1, 1, -1, -1,
    1, 1, -1, -1 }, { -1, 1, 1, -1, -1, 1, 1, -1 },
    { -1, -1, -1, -1, 1, 1, 1, 1 }, { -1, 1, -1, 1, -1, 1, -1, 1 }, { 1, -1, -1,
        1, -1, 1, 1, -1 }, { -1, -1, 1, 1, 1, 1, -1, -1 }, { 1, -1, 1, -1, -1,
        1, -1, 1 } };

/*
 u[0] = m_ale_geom.m_node_coord[cell.node (3)];
 u[1] = m_ale_geom.m_node_coord[cell.node (4)] - u[0];
 u[2] = m_ale_geom.m_node_coord[cell.node (5)] - u[0];

 u[3] = m_ale_geom.m_node_coord[cell.node (0)] - u[0];
 u[4] =
 m_ale_geom.m_node_coord[cell.node (1)] -
 m_ale_geom.m_node_coord[cell.node (4)] - u[3];
 u[5] =
 m_ale_geom.m_node_coord[cell.node (2)] -
 m_ale_geom.m_node_coord[cell.node (5)] - u[3];
 */

static const double kTPenta[6][6] = { { 0, 0, 0, 1, 0, 0 },
    { 0, 0, 0, -1, 1, 0 }, { 0, 0, 0, -1, 0, 1 }, { 1, 0, 0, -1, 0, 0 }, { -1,
        1, 0, 1, -1, 0 }, { -1, 0, 1, 1, 0, -1 } };

inline void computeAijT3D8(double alpha, double beta, double gamma,
    const point_t U[8], double a[3][3]) {
  point_t df[3];
  double g[3][3];

  int i, j;

  df[0].x = U[1].x * 0.5 + U[4].x * beta + U[6].x * gamma
      + 2.0 * U[7].x * beta * gamma;
  df[0].y = U[1].y * 0.5 + U[4].y * beta + U[6].y * gamma
      + 2.0 * U[7].y * beta * gamma;
  df[0].z = U[1].z * 0.5 + U[4].z * beta + U[6].z * gamma
      + 2.0 * U[7].z * beta * gamma;

  df[1].x = U[2].x * 0.5 + U[4].x * alpha + U[5].x * gamma
      + 2.0 * U[7].x * alpha * gamma;
  df[1].y = U[2].y * 0.5 + U[4].y * alpha + U[5].y * gamma
      + 2.0 * U[7].y * alpha * gamma;
  df[1].z = U[2].z * 0.5 + U[4].z * alpha + U[5].z * gamma
      + 2.0 * U[7].z * alpha * gamma;

  df[2].x = U[3].x * 0.5 + U[5].x * beta + U[6].x * alpha
      + 2.0 * U[7].x * beta * alpha;
  df[2].y = U[3].y * 0.5 + U[5].y * beta + U[6].y * alpha
      + 2.0 * U[7].y * beta * alpha;
  df[2].z = U[3].z * 0.5 + U[5].z * beta + U[6].z * alpha
      + 2.0 * U[7].z * beta * alpha;

  for (i = 0; i < 3; ++i)
    for (j = i; j < 3; ++j)
      g[i][j] = g[j][i] = df[i].x * df[j].x + df[i].y * df[j].y
	  + df[i].z * df[i].z;

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
 * u2 =  (-n1 + n2) / rac(2)
 * u3 = (-n1 + n4) / rac(2)
 * u4 = ( -n1 + n5) / rac(2)
 * u5 =  (n1 - n2 + n3 - n4) / 2 
 *
 * On trouve alors
 * f(a,b,c) = u1 +  rac(2) a u2 + rac(2) b u3 + rac(2) c u4 +  2 ( a b + c * Min[a,b]) u5
 *
 * Dans la fonction suivante a == alpha
 *                           b == beta
 *                           c == gamma
 */

inline void computeAijPyramide(double alpha, double beta, double gamma,
    const point_t U[5], double a[3][3]) {
  point_t df[3];
  double g[3][3];
  int i, j;
  // Calcul des derivees partielles de F / alpha,beta,gamma
  if (alpha < beta) {
    df[0].x = M_SQRT2 * U[1].x + 2.0 * (beta + gamma) * U[4].x;
    df[0].y = M_SQRT2 * U[1].y + 2.0 * (beta + gamma) * U[4].y;
    df[0].z = M_SQRT2 * U[1].z + 2.0 * (beta + gamma) * U[4].z;

    df[1].x = M_SQRT2 * U[2].x + 2.0 * alpha * U[4].x;
    df[1].y = M_SQRT2 * U[2].y + 2.0 * alpha * U[4].y;
    df[1].z = M_SQRT2 * U[2].z + 2.0 * alpha * U[4].z;

    df[2].x = M_SQRT2 * U[3].x + 2.0 * alpha * U[4].x;
    df[2].y = M_SQRT2 * U[3].y + 2.0 * alpha * U[4].y;
    df[2].z = M_SQRT2 * U[3].z + 2.0 * alpha * U[4].z;
  } else {
    df[0].x = M_SQRT2 * U[1].x + 2.0 * beta * U[4].x;
    df[0].y = M_SQRT2 * U[1].y + 2.0 * beta * U[4].y;
    df[0].z = M_SQRT2 * U[1].z + 2.0 * beta * U[4].z;

    df[1].x = M_SQRT2 * U[2].x + 2.0 * (alpha + gamma) * U[4].x;
    df[1].y = M_SQRT2 * U[2].y + 2.0 * (alpha + gamma) * U[4].y;
    df[1].z = M_SQRT2 * U[2].z + 2.0 * (alpha + gamma) * U[4].z;

    df[2].x = M_SQRT2 * U[3].x + 2.0 * beta * U[4].x;
    df[2].y = M_SQRT2 * U[3].y + 2.0 * beta * U[4].y;
    df[2].z = M_SQRT2 * U[3].z + 2.0 * beta * U[4].z;
  }

  // Calcul de la matrice G et de son inverse
  for (i = 0; i < 3; ++i)
    for (j = i; j < 3; ++j)
      g[i][j] = g[j][i] = df[i].x * df[j].x + df[i].y * df[j].y
	  + df[i].z * df[i].z;

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

static void computeAijPenta(double alpha, double beta, double gamma,
    const point_t U[6], double a[3][3]) {
  point_t df[3];
  double g[3][3];
  int i, j;
  df[0].x = U[1].x + gamma * U[4].x;
  df[0].y = U[1].y + gamma * U[4].y;
  df[0].z = U[1].z + gamma * U[4].z;

  df[1].x = U[2].x + gamma * U[5].x;
  df[1].y = U[2].y + gamma * U[5].y;
  df[1].z = U[2].z + gamma * U[5].z;

  df[2].x = U[3].x + alpha * U[4].x + beta * U[5].x;
  df[2].y = U[3].y + alpha * U[4].y + beta * U[5].y;
  df[2].z = U[3].z + alpha * U[4].z + beta * U[5].z;

  // Calcul de la matrice G et de son inverse
  for (i = 0; i < 3; ++i)
    for (j = i; j < 3; ++j)
      g[i][j] = g[j][i] = df[i].x * df[j].x + df[i].y * df[j].y
	  + df[i].z * df[i].z;

  a[0][0] = (g[1][1] * g[2][2] - g[2][1] * g[2][1]);
  a[1][1] = (g[0][0] * g[2][2] - g[2][0] * g[2][0]);
  a[2][2] = (g[0][0] * g[1][1] - g[1][0] * g[1][0]);
  a[0][1] = a[1][0] = (g[2][1] * g[2][0] - g[0][1] * g[2][2]);
  a[1][2] = a[2][1] = -(g[0][0] * g[2][1] - g[2][0] * g[0][1]);
  a[0][2] = a[2][0] = (g[0][1] * g[1][2] - g[1][1] * g[0][2]);

}

/* Pour un noeud (alpha)  donne on collectionne les M(alpha,beta) */

/* Constante 1 / sqrt(8.0) */

#define INVSQRT8  ( M_SQRT1_2 * 0.5)

int S3_Tipton(int_type nb_cells, int_type nb_nodes,
    const int_type * nb_nodes_per_cell, const int_type * nodes_number,
    double *x, double *y, double *z, const double *weights, const double *relax,
    int_type n_iter) {

#ifdef DEBUG_NP
  double max_depl, depl;
#endif

  int err = 0;

  int_type vtx, cell;
  const int_type *f_node_cell;
  int iter, i, j;

  int_type *nb_nodes_node = NULL;
  double *x0, *y0, *z0, *m_sum;

  nb_nodes_node = (int_type *) s3Malloc(nb_nodes * sizeof(nb_nodes_node[0]));
  if (nb_nodes_node == NULL)
    goto epilogue;

  x0 = y0 = z0 = m_sum = NULL;
  x0 = (double *) s3Malloc(nb_nodes * sizeof(x0[0]));
  if (x0 == NULL)
    goto epilogue;
  y0 = (double *) s3Malloc(nb_nodes * sizeof(y0[0]));
  if (y0 == NULL)
    goto epilogue;
  z0 = (double *) s3Malloc(nb_nodes * sizeof(z0[0]));
  if (z0 == NULL)
    goto epilogue;
  m_sum = (double *) s3Malloc(nb_nodes * sizeof(m_sum[0]));
  if (m_sum == NULL)
    goto epilogue;

  for (iter = 0; iter < n_iter; ++iter) {
    for (vtx = 0; vtx < nb_nodes; ++vtx) {
      x0[vtx] = y0[vtx] = z0[vtx] = m_sum[vtx] = 0.0;
    }

    for (f_node_cell = nodes_number, cell = 0; cell < nb_cells; ++cell) {
      s3_cell_3d_t my_cell;
      s3_edge_t my_edges[12];
      int nb_edges, edge;

      my_cell.nb_nodes = nb_nodes_per_cell[cell];
      for (vtx = 0; vtx < nb_nodes_per_cell[cell]; ++vtx)
	my_cell.vtx[vtx] = *(f_node_cell++);

      switch (my_cell.nb_nodes) {
      case 4:
	/* Cas du tetraedre */
      {
	point_t u[4];
	double g[3][3];
	double a11, a12, a23, a22, a13, a33;
	double matrix[4][4];

	u[0].x = x[my_cell.vtx[3]];
	u[0].y = y[my_cell.vtx[3]];
	u[0].z = z[my_cell.vtx[3]];
	for (i = 1; i < 4; ++i) {
	  u[i].x = M_SQRT1_2 * (x[my_cell.vtx[i - 1]] - u[0].x);
	  u[i].y = M_SQRT1_2 * (y[my_cell.vtx[i - 1]] - u[0].y);
	  u[i].z = M_SQRT1_2 * (z[my_cell.vtx[i - 1]] - u[0].z);
	}

	for (i = 0; i < 3; ++i)
	  for (j = i; j < 3; ++j)
	    g[i][j] = g[j][i] = 2.0
		* (u[i + 1].x * u[j + 1].x + u[i + 1].y * u[j + 1].y
		    + u[i + 1].z * u[j + 1].z);

	a11 = (g[1][1] * g[2][2] - g[2][1] * g[2][1]);
	a22 = (g[0][0] * g[2][2] - g[2][0] * g[2][0]);
	a33 = (g[0][0] * g[1][1] - g[1][0] * g[1][0]);
	a12 = (g[2][1] * g[2][0] - g[0][1] * g[2][2]);
	a23 = -(g[0][0] * g[2][1] - g[2][0] * g[0][1]);
	a13 = (g[0][1] * g[1][2] - g[1][1] * g[0][2]);

	matrix[0][0] = a11;
	matrix[0][1] = a12;
	matrix[0][2] = a13;
	matrix[0][3] = -a11 - a12 - a13;

	matrix[1][0] = a12;
	matrix[1][1] = a22;
	matrix[1][2] = a23;
	matrix[1][3] = -a12 - a22 - a23;

	matrix[2][0] = a13;
	matrix[2][1] = a23;
	matrix[2][2] = a33;
	matrix[2][3] = -a13 - a23 - a33;

	matrix[3][0] = -a11 - a12 - a13;
	matrix[3][1] = -a12 - a22 - a23;
	matrix[3][2] = -a13 - a23 - a33;
	matrix[3][3] = a11 + a22 + a33 + 2.0 * (a12 + a13 + a23);

	for (i = 0; i < 4; ++i) {

	  int id1 = my_cell.vtx[i];
	  for (j = i; j < 4; ++j) {
	    double w = weights[cell] * (1.0 / 6.0);
	    int id2 = my_cell.vtx[j];

	    if (id1 == id2)
	      m_sum[id1] += w * matrix[i][i];
	    else {
	      x0[id1] += w * x[id2] * matrix[i][j];
	      x0[id2] += w * x[id1] * matrix[j][i];

	      y0[id1] += w * y[id2] * matrix[i][j];
	      y0[id2] += w * y[id1] * matrix[j][i];

	      z0[id1] += w * z[id2] * matrix[i][j];
	      z0[id2] += w * z[id1] * matrix[j][i];
	    }
	  }
	}
      }
	break;

      case 5:
	break;
	/* Cas de la pyramide */
	{
	  point_t u[5];
	  double matrix[5][5], a[3][3];
	  double a11, a12, a23, a22, a13, a33;

	  u[0].x = x[my_cell.vtx[0]];
	  u[0].y = y[my_cell.vtx[0]];
	  u[0].z = z[my_cell.vtx[0]];

	  u[1].x = M_SQRT1_2 * (x[my_cell.vtx[1]] - u[0].x);
	  u[1].y = M_SQRT1_2 * (y[my_cell.vtx[1]] - u[0].y);
	  u[1].z = M_SQRT1_2 * (z[my_cell.vtx[1]] - u[0].z);

	  u[2].x = M_SQRT1_2 * (x[my_cell.vtx[3]] - u[0].x);
	  u[2].y = M_SQRT1_2 * (y[my_cell.vtx[3]] - u[0].y);
	  u[2].z = M_SQRT1_2 * (z[my_cell.vtx[3]] - u[0].z);

	  u[3].x = M_SQRT1_2 * (x[my_cell.vtx[4]] - u[0].x);
	  u[3].y = M_SQRT1_2 * (y[my_cell.vtx[4]] - u[0].y);
	  u[3].z = M_SQRT1_2 * (z[my_cell.vtx[4]] - u[0].z);

	  u[4].x = 0.5
	      * (u[0].x - x[my_cell.vtx[1]] + x[my_cell.vtx[2]]
		  - x[my_cell.vtx[3]]);
	  u[4].y = 0.5
	      * (u[0].y - y[my_cell.vtx[1]] + y[my_cell.vtx[2]]
		  - y[my_cell.vtx[3]]);
	  u[4].z = 0.5
	      * (u[0].z - z[my_cell.vtx[1]] + z[my_cell.vtx[2]]
		  - z[my_cell.vtx[3]]);

	  // Calcul de la matrice de rigidite
	  for (i = 0; i < 5; ++i)
	    for (j = 0; j < 5; ++j)
	      matrix[i][j] = 0.0;

	  // Calcul des coefficient AIJ pour le point 1/4,1/2,1/4
	  computeAijPyramide(0.25, 0.5, 0.25, u, a);

	  a11 = a[0][0];
	  a12 = a[0][1];
	  a13 = a[0][2];
	  a22 = a[1][1];
	  a23 = a[1][2];
	  a33 = a[2][2];

	  matrix[0][0] = a11 / 16. + (3 * a12) / 8. + (3 * a13) / 8.
	      + (9 * a22) / 16. + (9 * a23) / 8. + (9 * a33) / 16.;
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
	  matrix[1][2] = (3 * a11) / 16. - a12 / 8. - a13 / 8. - a22 / 16.
	      - a23 / 8. - a33 / 16.;
	  matrix[1][3] = (-3 * a11) / 16. + (3 * a12) / 8. + a13 / 8.
	      - (3 * a22) / 16. - a23 / 8. + a33 / 16.;
	  matrix[1][4] = a13 / 4. - a23 / 4. - a33 / 4.;

	  matrix[2][0] = matrix[0][2];
	  matrix[2][1] = matrix[1][2];
	  matrix[2][2] = (9 * a11) / 16. + (3 * a12) / 8. + (3 * a13) / 8.
	      + a22 / 16. + a23 / 8. + a33 / 16.;
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

	  matrix[0][0] += (9 * a11) / 16. + (3 * a12) / 8. + (9 * a13) / 8.
	      + a22 / 16. + (3 * a23) / 8. + (9 * a33) / 16.;
	  matrix[0][1] += (-9 * a11) / 16. + (3 * a12) / 8. - (3 * a13) / 8.
	      + (3 * a22) / 16. + (5 * a23) / 8. + (3 * a33) / 16.;
	  matrix[0][2] += (-3 * a11) / 16. - (5 * a12) / 8. - (3 * a13) / 8.
	      - (3 * a22) / 16. - (5 * a23) / 8. - (3 * a33) / 16.;
	  matrix[0][3] += (3 * a11) / 16. - a12 / 8. + (3 * a13) / 8.
	      - a22 / 16. - a23 / 8. + (3 * a33) / 16.;
	  matrix[0][4] += (-3 * a13) / 4. - a23 / 4. - (3 * a33) / 4.;

	  matrix[1][0] += (-9 * a11) / 16. + (3 * a12) / 8. - (3 * a13) / 8.
	      + (3 * a22) / 16. + (5 * a23) / 8. + (3 * a33) / 16.;
	  matrix[1][1] += (9 * a11) / 16. - (9 * a12) / 8. - (3 * a13) / 8.
	      + (9 * a22) / 16. + (3 * a23) / 8. + a33 / 16.;
	  matrix[1][2] += (3 * a11) / 16. + (3 * a12) / 8. + a13 / 8.
	      - (9 * a22) / 16. - (3 * a23) / 8. - a33 / 16.;
	  matrix[1][3] += (-3 * a11) / 16. + (3 * a12) / 8. - a13 / 8.
	      - (3 * a22) / 16. + a23 / 8. + a33 / 16.;
	  matrix[1][4] += (3 * a13) / 4. - (3 * a23) / 4. - a33 / 4.;

	  matrix[2][0] += (-3 * a11) / 16. - (5 * a12) / 8. - (3 * a13) / 8.
	      - (3 * a22) / 16. - (5 * a23) / 8. - (3 * a33) / 16.;
	  matrix[2][1] += (3 * a11) / 16. + (3 * a12) / 8. + a13 / 8.
	      - (9 * a22) / 16. - (3 * a23) / 8. - a33 / 16.;
	  matrix[2][2] += a11 / 16. + (3 * a12) / 8. + a13 / 8.
	      + (9 * a22) / 16. + (3 * a23) / 8. + a33 / 16.;
	  matrix[2][3] += -a11 / 16. - a12 / 8. - a13 / 8. + (3 * a22) / 16.
	      - a23 / 8. - a33 / 16.;
	  matrix[2][4] += a13 / 4. + (3 * a23) / 4. + a33 / 4.;

	  matrix[3][0] += (3 * a11) / 16. - a12 / 8. + (3 * a13) / 8.
	      - a22 / 16. - a23 / 8. + (3 * a33) / 16.;
	  matrix[3][1] += (-3 * a11) / 16. + (3 * a12) / 8. - a13 / 8.
	      - (3 * a22) / 16. + a23 / 8. + a33 / 16.;
	  matrix[3][2] += -a11 / 16. - a12 / 8. - a13 / 8. + (3 * a22) / 16.
	      - a23 / 8. - a33 / 16.;
	  matrix[3][3] += a11 / 16. - a12 / 8. + a13 / 8. + a22 / 16. - a23 / 8.
	      + a33 / 16.;
	  matrix[3][4] += -a13 / 4. + a23 / 4. - a33 / 4;

	  matrix[4][0] += (-3 * a13) / 4. - a23 / 4. - (3 * a33) / 4.;
	  matrix[4][1] += (3 * a13) / 4. - (3 * a23) / 4. - a33 / 4.0;
	  matrix[4][2] += a13 / 4. + (3 * a23) / 4. + a33 / 4.;
	  matrix[4][3] += -a13 / 4. + a23 / 4. - a33 / 4;
	  matrix[4][4] += a33;
	  for (i = 0; i < 5; ++i) {

	    int id1 = my_cell.vtx[i];
	    for (j = i; j < 5; ++j) {
	      double w = weights[cell] * (1.0 / 12.0);
	      int id2 = my_cell.vtx[j];

	      if (id1 == id2)
		m_sum[id1] += w * matrix[i][i];
	      else {
		x0[id1] += w * x[id2] * matrix[i][j];
		x0[id2] += w * x[id1] * matrix[j][i];

		y0[id1] += w * y[id2] * matrix[i][j];
		y0[id2] += w * y[id1] * matrix[j][i];

		z0[id1] += w * z[id2] * matrix[i][j];
		z0[id2] += w * z[id1] * matrix[j][i];
	      }
	    }
	  }
	}
	break;

      case 6:
	/* Allons y pour un prisme */
      {
	point_t u[6];
	double matrix[6][6];
	double a[3][3];
	double a11, a12, a23, a22, a13, a33;

	/* U = T Nodes */

	/* Changement de variables */
	for (i = 0; i < 6; ++i) {
	  double tx, ty, tz;
	  tx = ty = tz = 0.0;

	  for (j = 0; j < 6; ++j) {
	    double tij = kTPenta[i][j];
	    tx += x[my_cell.vtx[j]] * tij;
	    ty += y[my_cell.vtx[j]] * tij;
	    tz += z[my_cell.vtx[j]] * tij;
	  }
	  u[i].x = tx;
	  u[i].y = ty;
	  u[i].z = tz;
	}

	for (i = 0; i < 6; ++i)
	  for (j = 0; j < 6; ++j)
	    matrix[i][j] = 0.0;

	// Calcul des coefficient AIJ pour le point 1/3,1/3,1/2
	computeAijPenta(1.0 / 3.0, 1.0 / 3.0, 0.5, u, a);

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
	matrix[0][4] = -a11 * 0.25 - a12 * 0.25 + a13 / 3. + a23 / 6.
	    - a33 / 9.;
	matrix[0][5] = -a12 * 0.25 + a13 / 6. - a22 * 0.25 + a23 / 3.
	    - a33 / 9.;

	matrix[1][0] = matrix[0][1];
	matrix[1][1] = a11 * 0.25 + a13 / 3. + a33 / 9.;
	matrix[1][2] = a12 * 0.25 + a13 / 6. + a23 / 6. + a33 / 9.;
	matrix[1][3] = -a11 * 0.25 - a12 * 0.25 - a13 / 3. - a23 / 6.
	    - a33 / 9.;
	matrix[1][4] = a11 * 0.25 - a33 / 9.;
	matrix[1][5] = a12 * 0.25 - a13 / 6. + a23 / 6. - a33 / 9.;

	matrix[2][0] = matrix[0][2];
	matrix[2][1] = matrix[1][2];
	matrix[2][2] = a22 * 0.25 + a23 / 3. + a33 / 9.;
	matrix[2][3] = -a12 * 0.25 - a13 / 6. - a22 * 0.25 - a23 / 3.
	    - a33 / 9.;
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

	for (i = 0; i < 6; ++i) {

	  int id1 = my_cell.vtx[i];
	  for (j = i; j < 6; ++j) {
	    double w = weights[cell] * (1.0 / 2.0) * matrix[i][j];
	    int id2 = my_cell.vtx[j];

	    if (id1 == id2)
	      m_sum[id1] += w;
	    else {
	      x0[id1] += w * x[id2];
	      x0[id2] += w * x[id1];

	      y0[id1] += w * y[id2];
	      y0[id2] += w * y[id1];

	      z0[id1] += w * z[id2];
	      z0[id2] += w * z[id1];
	    }
	  }
	}
      }
	break;

      case 8:
	/* Cas de l'hexaedre */
      {

	point_t u[8];
	double g[3][3], a[3][3], matrix[8][8];

	for (i = 0; i < 8; ++i) {
	  double tx, ty, tz;
	  tx = ty = tz = 0.0;

	  for (j = 0; j < 8; ++j) {
	    double tij = kTHexa[i][j];
	    tx += x[my_cell.vtx[j]] * tij;
	    ty += y[my_cell.vtx[j]] * tij;
	    tz += z[my_cell.vtx[j]] * tij;
	  }
	  u[i].x = tx * INVSQRT8;
	  u[i].y = ty * INVSQRT8;
	  u[i].z = tz * INVSQRT8;
	}
	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++)
	    g[i][j] = 0.5
		* (u[i + 1].x * u[j + 1].x + u[i + 1].y * u[j + 1].y
		    + u[i + 1].z * u[j + 1].z);

	a[0][0] = (g[1][1] * g[2][2] - g[2][1] * g[2][1]);
	a[1][1] = (g[0][0] * g[2][2] - g[2][0] * g[2][0]);
	a[2][2] = (g[0][0] * g[1][1] - g[1][0] * g[1][0]);
	a[0][1] = a[1][0] = (g[2][1] * g[2][0] - g[0][1] * g[2][2]);
	a[1][2] = a[2][1] = -(g[0][0] * g[2][1] - g[2][0] * g[0][1]);
	a[0][2] = a[2][0] = (g[0][1] * g[1][2] - g[1][1] * g[0][2]);

	/* Ici on utilise le fait que l'on a une quadrature a un point */

	/* matrix = Transpose[T]PT */

	matrix[0][0] = 0.25
	    * (a[0][0] - a[0][1] - a[0][2] + a[1][1] + a[1][2] + a[2][2]);
	matrix[1][0] = matrix[0][1] = -a[1][1] * 0.25;
	matrix[2][0] = matrix[0][2] = a[0][1] * 0.25;
	matrix[3][0] = matrix[0][3] = -a[0][0] * 0.25;
	matrix[4][0] = matrix[0][4] = -a[2][2] * 0.25;
	matrix[5][0] = matrix[0][5] = -a[1][2] * 0.25;
	matrix[6][0] = matrix[0][6] = 0.0;
	matrix[7][0] = matrix[0][7] = a[0][2] * 0.25;

	matrix[1][1] = 0.25
	    * (a[0][0] + a[0][1] - a[0][2] + a[1][1] - a[1][2] + a[2][2]);
	matrix[2][1] = matrix[1][2] = -a[0][0] * 0.25;
	matrix[3][1] = matrix[1][3] = -a[0][1] * 0.25;
	matrix[4][1] = matrix[1][4] = a[1][2] * 0.25;
	matrix[5][1] = matrix[1][5] = -a[2][2] * 0.25;
	matrix[6][1] = matrix[1][6] = a[0][2] * 0.25;
	matrix[7][1] = matrix[1][7] = 0.0;

	matrix[2][2] = 0.25
	    * (a[0][0] - a[0][1] + a[0][2] + a[1][1] - a[1][2] + a[2][2]);
	matrix[2][3] = matrix[3][2] = -a[1][1] * 0.25;
	matrix[2][4] = matrix[4][2] = 0.0;
	matrix[2][5] = matrix[5][2] = -a[0][2] * 0.25;
	matrix[2][6] = matrix[6][2] = -a[2][2] * 0.25;
	matrix[2][7] = matrix[7][2] = a[1][2] * 0.25;

	matrix[3][3] = 0.25
	    * (a[0][0] + a[0][1] + a[0][2] + a[1][1] + a[1][2] + a[2][2]);
	matrix[3][4] = matrix[4][3] = -a[0][2] * 0.25;
	matrix[3][5] = matrix[5][3] = 0.0;
	matrix[3][6] = matrix[6][3] = -a[1][2] * 0.25;
	matrix[3][7] = matrix[7][3] = -a[2][2] * 0.25;

	matrix[4][4] = 0.25
	    * (a[0][0] - a[0][1] + a[0][2] + a[1][1] - a[1][2] + a[2][2]);
	matrix[4][5] = matrix[5][4] = -a[1][1] * 0.25;
	matrix[4][6] = matrix[6][4] = a[0][1] * 0.25;
	matrix[4][7] = matrix[7][4] = -a[0][0] * 0.25;

	matrix[5][5] = 0.25
	    * (a[0][0] + a[0][1] + a[0][2] + a[1][1] + a[1][2] + a[2][2]);
	matrix[5][6] = matrix[6][5] = -a[0][0] * 0.25;
	matrix[5][7] = matrix[7][5] = -a[0][1] * 0.25;

	matrix[6][6] = 0.25
	    * (a[0][0] - a[0][1] - a[0][2] + a[1][1] + a[1][2] + a[2][2]);
	matrix[6][7] = matrix[7][6] = -a[1][1] * 0.25;

	matrix[7][7] = 0.25
	    * (a[0][0] + a[0][1] - a[0][2] + a[1][1] - a[1][2] + a[2][2]);

	for (i = 0; i < 8; ++i) {
	  int id1 = my_cell.vtx[i];
	  for (j = i; j < 8; ++j) {
	    int id2 = my_cell.vtx[j];
	    double w = weights[cell] * matrix[i][j];

	    if (id1 == id2)
	      m_sum[id1] += w;
	    else {
	      x0[id1] += w * x[id2];
	      x0[id2] += w * x[id1];

	      y0[id1] += w * y[id2];
	      y0[id2] += w * y[id1];

	      z0[id1] += w * z[id2];
	      z0[id2] += w * z[id1];
	    }
	  }
	}
      }
	break;
      default:
	break;
      }
    }

#ifdef DEBUG_NP
    max_depl = 0.0;
#endif

    for (vtx = 0; vtx < nb_nodes; ++vtx) {
      double rel = relax[vtx];
      int seen = 0;
      if (rel != 0.0) {
	double m_x0, m_y0, m_z0, my_self;
	my_self = -m_sum[vtx];

	if (my_self != 0.0) {
	  m_x0 = x0[vtx] / my_self;
	  m_y0 = y0[vtx] / my_self;
	  m_z0 = z0[vtx] / my_self;
#ifdef DEBUG_NP
	  depl = (x[vtx] - m_x0) * (x[vtx] - m_x0)
	  + (y[vtx] - m_y0) * (y[vtx] - m_y0)
	  + (z[vtx] - m_z0) * (z[vtx] - m_z0);
	  if (depl > max_depl)
	  max_depl = depl;
#endif
	  x[vtx] = (1 - rel) * x[vtx] + rel * (m_x0);
	  y[vtx] = (1 - rel) * y[vtx] + rel * (m_y0);
	  z[vtx] = (1 - rel) * z[vtx] + rel * (m_z0);
	}
      }
    }

#ifdef DEBUG_NP
    fprintf (stderr, "Iteration %d max depl %e\n", iter, sqrt (max_depl));
#endif
  }

  err = 1; /* All went well */
  epilogue: if (m_sum != NULL)
    s3Free(m_sum);
  if (z0 != NULL)
    s3Free(z0);
  if (y0 != NULL)
    s3Free(y0);
  if (x0 != NULL)
    s3Free(x0);
  if (nb_nodes_node != NULL)
    s3Free(nb_nodes_node);

  return err;
}
