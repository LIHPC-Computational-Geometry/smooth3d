/*
 * MeanRatio.h
 *
 *  Created on: 20 nov. 2013
 *      Author: weilljc
 */
#ifndef MEANRATIO_H
#define MEANRATIO_H

#include <cfloat> // for DBL_MIN

#include "Real3.h"
namespace Smooth3D {

#define M_SQRT1_3 5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0)*/

/*****************************************************************************/
/* The following set of functions reference triangular elements to an        */
/* equilateral triangle in the plane defined by the normal.  They are        */
/* used when assessing the quality of a triangular element.  A zero          */
/* return value indicates success, while a nonzero value indicates failure.  */
/*****************************************************************************/

/*****************************************************************************/
/* Function evaluation requires 44 flops.                                    */
/*   Reductions possible when b == 1 or c == 1                               */
/*****************************************************************************/
inline bool mFvn2e(Real &obj, const Real3 x[3], const Real3 &n) {
  Real matr[9], f;
  Real g;

  /* Calculate M = [A*inv(W) n] */
  matr[0] = x[1].x() - x[0].x();
  matr[1] = (2.0 * x[2].x() - x[1].x() - x[0].x()) * M_SQRT1_3;
  matr[2] = n.x();

  matr[3] = x[1].y() - x[0].y();
  matr[4] = (2.0 * x[2].y() - x[1].y() - x[0].y()) * M_SQRT1_3;
  matr[5] = n.y();

  matr[6] = x[1].z() - x[0].z();
  matr[7] = (2.0 * x[2].z() - x[1].z() - x[0].z()) * M_SQRT1_3;
  matr[8] = n.z();

  /* Calculate det(M). */
  g = matr[0] * (matr[4] * matr[8] - matr[5] * matr[7])
      + matr[3] * (matr[2] * matr[7] - matr[1] * matr[8])
      + matr[6] * (matr[1] * matr[5] - matr[2] * matr[4]);

  if (g < 0.0)
    return 0.0; // Element is inverted

  /* Calculate norm(M). */
  f = matr[0] * matr[0] + matr[1] * matr[1] + matr[3] * matr[3]
      + matr[4] * matr[4] + matr[6] * matr[6] + matr[7] * matr[7];

  /* Calculate objective function. */
  obj = 2.0 * g / f;
  return (true);
}

/*****************************************************************************/
/* The following set of functions reference triangular elements to an        */
/* right triangle in the plane defined by the normal.  They are used when    */
/* assessing the quality of a quadrilateral elements.  A zero return value   */
/* indicates success, while a nonzero value indicates failure.               */
/*****************************************************************************/

/*****************************************************************************/
/* Function evaluation -- requires 41 flops.                                 */
/*   Reductions possible when b == 1, c == 1, or d == 1                      */
/*****************************************************************************/
inline bool mFcn2i(Real &obj, const Real3 x[3], const Real3 &n,
    const Real3 &d) {
  Real matr[9];
  Real f;
  Real g;

  /* Calculate M = A*inv(W). */
  matr[0] = d.x() * (x[1].x() - x[0].x());
  matr[1] = d.y() * (x[2].x() - x[0].x());
  matr[2] = n.x();

  matr[3] = d.x() * (x[1].y() - x[0].y());
  matr[4] = d.y() * (x[2].y() - x[0].y());
  matr[5] = n.z();

  matr[6] = d.x() * (x[1].z() - x[0].z());
  matr[7] = d.y() * (x[2].z() - x[0].z());
  matr[8] = n.z();

  /* Calculate det(M). */
  g = matr[0] * (matr[4] * matr[8] - matr[5] * matr[7])
      + matr[3] * (matr[2] * matr[7] - matr[1] * matr[8])
      + matr[6] * (matr[1] * matr[5] - matr[2] * matr[4]);


  /* Calculate norm(M). */
  f = matr[0] * matr[0] + matr[1] * matr[1] + matr[3] * matr[3]
      + matr[4] * matr[4] + matr[6] * matr[6] + matr[7] * matr[7];

  /* Calculate objective function. */
  obj = 2.0 * math::abs(g) / f;
  return (true);
}

}

#endif
