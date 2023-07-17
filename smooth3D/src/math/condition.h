// condition.h
//
// Calcul du conditionnement des matrices jacobiennes des sommets
// Exprimes pour des maillages surfaciques et volumiques
// Les Gradients sont aussi calcules
// $Id: condition.h,v 1.4 2003/06/25 14:19:37 weilljc Exp $

#ifndef CONDITION_H
#define CONDITION_H

#include "Real3.h"
namespace Smooth3D {

// Conditionnement de la matrice jacobienne dans un maillage surfacique
// le sommet est v1 et deux bras partent vers les points v2 et v3

inline Real condition2D(Real3 v1, Real3 v2, Real3 v3) {
  Real3 u = v2 - v1, v = v3 - v1;
  return (u.abs2() + v.abs2()) / math::vecMul(u, v).abs();
}

// Conditionnement de la matrice jacobienne dans un maillage volumique
// le sommet est v1 et deux bras partent vers les points v2, v3 et v4.

inline Real condition3D(Real3 v1, Real3 v2, Real3 v3, Real3 v4) {
  Real3 u = v2 - v1, v = v3 - v1, w = v4 - v1;

  Real a = (u.abs2() + v.abs2() + w.abs2());
  a *= (math::vecMul(u, v).abs2() + math::vecMul(u, w).abs2()
      + math::vecMul(v, w).abs2());

  return std::sqrt(a) / math::abs(math::mixteMul(u, v, w));

}


// gradient du conditionnement de la matrice jacobienne par rapport a v1

inline Real3 gradCond2DV1(Real3 v1, Real3 v2, Real3 v3) {
  Real3 u = v2 - v1, v = v3 - v1;
  Real3 prod_vect = math::vecMul(u, v);
  Real den = prod_vect.abs();
  Real num = u.abs2() + v.abs2();
  Real3 grad = (-2.0 * (u + v)
      + num * math::vecMul(v - u, prod_vect) / den / den) / den;

  return grad;
}

// gradient du conditionnement de la matrice jacobienne par rapport a v2

inline Real3 gradCond2DV2(Real3 v1, Real3 v2, Real3 v3) {
  Real3 u = v2 - v1, v = v3 - v1;
  Real3 prod_vect = math::vecMul(u, v);
  Real den = prod_vect.abs();
  Real num = u.abs2() + v.abs2();

  Real3 grad = (2.0 * u - num * math::vecMul(v, prod_vect) / den / den) / den;
  return grad;
}

inline Real3 gradCond3DV1(Real3 v1, Real3 v2, Real3 v3, Real3 v4) {
  Real3 u = v2 - v1, v = v3 - v1, w = v4 - v1;
  Real3 b1 = math::vecMul(u, v);
  Real3 b2 = math::vecMul(u, w);
  Real3 b3 = math::vecMul(v, w);

  // Le numerateur est la racine carre d'un produit
  // On calcule chaque terme du produit et son gradient

  // premier terme
  Real a = (u.abs2() + v.abs2() + w.abs2());
  Real3 da = -2.0 * (u + v + w);

  // second terme
  Real b = b1.abs2() + b2.abs2() + b3.abs2();
  Real3 db = 2.0
      * (math::vecMul(u - v, b1) + math::vecMul(u - w, b2)
	  + math::vecMul(v - w, b3));

  // Calcul du numerateur et de son gradient
  Real num = std::sqrt(a * b);
  Real3 dnum = 0.5 * (da * b + a * db) / num;

  // Calcul du denominateur et de son gradient
  Real den = math::mixteMul(u, v, w);
  Real3 dden = math::vecMul(v3, v2) + math::vecMul(v2, v4)
      + math::vecMul(v4, v3);

  if (den < 0) {
    den = -den;
    dden = -dden;
  }

  Real3 grad = (dnum * den - num * dden) / den / den;

  return grad;

}

inline Real3 gradCond3DV2(Real3 v1, Real3 v2, Real3 v3, Real3 v4) {
  Real3 u = v2 - v1, v = v3 - v1, w = v4 - v1;
  Real3 b1 = math::vecMul(u, v);
  Real3 b2 = math::vecMul(u, w);
  Real3 b3 = math::vecMul(v, w);

  // Le numerateur est la racine carre d'un produit
  // On calcule chaque terme du produit et son gradient

  // premier terme
  Real a = (u.abs2() + v.abs2() + w.abs2());
  Real3 da = 2.0 * (u);

  // second terme
  Real b = b1.abs2() + b2.abs2() + b3.abs2();
  Real3 db = 2.0 * (math::vecMul(v, b1) + math::vecMul(w, b2));

  // Calcul du numerateur et de son gradient
  Real num = std::sqrt(a * b);
  Real3 dnum = 0.5 * (da * b + a * db) / num;

  // Calcul du denominateur et de son gradient
  Real den = math::mixteMul(u, v, w);
  Real3 dden = math::vecMul(v1, v3) + math::vecMul(v4, v1)
      + math::vecMul(v3, v4);
  if (den < 0) {
    den = -den;
    dden = -dden;
  }

  Real3 grad = (dnum * den - num * dden) / den / den;

  return grad;
}

}

#endif

