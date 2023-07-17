//
// Un exemple d'algorithme de minimisation
// le gradient conjugue version Polak Ribiere
//

#include "Optimize.h"
#include "LineSearch.h"

namespace Optimize {

template<class T, class MyLineSearch = LineSearch<T> > class ConjugateGradientPr: public Minimizer<
    T> {
private:
  int m_iter;
  Real m_step;
  Real m_max_step;
  Real m_tol;
  T m_x1;
  T m_dx1;
  T m_x2;
  Real m_pnorm;
  T m_p;
  Real m_g0norm;
  T m_g0;

  Fonction<T> *m_fdf;
  int m_size;

public:
  ConjugateGradientPr(Fonction<T> *fdf, int size) :
      m_fdf(fdf), m_size(size) {
  }
  virtual ~ ConjugateGradientPr() {
  }

protected:
  virtual void _restart() {
    m_iter = 0;
  }

  virtual void _set(const T & x, Real & val, T & gradient, Real step_size,
      Real tol) {
    m_iter = 0;
    m_step = step_size;
    m_max_step = step_size;
    m_tol = tol;

    val = m_fdf->evalFDF(x, gradient);
    m_p = m_g0 = gradient;

    m_pnorm = m_g0norm = gradient.abs();
  }

  virtual int _iterate(T & x, Real & val, T & gradient, T & dx) {
    Real pnorm = m_pnorm;
    Real g0norm = m_g0norm;

    if (pnorm == 0.0 || g0norm == 0.0 || m_step == 0.0) {
      dx = T::null();
      return OPT_ENOPROG;
    }

    Real fa = val;
    Real fb, fc;
    Real stepa = 0.0, stepb, stepc = m_step, tol = m_tol;
    Real g1norm;

    Real pg = math::scaMul(m_p, gradient);
    Real dir = (pg >= 0.0) ? 1.0 : -1.0;

    dx = m_p * (-dir * stepc / pnorm);
    m_x1 = x + dx;

    fc = m_fdf->evalF(m_x1);

    if (fc < fa) {
      m_step = stepc * 2.0;
      val = fc;
      x = m_x1;
      m_fdf->evalDF(m_x1, gradient);
      m_g0norm = gradient.abs();
      return OPT_SUCCESS;
    }

    MyLineSearch linesearch;
    linesearch(m_fdf, x, m_p, dir / pnorm, pg, stepa, stepc, fa, fc, tol, m_x1,
	m_dx1, m_x2, dx, gradient, m_step, val, g1norm);

    if (m_step == 0.0)
      return OPT_ENOPROG;

    x = m_x2;
    m_iter = (m_iter + 1) % m_size;

    if (m_iter == 0) {
      m_p = gradient;
      m_pnorm = g1norm;
    } else {
      Real beta = math::scaMul(m_g0 - gradient, gradient) / m_g0.abs2();
      m_p = gradient - beta * m_p;
      m_pnorm = m_p.abs();

    }
    m_g0 = gradient;
    m_g0norm = g1norm;

    return OPT_SUCCESS;
  }
};
}
;
