
//
// Minimization directionel
// LineSearch Algorithms

#include "Optimize.h"

#ifndef LINESEARCH_H
#define LINESEARCH_H

namespace Optimize
{


  template <class T> class LineSearch {
        
  // fa < fc
  // Cherche un point entre stepa et stepc plus petit que fa
  
  inline void
    _intermediatePoint (Fonction < T > *fdf,
			const T & x, const T & p,
			Real lambda, Real pg, Real stepa, Real stepc,
			Real fa, Real fc,
			T & x1, T & dx, T & gradient,
			Real & step, Real & f)
  {
    Real stepb = 0, fb;
    do
      {
	Real u = math::abs (pg * lambda * stepc);
	  stepb = 0.5 * stepc * u / ((fc - fa) + u);

	  x1 = x - stepb * lambda * p;

	  fb = fdf->evalF (x1);

	if (fb < fa || stepb <= 1.0e-13)
	    break;
	  fc = fb;
	  stepc = stepb;
      }
    while (true);

      step = stepb;
      f = fb;
      fdf->evalDF (x1, gradient);
  }

   public:
// Starting at (x0, f0) move along the direction p to find
// a minimum f (x0 - lambda * p), returning the new point x1 = x0 - lambda * p
// f1 = f(x1) et g1 = grad(f)(x1)

  inline void operator () (Fonction < T > *fdf, const T & x, const T & p,
		Real lambda, Real pg, 
                Real stepa,  Real stepc,
		Real fa,  Real fc, Real tol,
		T & x1, T & dx1,
		T & x2, T & dx2,
		T & gradient,
		Real & step, Real & f, Real & gnorm)
  {
    
   Real fb, stepb;
     _intermediatePoint  (fdf, x, p, lambda, pg,
                         stepa, stepc, 
                         fa, fc,
                         x1, dx1, gradient,
                         stepb, fb);
    x2 = x1;
    dx2 = dx1;
   
    gnorm = gradient.abs ();
    f = fb;
    step = stepb;
    
    if (step == 0)
        return;

    Real u = stepb;
    Real v = stepa;
    Real w = stepc;

    Real fu = fb;
    Real fv = fa;
    Real fw = fc;

    Real old1 = math::abs (v - u);
    Real old2 = math::abs (w - v);

    for (int iter = 0; iter < 10; iter++)
      {
	Real dw = w - u;
	Real dv = v - u;
	Real du = 0.0;

	Real e1 = ((fv - fu) * dw * dw + (fu - fw) * dv * dv);
	Real e2 = 2.0 * ((fv - fu) * dw + (fu - fw) * dv);

	if (e2 != 0.0)
	    du = e1 / e2;

	Real stepm;

	if (du > 0 && du < (stepc - stepb)
	    && math::abs (du) < 0.5 * old2)
	    stepm = u + du;
	else if (du < 0 && du > (stepa - stepb)
		 && math::abs (du) < 0.5 * old2)
	    stepm = u + du;
	else if ((stepc - stepb) > (stepb - stepa))
	    stepm = 0.38 * (stepc - stepb) + stepb;
       
	else
	    stepm = stepb - 0.38 * (stepb - stepa);

	  dx1 = (-stepm * lambda) * p;
	  x1 = x + dx1;

	Real fm = fdf->evalF (x1);

	if (fm > fb)
	  {
	    if (fm < fv)
	      {
		w = v;
		v = stepm;
		fw = fv;
		fv = fm;
	      }
	    else if (fm < fw)
	      {
		w = stepm;
		fw = fm;
	      }
	    if (stepm < stepb)
	      {
		stepa = stepm;
		fa = fm;
	      }
	    else
	      {
		stepc = stepm;
		fc = fm;
	      }
	    continue;
	  }
	else			// fm <= fb
	  {
	    old2 = old1;
	    old1 = math::abs (u - stepm);
	    w = v;
	    v = u;
	    u = stepm;

	    fw = fv;
	    fv = fu;
	    fu = fm;

	    x2 = x1;
	    dx2 = dx1;
	    fdf->evalDF (x1, gradient);
	    Real pg = math::scaMul (p, gradient);
	    Real gnorm1 = gradient.abs ();

            f = fm;
	    step = stepm;
	    gnorm = gnorm1 ;

	    if (gnorm1 == 0.0)
	      return;
	    if (math::abs (pg * lambda / gnorm1) < tol)
	      return;
	    if (stepm < stepb)
	      {
		stepc = stepb;
		fc = fb;
		stepb = stepm;
		fb = fm;
	      }
	    else
	      {
		stepa = stepb;
		fa = fb;
		stepb = stepm;
		fb = fm;
	      }
	  }
      }				// iteration
    // Arrive au maximum des iterations
  }				// fin de la fonction
  };

};
#endif
