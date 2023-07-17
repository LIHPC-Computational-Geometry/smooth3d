//
// Limited Memory Broyden-Fletcher-Goldfarb-Shanno
// A Quasi Newton Optimization Algorithm
//


#include "Optimize.h"
#include "LineSearch.h"

namespace Optimize
{
  template < class T ,class MyLineSearch = LineSearch<T> > class BFGS:public Minimizer < T >
  {

    private:
    int m_iter;
    int m_size;
    Real m_step;
    Real m_max_step;
    Real m_tol;
    T m_x1;
    T m_dx1;
    T m_x2;
    T m_g0;
    Real m_g0norm;
    Real m_pnorm;
    T m_p;
    T m_x0;

      Fonction < T > *m_fdf;

    protected:
      virtual int _iterate (T & x, Real & val, T & gradient, T & dx)
    {
      Real pnorm = m_pnorm;
      Real g0norm = m_g0norm;

      if (m_step == 0.0 || pnorm == 0.0 || g0norm == 0.0)
	{
	  dx = T::null ();
	  return OPT_ENOPROG;
	}

      Real fa = val;
      Real fb, fc;
      Real stepa = 0.0, stepb, stepc = m_step, tol = m_tol;
      Real g1norm;

      Real pg = math::scaMul (m_p, gradient);
      Real dir = (pg >= 0.0) ? 1.0 : -1.0;

        dx = m_p * (-dir * stepc / pnorm);
        m_x1 = x + dx;

        fc = m_fdf->evalF (m_x1);

      if (fc < fa)
	{
	  m_step = stepc * 2.0;
	  val = fc;
	  x = m_x1;
	  m_fdf->evalDF (m_x1, gradient);
	  m_g0norm = gradient.abs ();
	  return OPT_SUCCESS;
	}
     
        MyLineSearch linesearch;
        linesearch(m_fdf, x, m_p, dir / pnorm, pg,
                        stepa, stepc, fa, fc, tol,
                        m_x1, m_dx1,
                        m_x2, dx, gradient, m_step, val, g1norm);
      
        if (m_step == 0.0)
        return OPT_ENOPROG;
      

        x = m_x2;
        m_iter = (m_iter + 1) % m_size;

      if (m_iter == 0)
	{
	  m_p = gradient;
	  m_pnorm = g1norm;
	}
      else
	{

	  T dx0 = x - m_x0;
	  T dg0 = gradient - m_g0;

	  Real dxg = math::scaMul (gradient, dx0);
	  Real dgg = math::scaMul (gradient, dg0);
	  Real dxdg = math::scaMul (dx0, dg0);
          
          if (dxdg != 0)
          {

	    Real B = dxg / dxdg;


	    Real A = -(1.0 + dg0.abs2 () / dxdg) * B + dgg / dxdg;
	    m_p = gradient - A * dx0 - B * dg0;
	    m_pnorm = m_p.abs ();
          } else
          {
            m_p = gradient;
            m_pnorm = g1norm;
          }

	}
      m_g0 = gradient;
      m_x0 = x;
      m_g0norm = g1norm;

      return OPT_SUCCESS;
    }

    virtual void _set (const T & x, Real & val, T & gradient,
		      Real step_size, Real tol)
    {
      m_iter = 0;
      m_step = step_size;
      m_tol = tol;
      m_max_step = step_size;

      m_x0 = x;
      val = m_fdf->evalFDF (m_x0, gradient);

      m_p = gradient;
      m_g0 = gradient;
      m_g0norm = m_pnorm = gradient.abs ();

    }

    virtual void _restart ()
    {
      m_iter = 0;
    }

    public:
      BFGS (Fonction < T > *fdf, int size):m_fdf (fdf), m_size (size)
    {
    }
     ~BFGS ()
    {
    }


  };


}
