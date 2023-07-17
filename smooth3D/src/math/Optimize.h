// Gestion des algorithmes d'optimisations
//

#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include "Real3.h"


namespace Optimize
{


using Smooth3D::Real;
  // Interface des fonctions a minimiser

  template < class T >class Fonction
    {
      public:

        Fonction ()
      {
      }

      virtual ~ Fonction ()
      {
      }

      virtual Real evalF (const T & x) = 0;
      virtual void evalDF (const T & x, T & gradient) = 0;
      virtual Real evalFDF (const T & x, T & gradient) = 0;
    };

  enum
    {
      OPT_SUCCESS, OPT_ENOMEM, OPT_ENOPROG
    };

  template < class T >class Minimizer
    {

      protected:
        virtual void _restart () = 0;

      virtual void _set (const T & x,Real & val, T & gradient,
			Real step_size, Real tol) = 0;

      virtual int _iterate (T & x, Real & val, T & gradient, T & dx) = 0;

      T m_x;
      T m_gradient;
      Real m_val;
      T m_dx;

      public:


        Minimizer ()
      {
      }
      virtual ~ Minimizer ()
      {
      }

      void globalSet (const T & x, Real step_size, Real tol)
      {
	m_x = x;
	m_dx = T::null();
	_set (m_x, m_val, m_gradient, step_size, tol);
      }

      int globalIterate ()
      {
	return _iterate (m_x, m_val, m_gradient, m_dx);
      }

      T solution () const
      {
	return m_x;
      }

      T gradient () const 
      {
	return m_gradient;
      }

      T  dx () const
      {
	return m_dx;
      }
      
      Real value () const
      {
        return m_val;
      }

      bool convergenceP (Real epsilon)
      {
	return m_gradient.abs () < epsilon;
      }
    };
};
#endif
