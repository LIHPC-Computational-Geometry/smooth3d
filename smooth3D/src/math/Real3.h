//
// Classe Real3 Compatible avec celle du code T
//

#ifndef Real3_H
#define Real3_H

#include <cmath>
#include <iostream>

namespace Smooth3D {

typedef double Real;

class Real3 {
private:
  Real m_x, m_y, m_z;

public:

  Real3(Real x0, Real y0, Real z0) :
      m_x(x0), m_y(y0), m_z(z0) {
  }
  Real3() :
      m_x(0.0), m_y(0.0), m_z(0.0) {
  }
  Real3(const Real3 & f) :
      m_x(f.m_x), m_y(f.m_y), m_z(f.m_z) {
  }


  double x() const {
    return m_x;
  }

  double y() const {
    return m_y;
  }

  double z() const {
    return m_z;
  }

  static Real3 null() {
    return Real3(0.0, 0.0, 0.0);
  }

  Real3 operator =(Real3 f) {
    m_x = f.m_x, m_y = f.m_y, m_z = f.m_z;
    return *this;
  }
  Real3 operator =(Real v) {
    m_x = m_y = m_z = v;
    return *this;
  }


  Real abs2() const {
    return m_x * m_x + m_y * m_y + m_z * m_z;
  }

  Real abs() const {
    return std::sqrt(abs2());
  }

  Real3 & operator +=(Real3 b) {
    m_x += b.m_x;
    m_y += b.m_y;
    m_z += b.m_z;
    return *this;
  }
  Real3 & operator -=(Real3 b) {
    m_x -= b.m_x;
    m_y -= b.m_y;
    m_z -= b.m_z;
    return *this;
  }

  Real3 operator +(Real3 b) const {
    return Real3(m_x + b.m_x, m_y + b.m_y, m_z + b.m_z);
  }

  Real3 operator -(Real3 b) const {
    return Real3(m_x - b.m_x, m_y - b.m_y, m_z - b.m_z);
  }

  Real3 operator -() const {
    return Real3(-m_x, -m_y, -m_z);
  }
  bool operator ==(Real3 b) const {
    return m_x == b.m_x && m_y == b.m_y && m_z == b.m_z;
  }
  bool operator !=(Real3 b) const {
    return !operator ==(b);
  }

  Real3 operator *(Real b) const {
    return Real3(m_x * b, m_y * b, m_z * b);
  }

  Real3 operator /(Real b) const {
    return Real3(m_x * (1 / b), m_y * (1 / b), m_z * (1 / b));
  }

  inline bool lexicographicLessThan(Real3 other, Real EPS) const;
};

inline bool Real3::lexicographicLessThan(Real3 other, Real EPS) const {
  if (m_x < other.m_x - EPS)
    return 1;
  if (m_x > other.m_x + EPS)
    return 0;
  if (m_y < other.m_y - EPS)
    return 1;
  if (m_y > other.m_y + EPS)
    return 0;
  if (m_z < other.m_z - EPS)
    return 1;
  if (m_z > other.m_z + EPS)
    return 0;
  return 0;
}

inline std::ostream & operator <<(std::ostream & o, Real3 v) {
  o << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
  return o;

}

inline Real3 operator *(Real c, Real3 b) {
  return Real3(c * b.x(), c * b.y(), c * b.z());
}

class VertexLexicographicLessThan {
  Real m_eps;
public:
  VertexLexicographicLessThan(Real eps = 1.e-6) :
      m_eps(eps) {
  }
  ;
  bool operator ()(Real3 * ent1, Real3 * ent2) const {
    return ent1->lexicographicLessThan(*ent2, m_eps);
  }

};

class Real2 {

private:
  Real m_x, m_y;

public:

  double x() const { return m_x; }
  double y() const { return m_y; }

  Real2(Real x0, Real y0) :
      m_x(x0), m_y(y0) {
  }
  Real2() :
      m_x(0.0), m_y(0.0) {
  }
  Real2(const Real2 & f) :
      m_x(f.m_x), m_y(f.m_y) {
  }



  static Real2 null() {
    return Real2(0.0, 0.0);
  }

  Real2 operator =(Real2 f) {
    m_x = f.m_x, m_y = f.m_y;
    return *this;
  }
  Real2 operator =(Real v) {
    m_x = m_y = v;
    return *this;
  }

  Real abs2() const {
    return m_x * m_x + m_y * m_y;
  }

  Real abs() const {
    return std::sqrt(abs2());
  }

  Real2 & operator +=(Real2 b) {
    m_x += b.m_x;
    m_y += b.m_y;
    return *this;
  }
  Real2 & operator -=(Real2 b) {
    m_x -= b.m_x;
    m_y -= b.m_y;
    return *this;
  }

  Real2 operator +(Real2 b) const {
    return Real2(m_x + b.m_x, m_y + b.m_y);
  }

  Real2 operator -(Real2 b) const {
    return Real2(m_x - b.m_x, m_y - b.m_y);
  }

  Real2 operator -() const {
    return Real2(-m_x, -m_y);
  }
  bool operator ==(Real2 b) const {
    return m_x == b.m_x && m_y == b.m_y;
  }
  bool operator !=(Real2 b) const {
    return !operator ==(b);
  }

  Real2 operator *(Real b) const {
    return Real2(m_x * b, m_y * b);
  }

  Real2 operator /(Real b) const {
    return Real2(m_x * (1 / b), m_y * (1 / b));
  }

};

inline std::ostream & operator <<(std::ostream & o, Real2 v) {
  o << "(" << v.x() << ", " << v.y() << ")";
  return o;

}

inline Real2 operator *(Real c, Real2 b) {
  return Real2(c * b.x(), c * b.y());
}

class Complex {
  Real m_real;
  Real m_imag;

public:
  Complex() :
      m_real(0.0), m_imag(0.0) {
  }
  Complex(Real real, Real imag) :
      m_real(real), m_imag(imag) {
  }

  Complex(const Complex &v) :
      m_real(v.m_real), m_imag(v.m_imag) {
  }

  Complex(Real real) : m_real(real), m_imag(0.0)
      {}

  Real abs2() const {
    return m_real*m_real + m_imag*m_imag;
  }

  Real abs() const {
    return std::sqrt(abs2());
  }

  Complex & operator +=(Complex b) {
    m_real += b.m_real;
    m_imag += b.m_imag;
    return *this;
  }
  Complex & operator -=(Complex b) {
    m_real -= b.m_real;
    m_imag -= b.m_imag;
    return *this;
  }

  Complex operator +( Complex b) const {
    return Complex(m_real + b.m_real, m_imag + b.m_imag);
  }

  Complex operator -( Complex b) const {
    return Complex(m_real - b.m_real, m_imag - b.m_imag);
  }

  Complex operator -() const {
    return Complex(-m_real, -m_imag);
  }
  bool operator ==(const Complex b) const {
    return m_real == b.m_real && m_imag == b.m_imag;
  }
  bool operator !=(const Complex b) const {
    return !operator ==(b);
  }

  Complex operator *(Real b) const {
    return Complex(m_real * b, m_imag * b);
  }

  Complex operator /(Real b) const {
    return Complex(m_real * (1 / b), m_imag * (1 / b));
  }
  Complex operator *(const Complex b) const {
    return Complex(m_real * b.m_real - m_imag * b.m_imag,
	m_real * b.m_imag + m_imag * b.m_real);
  }

  Complex conjug() const {
    return Complex(m_real, -m_imag);
  }

  Real realPart() const {
    return m_real;
  }
  Real imagPart() const {
    return m_imag;
  }
};

}

namespace math {
// Vector Product
using Smooth3D::Real;
using Smooth3D::Real2;
using Smooth3D::Real3;

inline Real abs(const Real & u) {
  return std::fabs(u);
}

inline Real sqr(const Real & u) {
  return (u*u);
}

inline Real3 vecMul(Real3 v1, Real3 v2) {
  return Real3(v1.y() * v2.z() - v1.z() * v2.y(), -v1.x() * v2.z() + v1.z() * v2.x(),
      v1.x() * v2.y() - v1.y() * v2.x());
}
inline Real scaMul(Real3 u, Real3 v) {
  return (u.x() * v.x() + u.y() * v.y() + u.z() * v.z());
}

inline Real scaMul(Real2 u, Real2 v) {
  return (u.x() * v.x() + u.y() * v.y());
}
inline Real mixteMul(Real3 v1, Real3 v2, Real3 v3) {
  return v1.x() * v2.y() * v3.z() + v2.x() * v3.y() * v1.z() + v3.x() * v1.y() * v2.z()
      - v1.z() * v2.y() * v3.x() - v2.z() * v3.y() * v1.x() - v3.z() * v1.y() * v2.x();
}
}

#endif
