#ifndef LAGRANGE_MAPPING_H
#define LAGRANGE_MAPPING_H

namespace Smooth3D {

// Tentative de rationnaliser les definitions isoparametriques des elements

class LagrangeMapping2D {
protected:
  LagrangeMapping2D() {
  }
public:
  virtual ~ LagrangeMapping2D() {
  }

  virtual Real3 eval(const Real2 &) const = 0;
  virtual void gradient(const Real2 &, Real3 &, Real3 &) const = 0;
  virtual bool isvalid(const Real2 &) const = 0;

};

class TriangleLagrangeMapping2D: public LagrangeMapping2D {
  Real3 m_a[3];
public:
  TriangleLagrangeMapping2D(const Real3 & m_a0, const Real3 & m_a1,
      const Real3 & m_a2) {
    m_a[0] = m_a0;
    m_a[1] = m_a1;
    m_a[2] = m_a2;
  }

  virtual Real3 eval(const Real2 & uv) const {
    Real u = uv.x(), v = uv.y();
    return (1.0 - u - v) * m_a[0] + u * m_a[1] + v * m_a[2];
  }

  virtual void gradient(const Real2 & uv, Real3 & dmdu, Real3 & dmdv) const {
    dmdu = m_a[1] - m_a[0];
    dmdv = m_a[2] - m_a[0];
  }

  virtual bool isvalid(const Real2 & uv) const {
    Real u = uv.x(), v = uv.y();
    return u >= 0.0 && v >= 0.0 && u + v <= 1.0;
  }
};

class QuadrangleLagrangeMapping2D: public LagrangeMapping2D {
  Real3 m_a[4];
public:
  QuadrangleLagrangeMapping2D(const Real3 & m_a0, const Real3 & m_a1,
      const Real3 & m_a2, const Real3 & m_a3) {
    m_a[0] = m_a0;
    m_a[1] = m_a1;
    m_a[2] = m_a2;
    m_a[3] = m_a3;
  }

  virtual Real3 eval(const Real2 & uv) const {
    Real u = uv.x(), v = uv.y();
    return (1.0 - u) * ((1.0 - v) * m_a[0] + v * m_a[2])
	+ u * ((1.0 - v) * m_a[1] + v * m_a[3]);
  }

  virtual void gradient(const Real2 & uv, Real3 & dmdu, Real3 & dmdv) const {

    Real u = uv.x(), v = uv.y();
    dmdu = (1.0 - v) * (m_a[1] - m_a[0]) + v * (m_a[3] - m_a[2]);
    dmdv = (1.0 - u) * (m_a[2] - m_a[0]) + u * (m_a[3] - m_a[1]);
  }

  virtual bool isvalid(const Real2 & uv) const {
    Real u = uv.x(), v = uv.y();
    return u >= 0.0 && v >= 0.0 && u <= 1.0 && v <= 1.0;
  }
};

}

#endif
