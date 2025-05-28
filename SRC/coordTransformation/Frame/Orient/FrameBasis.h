#pragma once
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <Versor.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <MatrixND.h>
#include <Matrix3D.h>

namespace OpenSees {
class FrameBasis
{
public:
  virtual int           initialize() =0;
  virtual int           update() =0;

  virtual double getLength() const =0;
  // x, \Lambda
  virtual Matrix3D      getRotation() const =0;
  virtual Vector3D      getPosition() =0;
  // \psi
  virtual Vector3D      getPositionVariation() =0; 
  virtual Vector3D      getRotationVariation() =0;
  virtual Matrix3D      getRotationDelta() =0;
  //
  virtual MatrixND<3,6> getRotationGradient(int node) =0;

  // virtual int getBasis(Vector3D&, Vector3D&, Vector3D&) =0;
};


template <int nn>
class RankineBasis : public FrameBasis
{
public:
  RankineBasis(std::array<Node*,nn>& nodes, const Vector3D& vecxz)
  : nodes{nodes}, vz(vecxz), c{}, R{} {
  };

  virtual int 
  initialize() {

    for (int i=0; i<nn; i++)
      nodes[i]->getTrialRotation();

    const Vector &XI = nodes[   0]->getCrds();
    const Vector &XJ = nodes[nn-1]->getCrds();
  

    for (int i=0; i<3; i++)
      dX[i] = XJ[i] - XI[i];

    L = dX.norm();
    Ln = L;
    Vector3D e1 = dX/L;

    //
    Vector3D e2 = vz.cross(e1);

    const double ynorm = e2.norm();

    if (ynorm == 0.0)
        return -1;

    e2 /= ynorm;

    Vector3D e3 = e1.cross(e2);

    e2 = e3.cross(e1);

    for (int i = 0; i < 3; i++) {
      R[init](i,0) = e1[i];
      R[init](i,1) = e2[i];
      R[init](i,2) = e3[i];
    }
#if 1
    c[init] = R[init]^(nodes[ic]->getCrds());
#endif
    update();
    return 0;
  }

  virtual int
  update() {

    Vector3D e1 = dX;
    {
      //
      // Update state
      //
      {
        // Relative translation
#if 1
        const Vector& uI = nodes[   0]->getTrialDisp();
        const Vector& uJ = nodes[nn-1]->getTrialDisp();
        for (int k = 0; k < 3; k++)
          e1[k] += uJ(k) - uI(k);
#endif
        // Calculate the deformed length
        Ln = e1.norm();

        if (Ln == 0.0) [[unlikely]] {
          opserr << "\nSouzaFrameTransf: deformed length is 0.0\n";
          return -2;
        }

        e1 /= Ln;
      }
    }

    {
#if 1
      Matrix3D Ri = MatrixFromVersor(nodes[0]->getTrialRotation()); //*R[init];
      Ri *= 0.5;
      Ri.addMatrix(MatrixFromVersor(nodes[0]->getTrialRotation()), 0.5);
      Vector3D e2 = Ri^(vz.cross(e1));
      // Vector3D e2 = (R[pres]^Ri)*(vz.cross(e1));
      // Vector3D e2 = Ri*R[pres]^(vz.cross(e1));
#else 
      Vector3D e2 = vz.cross(e1);
#endif
      Vector3D e3 = e1.cross(e2);
      e3 /= e3.norm();

      e2 = e3.cross(e1);

      for (int i = 0; i < 3; i++) {
        R[pres](i,0) = e1[i];
        R[pres](i,1) = e2[i];
        R[pres](i,2) = e3[i];
      }
    }

#if 0
    Vector3D uc = nodes[ic]->getTrialDisp();
#else 
    Vector3D uc{};
#endif
    c[pres] = R[pres]^(R[init]*c[init] + uc);
    return 0;
  };

  virtual double 
  getLength() const override {
    return Ln;
  }

  // virtual Vector3D
  // getMoment() {
  // }

  virtual Vector3D 
  getPositionVariation() {
    // psi_x = c[pres] - ();
    Vector3D psi_x{};
    return psi_x;
  }

  virtual Vector3D
  getRotationVariation() {
    // psi_r = omega
    Vector3D w{};
    for (int i=0; i<nn; i++) {
      const Vector &du = nodes[i]->getIncrDeltaDisp();
      w += this->getRotationGradient(i) * du;
    }
    return w;
  }


  virtual MatrixND<3,6> 
  getRotationGradient(int node) {
    MatrixND<3,6> Gb{};

    constexpr Vector3D axis{1, 0, 0};
    constexpr Matrix3D ix = Hat(axis);
    constexpr Matrix3D ioi = axis.bun(axis);
    Gb.template insert<0, 3>(ioi, 0.5);
    if (node == 0)
      Gb.template insert<0,0>(ix, -1/L);
    else if (node == nn-1)
      Gb.template insert<0,0>(ix,  1/L);
  
    return Gb;
  }


  virtual Matrix3D
  getRotation() const override {
    return R[pres];
  }


  virtual Matrix3D 
  getRotationDelta() {
    return R[pres] - R[init];
  }


  virtual Vector3D
  getPosition() {
    // Return Delta c
    Vector3D Dc =  c[pres] - (R[pres]^c[init]);
    return Dc;
  }

private:
  constexpr static int ic = 0; // std::floor(0.5*(nn+1));
  enum { pres, prev, init};
  double L, Ln;
  Vector3D vz, dX;
  Matrix3D R[3];
  Vector3D c[3];
  Matrix3D dR;
  std::array<Node*,nn>& nodes;
};

} // namespace OpenSees
