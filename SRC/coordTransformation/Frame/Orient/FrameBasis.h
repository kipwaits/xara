//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
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

  virtual double        getLength() const =0;
  // x, \Lambda
  virtual Matrix3D      getRotation() const =0;
  virtual Vector3D      getPosition() =0;
  // \psi
  virtual Vector3D      getPositionVariation(int ndf, double* du) =0; 
  virtual Vector3D      getRotationVariation(int ndf, double* du) =0;
  virtual Matrix3D      getRotationDelta() =0;
  //
  virtual MatrixND<3,6> getRotationGradient(int node) =0;

};


template <int nn>
class RankineBasis : public FrameBasis
{
public:
  RankineBasis(std::array<Node*,nn>& nodes, const Vector3D& vecxz)
  : nodes{nodes}, vz(vecxz), Xc{}, c{}, R{} {
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
    Xc = nodes[ic]->getCrds();
    c[init] = R[init]^(Xc);
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

    Vector3D uc = nodes[ic]->getTrialDisp();
    Vector3D X = nodes[ic]->getCrds(); // R[init]*c[init];
    c[pres] = R[pres]^(X + uc);
    return 0;
  };

  virtual double 
  getLength() const override {
    return Ln;
  }


  virtual Vector3D 
  getPositionVariation(int ndf, double* du) {
    return Vector3D {du[ndf*ic+0], du[ndf*ic+1], du[ndf*ic+2]};
  }

  virtual Vector3D
  getRotationVariation(int ndf, double* du) {
    // psi_r = omega
    Vector3D w{};
    for (int i=0; i<nn; i++) {
      // const Vector &du = nodes[i]->getIncrDeltaDisp();
      auto Wi = this->getRotationGradient(i);
      for (int j=0; j<3; j++)
        for (int k=0; k<6; k++)
          w[j] += Wi(j,k) * du[ndf*i + k];
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
      Gb.template insert<0,0>(ix, -1/Ln);
    else if (node == nn-1)
      Gb.template insert<0,0>(ix,  1/Ln);
  
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

  Vector3D
  getLocation() {
    return c[pres];
  }

  virtual Vector3D
  getPosition() {
    // Return Delta c
    Vector3D X = nodes[ic]->getCrds();
    Vector3D Dc =  c[pres] - (R[init]^X) ; // (R[pres]^c[init]);
    return Dc;
  }

private:
  constexpr static int ic = 0; // std::floor(0.5*(nn+1));
  enum { pres, prev, init};
  double L, Ln;
  Vector3D vz, dX, Xc;
  Matrix3D R[3];
  Vector3D c[3];
  Matrix3D dR;
  std::array<Node*,nn>& nodes;
};

} // namespace OpenSees
