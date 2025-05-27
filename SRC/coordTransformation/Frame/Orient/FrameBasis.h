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
    virtual int       initialize() =0;
    virtual int       update() =0;
    virtual Matrix3D  getRotation() =0;
    virtual Vector3D  getPosition() =0;

    virtual Vector3D      getRotationVariation() =0;
    virtual MatrixND<3,6> getRotationGradient(int node) =0;

    // virtual int getBasis(Vector3D&, Vector3D&, Vector3D&) =0;
};


template <int nn>
class RankineBasis : public FrameBasis
{
public:
    RankineBasis(std::array<Node*,nn>& nodes, const Vector3D& vecxz)
    : nodes{nodes}, vz(vecxz) {
    };

    virtual int 
    initialize() {

      const Vector &XI = nodes[   0]->getCrds();
      const Vector &XJ = nodes[nn-1]->getCrds();
    
      Vector3D xi, xj;
      for (int i=0; i<3; i++) {
        xi[i] = XI[i];
        xj[i] = XJ[i];
      }
      
      Vector3D dx = xj - xi;

      Vector3D e1 = dx/dx.norm();
  
      //
      Vector3D e2 = vz.cross(e1);
  
      const double ynorm = e2.norm();
  
      if (ynorm == 0.0)
          return -1;
  
      e2 /= ynorm;
  
      Vector3D e3 = e1.cross(e2);
  
      for (int i = 0; i < 3; i++) {
        R[init](i,0) = e1[i];
        R[init](i,1) = e2[i];
        R[init](i,2) = e3[i];
      }


      //   int ic = 0; // std::floor(0.5*(nn+1));
      c[init] = {0,0,0}; // nodes[ic]->getCrds();

      update();
      return 0;
    }

    virtual int
    update() {
      int ic = 0; // std::floor(0.5*(nn+1));
      R[pres] = R[init] * MatrixFromVersor(nodes[ic]->getTrialRotation());
      Vector3D uc = nodes[ic]->getTrialDisp();
      c[pres] = c[init] - (R[pres]^uc);
      return 0;
    };

    virtual Vector3D
    getRotationVariation() {
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
    getRotation() {
        return R[pres];
    }

    virtual Vector3D
    getPosition() {
        return c[pres];
    }

private:
  enum { pres, prev, init};
  double L; // TODO
  Vector3D vz;
  Matrix3D R[3];
  Vector3D c[3];
  Matrix3D dR;
  std::array<Node*,nn>& nodes;
};

} // namespace OpenSees
