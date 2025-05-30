//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation for the
// RigidFrameTransf class. RigidFrameTransf is a nonlinear
// transformation for a space frame between the global
// and basic coordinate systems
//
// Written: cmp
// Created: 04/2025
//
#pragma once
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>
#include <Node.h>
#include <Channel.h>
#include <Logging.h>
#include <Rotations.hpp>
#include <RigidFrameTransf.hpp>

using namespace OpenSees;


template <int nn, int ndf, typename BasisT>
RigidFrameTransf<nn,ndf,BasisT>::RigidFrameTransf(int tag, 
                                           const Vector3D &vecxz, 
                                           const std::array<Vector3D, nn> *offset,
                                           int offset_flags)
  : FrameTransform<nn,ndf>(tag),
    L(0),
    nodes{},
    ur{},
    offsets{nullptr},
    offset_flags(offset_flags),
    basis{nodes, vecxz}
{

  double nz = vecxz.norm();
  for (int i=0; i<3; i++)
    vz[i] = vecxz[i]/nz;

  // Rigid joint offsets
  if (offset != nullptr) {
    offsets = new std::array<Vector3D, nn>{};
    *offsets = *offset;
  }
}



template <int nn, int ndf, typename BasisT>
RigidFrameTransf<nn,ndf,BasisT>::~RigidFrameTransf()
{
  if (offsets != nullptr)
    delete offsets;
}

template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::commit()
{
  return 0;
}

template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::revertToLastCommit()
{
  return 0;
}

template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::revertToStart()
{
  return 0;
}


template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::initialize(std::array<Node*, nn>& new_nodes)
{

  for (int i=0; i<nn; i++) {
    nodes[i] = new_nodes[i];
    if (nodes[i] == nullptr) {
      opserr << "invalid pointers to the element nodes\n";
      return -1;
    }
  }

  int error;
  // set element length and orientation
  if ((error = this->computeElemtLengthAndOrient()))
    return error;

  return 0;
}


template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::computeElemtLengthAndOrient()
{

  const Vector &XI = nodes[   0]->getCrds();
  const Vector &XJ = nodes[nn-1]->getCrds();

  for (int i=0; i<3; i++) {
    xi[i] = XI[i];
    xj[i] = XJ[i];
  }
  
  Vector3D dx = xj - xi;

  if (offsets != nullptr) {
    for (int i=0; i<3; i++)
      dx(i) -= (*offsets)[   0][i];
    for (int i=0; i<3; i++)
      dx(i) += (*offsets)[nn-1][i];
  }

  // calculate the element length
  L = dx.norm();

  if (L == 0.0)
    return -2;

  return basis.initialize();
}


template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::getLocalAxes(Vector3D &e1, Vector3D &e2, Vector3D &e3) const
{
  Matrix3D R = basis.getRotation();
  for (int i = 0; i < 3; i++) {
    e1[i] = R(i,0);
    e2[i] = R(i,1);
    e3[i] = R(i,2);
  }
  return 0;
}

template <int nn, int ndf, typename BasisT>
double
RigidFrameTransf<nn,ndf,BasisT>::getInitialLength()
{
  return L;
}

template <int nn, int ndf, typename BasisT>
double
RigidFrameTransf<nn,ndf,BasisT>::getDeformedLength()
{
  return L;
}


//
// Pull
//
template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::update()
{
  if (basis.update() < 0) 
    return -1;

  Matrix3D R = basis.getRotation();
  for (int i=0; i<nn; i++) {
    Versor q = nodes[i]->getTrialRotation();
    ur[i] = LogSO3(R^MatrixFromVersor(q));
  }

  return 0;
}


template <int nn, int ndf, typename BasisT>
VectorND<nn*ndf> 
RigidFrameTransf<nn,ndf,BasisT>::pullVariation(const VectorND<nn*ndf>& ug, 
             const Matrix3D& R, 
             const std::array<Vector3D, nn> *offset,
             int offset_flags) 
{

  constexpr static int N = nn * ndf;

  // Initialize ul = ug
  VectorND<N> ul = ug;

  // (1)
  // Do ui -= ri x wi
  if constexpr (ndf >= 6)
    if (offset && !(offset_flags&OffsetLocal)) [[unlikely]] {
      const std::array<Vector3D, nn>& offsets = *offset;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, offsets[i].cross(w), -1.0);
      }
    }

  // (2) Rotations and translations

  for (int i=0; i<nn; i++) {
    const int j = i * ndf;
    ul.insert(j  , R^Vector3D{ul[j+0], ul[j+1], ul[j+2]}, 1.0);
    ul.insert(j+3, R^Vector3D{ul[j+3], ul[j+4], ul[j+5]}, 1.0);
  }

  double Ln = basis.getLength();
  {
    Vector3D wr = basis.getRotationVariation(ndf, &ul[0]);
    Vector3D dc = basis.getPositionVariation(ndf, &ul[0]);
    // Vector3D c  = basis.getPosition();

    for (int i=0; i<nn; i++) {

      #if 1
      Vector3D ui = this->getNodePosition(i);
      // ui -= c;
      ul.assemble(i*ndf+0, dc, -1.0);
      ul.assemble(i*ndf+0, ui.cross(wr), 1.0);
      #else
      int j = i * ndf;
      Vector3D xi = {double(i)/double(nn-1)*Ln, 0, 0}; // R^(nodes[i]->getCrds()); // 
      // xi += ui;
      // xi -= c;
      opserr << "u[" << i << "] = " << Vector(ul.template extract<3>(j));

      ul.assemble(i*ndf+0, dc, -1.0);
      opserr << "u[" << i << "] = " << Vector(ul.template extract<3>(j));
      ul.assemble(i*ndf+0,  c.cross(wr), -1.0);
      opserr << "u[" << i << "] = " << Vector(ul.template extract<3>(j));
      ul.assemble(i*ndf+0, DR^(nodes[i]->getCrds()), 1.0);
      opserr << "u[" << i << "] = " << Vector(ul.template extract<3>(j));
      ul.assemble(i*ndf+0, xi.cross(wr), 1.0);
      #endif
      ul.assemble(i*ndf+3, wr, -1.0);
    }
  }

  // 3) Offsets
  if constexpr (ndf >= 6)
    if (offset && (offset_flags&OffsetLocal)) [[unlikely]] {
      const std::array<Vector3D, nn>& offsets = *offset;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, offsets[i].cross(w), -1.0);
      }
    }

  // (5) Logarithm of rotations
  if (0) { // !(offset_flags & LogIter)) {
    for (int i=0; i<nn; i++) {
      const int j = i * ndf+3;
      Vector3D v {ul[j+0], ul[j+1], ul[j+2]};
      ul.insert(i*ndf+3, dLogSO3(ur[i])*v, 1.0);
    }
  }

  return ul;
}

template <int nn, int ndf, typename BasisT>
VectorND<nn*ndf>
RigidFrameTransf<nn,ndf,BasisT>::getStateVariation()
{

  static VectorND<nn*ndf> ug;
  for (int i=0; i<nn; i++) {
    const Vector &ddu = nodes[i]->getIncrDeltaDisp();
    for (int j = 0; j < ndf; j++) {
      ug[i*ndf+j] = ddu(j);
    }
  }

  Matrix3D R = basis.getRotation();
  return RigidFrameTransf<nn,ndf,BasisT>::pullVariation(ug, R, offsets, offset_flags);
}

template <int nn, int ndf, typename BasisT>
Vector3D
RigidFrameTransf<nn,ndf,BasisT>::getNodePosition(int node)
{
  Vector3D v = this->pullPosition<&Node::getTrialDisp>(node) 
             - basis.getPosition();

  v += basis.getRotationDelta()^(nodes[node]->getCrds());

  return v;
}


template <int nn, int ndf, typename BasisT>
Vector3D
RigidFrameTransf<nn,ndf,BasisT>::getNodeRotationLogarithm(int node)
{
  return ur[node];
}


//
// Push
//
template <int nn, int ndf, typename BasisT>
VectorND<nn*ndf>
RigidFrameTransf<nn,ndf,BasisT>::pushResponse(VectorND<nn*ndf>&p)
{
  VectorND<nn*ndf> pa = p;

  // 1) Logarithm
  if (0) { // !(offset_flags & LogIter)) {
    for (int i=0; i<nn; i++) {
      const int j = i * ndf+3;
      Vector3D m {p[j+0], p[j+1], p[j+2]};
      pa.insert(j, dLogSO3(ur[i])^m, 1.0);
    }
  }

  MatrixND<nn*ndf,nn*ndf> A = getProjection();
  pa = A^pa;

  // 3,4) Rotate and joint offsets
  auto pg = this->FrameTransform<nn,ndf>::pushConstant(pa);

  return pg;
}


template <int nn, int ndf, typename BasisT>
MatrixND<nn*ndf,nn*ndf>
RigidFrameTransf<nn,ndf,BasisT>::pushResponse(MatrixND<nn*ndf,nn*ndf>&kb, const VectorND<nn*ndf>&pb)
{
  MatrixND<nn*ndf,nn*ndf> Kb = kb;
  VectorND<nn*ndf> p = pb;

  if (0) {//!(offset_flags & LogIter)) {
    for (int i=0; i<nn; i++) {
      Vector3D m{pb[i*ndf+3], pb[i*ndf+4], pb[i*ndf+5]};
      const Matrix3D Ai = dLogSO3(ur[i]);
      p.insert(i*ndf+3, Ai^m, 1.0);

      Matrix3D kg = ddLogSO3(ur[i], m);
      for (int j=0; j<nn; j++) {
        const Matrix3D Aj = dLogSO3(ur[j]);
        // loop over 3x3 blocks for n and m
        for (int k=0; k<2; k++) {
          for (int l=0; l<2; l++) {
            Matrix3D Kab {{
              {Kb(i*ndf+3*k+0, j*ndf+3*l  ), Kb(i*ndf+3*k+1, j*ndf+3*l  ), Kb(i*ndf+3*k+2, j*ndf+3*l  )},
              {Kb(i*ndf+3*k+0, j*ndf+3*l+1), Kb(i*ndf+3*k+1, j*ndf+3*l+1), Kb(i*ndf+3*k+2, j*ndf+3*l+1)},
              {Kb(i*ndf+3*k+0, j*ndf+3*l+2), Kb(i*ndf+3*k+1, j*ndf+3*l+2), Kb(i*ndf+3*k+2, j*ndf+3*l+2)}
            }};
            if (k == 1)
              Kab = Kab*Aj;
            if (l == 1)
              Kab = Ai^Kab;

            Kb.insert(Kab, i*ndf+3*k, j*ndf+3*l, 1.0);
            if (i == j && k == 1 && l == 1)
              Kb.assemble(kg, i*ndf+3*k, j*ndf+3*l, 1.0);
          }
        }
      }
    }
  }

  // Kb = kb;

  MatrixND<nn*ndf,nn*ndf> Kl;
  MatrixND<nn*ndf,nn*ndf> A = getProjection();
  Kl.addMatrixTripleProduct(0, A, Kb, 1);

  //
  // Kl += -W'*Pn'*A
  //
  p = A^p;
  Kb.zero();
  for (int j=0; j<nn; j++) {
    MatrixND<3,6> Gj = basis.getRotationGradient(j);
    for (int i=0; i<nn; i++) {
      auto PnGj = Hat(&p[i*ndf+0])*Gj;
      Kb.assemble(PnGj,                i*ndf+0, j*ndf, -1.0);

      // Kl += -Pnm*W
      Kl.assemble(PnGj,                i*ndf+0, j*ndf, -1.0);
      Kl.assemble(Hat(&p[i*ndf+3])*Gj, i*ndf+3, j*ndf, -1.0);
    }
  }
  Kl.addMatrixTransposeProduct(1.0, Kb, A,  -1.0);

  // Kl = diag(R) * Kl * diag(R)^T
  return this->FrameTransform<nn,ndf>::pushConstant(Kl);
}


template <int nn, int ndf, typename BasisT>
FrameTransform<nn,ndf> *
RigidFrameTransf<nn,ndf,BasisT>::getCopy() const
{
  return new RigidFrameTransf<nn,ndf,BasisT>(this->getTag(), vz, offsets);
}


//
// Sensitivity
//
template <int nn, int ndf, typename BasisT>
bool
RigidFrameTransf<nn,ndf,BasisT>::isShapeSensitivity()
{
  int nodeParameterI = nodes[   0]->getCrdsSensitivity();
  int nodeParameterJ = nodes[nn-1]->getCrdsSensitivity();
  // TODO(sensitivity): implement dvz

  return (nodeParameterI != 0 || nodeParameterJ != 0);
}


template <int nn, int ndf, typename BasisT>
double
RigidFrameTransf<nn,ndf,BasisT>::getLengthGrad()
{
  const int di = nodes[0]->getCrdsSensitivity();
  const int dj = nodes[1]->getCrdsSensitivity();

  Vector3D dxi{0.0};
  Vector3D dxj{0.0};

  if (di != 0)
    dxi(di-1) = 1.0;
  if (dj != 0)
    dxj(dj-1) = 1.0;

  return 1/L*(xj - xi).dot(dxj - dxi);
}

template <int nn, int ndf, typename BasisT>
double
RigidFrameTransf<nn,ndf,BasisT>::getd1overLdh()
{
  return -getLengthGrad()/(L*L);
}


template <int nn, int ndf, typename BasisT>
void
RigidFrameTransf<nn,ndf,BasisT>::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"RigidFrameTransf\"";
    s << ", \"vecxz\": [" 
      << vz[0] << ", " 
      << vz[1] << ", "
      << vz[2] << "]";
    if (offsets != nullptr) {
      s << ", \"offsets\": [";
      for (int i=0; i<nn; i++) {
        s << "["
          << (*offsets)[i][0] << ", " 
          << (*offsets)[i][1] << ", "
          << (*offsets)[i][2] << "]";
        if (i < nn-1)
          s << ", ";
      }
      s << "]";
    }

    s << "}";

    return;
  }
}

