//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// PDeltaFrameTransf.h. PDeltaFrameTransf provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems
//
// Adapted: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
//
#ifndef PDeltaFrameTransf_h
#define PDeltaFrameTransf_h

#include <array>
#include <FrameTransform.h>
#include <Vector.h>
#include <Matrix.h>

template <int nn, int ndf>
class PDeltaFrameTransf: public FrameTransform<nn,ndf>
{
public:

    PDeltaFrameTransf(int tag, 
                      const Vector3D &vecxz,
                      const std::array<Vector3D, nn> *offset=nullptr,
                      int offset_flags = 0);

    ~PDeltaFrameTransf();
    
    const char *getClassType() const {return "PDeltaFrameTransf";}
    
    double getInitialLength();
    double getDeformedLength();

    virtual int initialize(std::array<Node*, nn>& new_nodes);
    virtual int update();
    virtual int commit();
    virtual int revertToLastCommit();        
    virtual int revertToStart();

    virtual VectorND<nn*ndf> getStateVariation() final;
    virtual Vector3D getNodePosition(int tag) final;
    virtual Vector3D getNodeRotationLogarithm(int tag) final;
    virtual const std::array<Vector3D,nn> *getRigidOffsets() const { return linear.getRigidOffsets();}

    virtual VectorND<nn*ndf>    pushResponse(VectorND<nn*ndf>&pl) final;
    virtual MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl) final;

    virtual FrameTransform<nn,ndf> *getCopy() const;

    virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const;

    // Sensitivity
    double getLengthGrad();

    // Tagged Object
    void Print(OPS_Stream &s, int flag = 0);

    private:
      int offset_flags;
      LinearFrameTransf<nn,ndf> linear;

};
#include "PDeltaFrameTransf3d.tpp"
#endif
