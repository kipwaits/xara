/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Description: This file contains the class definition for
// CorotCrdTransf2d.h. CorotCrdTransf2d provides the
// abstraction of a corotation transformation for a planar frame element
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
//
#ifndef CorotCrdTransf2d_h
#define CorotCrdTransf2d_h

#include <FrameTransform.h>
#include <Vector.h>
#include <Matrix.h>

class CorotCrdTransf2d: public CrdTransf
{
public:
    CorotCrdTransf2d(int tag, const Vector &rigJntOffsetI, const Vector &rigJntOffsetJ);
    
    CorotCrdTransf2d();
    ~CorotCrdTransf2d();
    
    const char *getClassType() const {
        return "CorotCrdTransf2d";
    }
    
    int initialize(Node *nodeIPointer, Node *nodeJPointer);
    int update();
    double getInitialLength();
    double getDeformedLength();
    
    int commitState();
    int revertToLastCommit();        
    int revertToStart();
    
    const Vector &getBasicTrialDisp();
    const Vector &getBasicIncrDisp();
    const Vector &getBasicIncrDeltaDisp();
    const Vector &getBasicTrialVel();
    const Vector &getBasicTrialAccel();
    
    const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &p0);
    const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce);
    const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff);
    
    // Sensitivity
    const Vector &getBasicDisplFixedGrad();
    const Vector &getBasicDisplTotalGrad(int gradNumber);
    const Vector &getGlobalResistingForceShapeSensitivity(const Vector &q,
							  const Vector &p0,
							  int gradNumber);
    bool isShapeSensitivity();
    double getLengthGrad();
    double getd1overLdh();

    CrdTransf *getCopy2d();
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);

    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps);    

    int getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);
  int getRigidOffsets(Vector &offsets);
  
private:
    int compElemtLengthAndOrient();
    int compElemtLengthAndOrientWRTLocalSystem(const Vector &ul);
    void compTransfMatrixLocalGlobal(Matrix &Tlg);
    void compTransfMatrixBasicLocal(Matrix &Tbl);
    void transfLocalDisplsToBasic(const Vector &ul);
    const Matrix &getGeomStiffMatrix(const Vector &pb) const;
    
    // internal data
    Node *nodeIPtr, *nodeJPtr; // pointers to the element two endnodes
    
    Vector nodeIOffset, nodeJOffset;  // rigid joint offsets
    
    double cosTheta, sinTheta; // direction cosines of undeformed element wrt to global system 
    double cosAlpha, sinAlpha; // direction cosines of deformed element wrt to local system
    double L;                  // undeformed element length
    double Ln;                 // deformed element length
    double Lx, Ly;             // components of the deformed member
    double Lxdot, Lydot;       // time derivatives of components of the deformed member
    double Lxdotdot, Lydotdot; // double time derivatives of components of the deformed member
    
    Vector ub;                 // basic displacements
    Vector ubcommit;           // committed basic displacements
    Vector ubpr;               // previous basic displacements
    
    static Matrix Tlg;         // matrix that transforms from global to local coordinates
    static Matrix Tbl;         // matrix that transforms from local  to basic coordinates
    static Matrix kg;          // global stiffness matrix
    static Vector uxg;     
    static Vector pg;     
    static Vector dub;     
    static Vector Dub;     
    
    double *nodeIInitialDisp, *nodeJInitialDisp;
    bool initialDispChecked;
    bool nodeOffsets;
};
#endif
