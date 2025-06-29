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
// Description: This file contains the implementation for FullGenLinSOE
//
// Written: fmk 
// Created: February 1997
// Revision: A
//
#include <FullGenLinSOE.h>
#include <FullGenLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <assert.h>
#include <math.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <iostream>

FullGenLinSOE::FullGenLinSOE(FullGenLinSolver &theSolvr)
: LinearSOE(theSolvr, LinSOE_TAGS_FullGenLinSOE),
 size(0), A(0), B(0), X(0), 
 matA(nullptr),
 factored(false)
{
    theSolvr.setLinearSOE(*this);
}


FullGenLinSOE::FullGenLinSOE(int N, FullGenLinSolver &theSolvr)
:LinearSOE(theSolvr, LinSOE_TAGS_FullGenLinSOE),
 size(0), A(0), B(N), X(N), 
 matA(nullptr),
 factored(false)
{
    size = N;
    A = new double[size*size]{};
    matA  = new Matrix(A, size, size);

    theSolvr.setLinearSOE(*this);
    
    // invoke setSize() on the Solver        
    if (theSolvr.setSize() < 0) {
	// opserr << "WARNING :FullGenLinSOE::FullGenLinSOE :";
	// opserr << " solver failed setSize() in constructor\n";
    }    
    
}

    
FullGenLinSOE::~FullGenLinSOE()
{
  if (A != 0) delete [] A;     
  if (matA != 0) delete matA;        
}


int
FullGenLinSOE::getNumEqn() const
{
  return size;
}

int 
FullGenLinSOE::setSize(Graph &theGraph)
{
    int result = 0;
    int oldSize = size;
    size = theGraph.getNumVertex();

    if (size > oldSize) { // we have to get another space for A

        if (A != nullptr) 
            delete [] A;

        A = new double[size*size];
    }

    factored = false;

    B.resize(size);
    B.Zero();
    X.resize(size);
    X.Zero();


    // create new Vectors
    if (size != oldSize || matA == nullptr) {
        if (matA != nullptr)
            delete matA;
        matA  = new Matrix(A, size, size);	
    }
    matA->Zero();

    // invoke setSize() on the Solver    
    LinearSOESolver *theSolvr = this->getSolver();
    int solverOK = theSolvr->setSize();
    if (solverOK < 0) {
        // opserr << "WARNING:FullGenLinSOE::setSize :";
        // opserr << " solver failed setSize()\n";
        return solverOK;
    }    
    
    return result;
}

int 
FullGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check that m and id are of similar size
    assert(id.Size() == m.noRows() && id.Size() == m.noCols());

    // check for a quick return 
    if (fact == 0.0)
      return 0;

    int idSize = id.Size();

    if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *startColiPtr = A + col*size;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0) {
			 double *APtr = startColiPtr + row;
			 *APtr += m(j,i);
		    }
		}
	    } 
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *startColiPtr = A + col*size;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0) {
			 double *APtr = startColiPtr + row;
			 *APtr += m(j,i) * fact;
		     }
		}  // for j
	    } 
	}  // for i
    }    
    return 0;
}



int 
FullGenLinSOE::addColA(const Vector &colData, int col, double fact)
{
  
  assert(colData.Size() == size);
  assert(col <= size && col >= 0);

  if (fact == 0.0)
    return 0; 

  if (fact == 1.0) { // do not need to multiply
    double *coliPtr = A + col*size;
    for (int row=0; row<size; row++) {
      *coliPtr += colData(row);
      coliPtr++;
    }

  } else {
    double *coliPtr = A + col*size;
    for (int row=0; row<size; row++) {
      *coliPtr += colData(row) * fact;
      coliPtr++;
    }

  }

  return 0;
}


int 
FullGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar size
    assert(idSize == v.Size() );

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	
    return 0;
}



int
FullGenLinSOE::setB(const Vector &v, double fact)
{
    assert (v.Size() == size);

    // check for a quick return 
    if (fact == 0.0)
      return 0;

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
        for (int i=0; i<size; i++) {
            B[i] = v(i);
        }
    } else if (fact == -1.0) {
        for (int i=0; i<size; i++) {
            B[i] = -v(i);
        }
    } else {
        for (int i=0; i<size; i++) {
            B[i] = v(i) * fact;
        }
    }	
    return 0;
}

void 
FullGenLinSOE::zeroA()
{
    double *Aptr = A;
    int theSize = size*size;
    for (int i=0; i<theSize; i++)
        *Aptr++ = 0;

    factored = false;
}
	
void 
FullGenLinSOE::zeroB()
{
    B.Zero();
}

int
FullGenLinSOE::formAp(const Vector &p, Vector &Ap)
{
  // Check that p and A are same size
  assert(size == p.Size() && size == Ap.Size() && p.Size() == Ap.Size());

  for (int row = 0; row < size; row++) {
    double sum = 0.0;
    double *APtr = A + row;
    for (int col = 0; col < size; col++) {
      APtr += size;
      sum += (*APtr) * p(col);
    }
    Ap(row) = sum;
  }

  return 0;
}

void 
FullGenLinSOE::setX(int loc, double value)
{
  if (loc < size && loc >=0)
	X[loc] = value;
}

void 
FullGenLinSOE::setX(const Vector &x)
{
  if (x.Size() == size)
    X = x;
}

const Vector &
FullGenLinSOE::getX()
{
  return X;
}

const Vector &
FullGenLinSOE::getB()
{
  return B;
}

const Matrix *
FullGenLinSOE::getA()
{
  assert(matA != nullptr);
  return matA;
}

double 
FullGenLinSOE::normRHS()
{
  double norm =0.0;
  for (int i=0; i<size; i++) {
      double Yi = B[i];
      norm += Yi*Yi;
  }
  return sqrt(norm);
}    


int
FullGenLinSOE::setFullGenSolver(FullGenLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);

    if (size != 0) {
        int solverOK = newSolver.setSize();
        if (solverOK < 0) {
                // opserr << "WARNING:FullGenLinSOE::setSolver :";
                // opserr << "the new solver could not setSize() - staying with old\n";
                return -1;
            }
    }
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
FullGenLinSOE::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int 
FullGenLinSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}


