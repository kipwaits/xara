// Attalla2D.cpp: implementation of the YieldSurfaceBC class.
//
//////////////////////////////////////////////////////////////////////

#include "Attalla2D.h"
#include <math.h>

#define ATTALLA_CLASS_TAG -1

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Attalla2D::Attalla2D(int tag, double xmax, double ymax,
                     YS_Evolution &model,
                     double a01, double a02, double a03,
                     double a04, double a05, double a06)
:YieldSurface_BC2D(tag, ATTALLA_CLASS_TAG, xmax, ymax, model),
 a1(a01), a2(a02), a3(a03), a4(a04), a5(a05), a6(a06), driftAlgo(0)

{

}

Attalla2D::~Attalla2D()
{

}

//////////////////////////////////////////////////////////////////////
// YS specific methods
//////////////////////////////////////////////////////////////////////
void Attalla2D::setExtent()
{
	// Extent along the axis
	xPos =  1;
	xNeg = -1;
	yPos =  0.98;
	yNeg = -0.98;
}


void Attalla2D::getGradient(double &gx, double &gy, double x, double y)
{
// check if the point is on the surface
    double drift =  getDrift(x, y);
    double loc   =  forceLocation(drift);
    double capx = capXdim;
    double capy = capYdim;

    if(loc != 0)
    {
     	opserr << "ERROR - Attalla2D::getGradient(double &gx, double &gy, double x, double y)\n";
        opserr << "Force point not on yield surface, drift = " << drift << " loc = " << loc <<"\n";
        opserr << "\a";
    }
    else
    {
		double a = 10.277;//3.043;//4.29293; //by compatibility, check err in gradient
		double yt = 0.95;

		if(y > yt)
		{
		 	gx = 2*a*x/capx;
			gy = 1;
		}
		else if(y < -yt)
		{
		 	gx = 2*a*x/capx;
			gy = -1; //-1*y/capy; <-- why did i write this? // -1 ?
		}
		else
		{
    		gx = 6*a2*pow(x,5)/capx + 4*a4*pow(x,3)/capx + 2*a6*x/capx;
        	gy = 6*a1*pow(y,5)/capy + 4*a3*pow(y,3)/capy + 2*a5*y/capy;
		}
    }

}

double Attalla2D::getSurfaceDrift(double x, double y)
{
double phi;
double a = 10.277; //3.043;//4.29293;
double yt = 0.95, xt = 0.054029;

	if(y > yt && fabs(x) < fabs(y)*xt/yt)
	{
		phi = a*x*x + y + 0.02;
	}
	else if(y < -yt && fabs(x) < fabs(y)*xt/yt)
	{
		phi = a*x*x - y + 0.02;
	}
	else
	{
		phi = a1*pow(y,6) + a2*pow(x,6) + a3*pow(y,4) + a4*pow(x,4) + a5*y*y + a6*x*x;
	}

double drift = phi - 1;
	return drift;
}

// Always call the base class method first
void Attalla2D::customizeInterpolate(double &xi, double &yi, double &xj, double &yj)
{
	this->YieldSurface_BC2D::customizeInterpolate(xi, yi, xj, yj);

double yt = 0.95, xt = 0.054029;

	if(fabs(yj) >= yt && fabs(xj) < fabs(yj)*xt/yt) //set to radial return
	{
		xi = 0;
		yi = 0;
	}
}


YieldSurface_BC *Attalla2D::getCopy(void)
{
    Attalla2D *theCopy = new Attalla2D(this->getTag(), capX, capY, *hModel,
                                       a1, a2, a3, a4, a5, a6);
    //later  copy all the state variables
    return theCopy;
}

void Attalla2D::Print(OPS_Stream &s, int flag)
{
    s << "\nYield Surface No: " << this->getTag() << " type: Attalla2D\n";
}

