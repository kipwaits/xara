//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Linearized Euler frame formulation with C1 displacement interpolation
//
// Adapted from DispBeamColumn3d?
//
#include <EulerFrame3d.h>
#include <Node.h>
#include <FrameSection.h>
#include <ID.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <Vector.h>
#include <Domain.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <FrameTransform.h>
#include <BeamIntegration.h>
#include <Parameter.h>

#include <array>
#include <math.h>
#include <string.h>

#include <BasicFrameTransf.h>
#include <runtime/commands/modeling/transform/FrameTransformBuilder.hpp>

#define ELE_TAG_EulerFrame3d 0 // TODO

using OpenSees::MatrixND;
using OpenSees::VectorND;


EulerFrame3d::EulerFrame3d(int tag, std::array<int,2>& nodes,
                           int numSec, FrameSection **s,
                           BeamIntegration &bi,
                           FrameTransformBuilder& tb,
                           double r, 
                           int cm)
 : FiniteElement<2, 3, 6>(tag, ELE_TAG_EulerFrame3d, nodes),
   basic_system(new BasicFrameTransf3d<6>(tb.template create<2,6>())),
   numSections(numSec),
   stencil(nullptr),
   density(r), 
   mass_flag(cm), 
   use_density(true),
   mass_initialized(false)
{
  q.zero();
  for (int i = 0; i < numSec; i++) {    
    // Get copies of the material for each integration point
    points.push_back({
        .point=0,
        .weight=0,
        .material=s[i]->getFrameCopy(scheme)
    });
  }
  
  stencil = bi.getCopy();
}


EulerFrame3d::~EulerFrame3d()
{    
  for (GaussPoint& point : points)
    if (point.material != nullptr)
      delete point.material;
  
  if (stencil != nullptr)
    delete stencil;

  if (basic_system != nullptr)
    delete basic_system;
}

int
EulerFrame3d::commitState()
{
  int status = 0;

  // Call element commitState to do any base class stuff
  if ((status = this->Element::commitState()) != 0)
    return status;

  for (GaussPoint& point : points) {
    if (point.material->commitState() != 0)
      return -1;
  }

  // Commit the transformation
  if ((status = basic_system->commitState()) != 0)
    return status;

  return status;
}

int
EulerFrame3d::revertToLastCommit()
{
  for (GaussPoint& point : points) {
    FrameSection& section = *point.material;

    if (section.revertToLastCommit() != 0)
      return -1;
  }

  if (basic_system->revertToLastCommit() != 0)
    return -2;

  return 0;
}

int
EulerFrame3d::revertToStart()
{

  for (GaussPoint& point : points) {
    if (point.material->revertToStart() != 0)
      return -1;
  }

  if (basic_system->revertToStart() != 0)
    return -2;

  return 0;
}


int
EulerFrame3d::setNodes()
{
  // call the parent class method
  int status = 0;

  if (basic_system->initialize(theNodes[0], theNodes[1]) != 0) {
      opserr << "BasicFrame3d::setDomain  tag: " 
            << this->getTag()
            << " -- Error initializing coordinate transformation\n";
      return -1;
  }

  double L = basic_system->getInitialLength();

  this->BasicFrame3d::setLength(L);

  if (L == 0.0)
    return -1;

  int numSections = points.size();
//double *xi = new double[numSections];
//double *wt = new double[numSections];
  stencil->getSectionLocations(numSections, L, xi);
  stencil->getSectionWeights(numSections, L, wt);
  for (int i=0; i<numSections; i++) {
    points[i].point  = xi[i];
    points[i].weight = wt[i];
  }
//delete[] xi;
//delete[] wt;

// NOTE(cmp) uncomment to match upstream behavior
//status += this->update();
//status += this->setState(State::Pres);
  return status;
}

int
EulerFrame3d::getIntegral(Field field, State state, double& total)
{

  if (this->setState(State::Init) != 0)
    return -1;

  total = 0.0;
  switch (field) {

    // Integrate density to compute total mass
    case Field::Density: {
      for (GaussPoint& sample : points) {
        double value = 0.0;
        // use element density if supplied
        if (use_density)
          total += sample.weight*density;

        // try using section's internal density
        else if (sample.material->getIntegral(Field::Density, state, value) == 0)
          total += sample.weight*value;

        else
          return -1;
      }
      return 0;
    }

    case Field::PolarInertia: {
      for (GaussPoint& sample : points) {
        double A;
        if (sample.material->getIntegral(Field::UnitYY, state, A) != 0)
          continue;

        // Get \int \rho y^2
        double Iz;
        if (sample.material->getIntegral(Field::DensityYY, state, Iz) != 0) {
          // Section does not support integrating density; try
          // integrating product of inertia and multiplying by rho
          if (use_density && sample.material->getIntegral(Field::UnitYY, state, Iz) == 0)
            Iz *= density/A;
          else
            continue;
        }
        // Get \int \rho z^2
        double Iy;
        if (sample.material->getIntegral(Field::DensityZZ, state, Iy) != 0) {
          if (use_density && sample.material->getIntegral(Field::UnitZZ, state, Iy) == 0)
            Iy *= density/A;
          else
            continue;
        }
        total += sample.weight*(Iy + Iz);
      }
      return 0;
    }

    default:
      return -1;
  }
}

int
EulerFrame3d::update()
{
  // Note that setNodes must be called prior to calling this
  // so that basic_system, xi and wt are initialized
  int err = 0;

  // Update the transformation
  basic_system->update();

  double jsx = 1.0/basic_system->getInitialLength();
  
  // Get basic deformations
  const Vector &v = basic_system->getBasicTrialDisp();
  
  //
  // Gauss points
  //
  for (int i = 0; i < points.size(); i++) {
      double xi6 = 6.0*xi[i];

      const VectorND<nsr> e {
           jsx*v[0],                                 // N
           jsx*(xi6-4.0)*v[1] + jsx*(xi6-2.0)*v[2],  // Mz
           jsx*(xi6-4.0)*v[3] + jsx*(xi6-2.0)*v[4],  // My
           jsx*v[5],                                 // T
      };
      
      // Set the section deformations
      err += points[i].material->setTrialState<nsr, scheme>(e);
  }

  return err;
}

VectorND<6>&
EulerFrame3d::getBasicForce()
{
  return q;
}


const Vector &
EulerFrame3d::getResistingForce()
{
  this->getBasicTangent(State::Pres, 0);

  double q0 = q[0];
  double q1 = q[1];
  double q2 = q[2];
  double q3 = q[3];
  double q4 = q[4];
  double q5 = q[5];

  double oneOverL = 1.0 / basic_system->getInitialLength();

  thread_local VectorND<12> pl;
  pl[0]  = -q0;                    // Ni
  pl[1]  =  oneOverL * (q1 + q2);  // Viy
  pl[2]  = -oneOverL * (q3 + q4);  // Viz
  pl[3]  = -q5;                    // Ti
  pl[4]  =  q3;
  pl[5]  =  q1;
  pl[6]  =  q0;                    // Nj
  pl[7]  = -pl[1];                 // Vjy
  pl[8]  = -pl[2];                 // Vjz
  pl[9]  = q5;                     // Tj
  pl[10] = q4;
  pl[11] = q2;

  thread_local VectorND<12> pf{0.0};
  pf[0] = p0[0];
  pf[1] = p0[1];
  pf[7] = p0[2];
  pf[2] = p0[3];
  pf[8] = p0[4];


  thread_local VectorND<12> pg;
  thread_local Vector wrapper(pg);
//    const Vector p0Vec(p0);
//    P = basic_system->getGlobalResistingForce(q, p0Vec);
  pg  = basic_system->t.pushResponse(pl);
  pg += basic_system->t.pushConstant(pf);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (total_mass != 0.0)
    wrapper.addVector(1.0, p_iner, -1.0);

  return wrapper;
}


MatrixND<6,6>&
EulerFrame3d::getBasicTangent(State state, int rate)
{

  // Zero for integral
  kb.zero();
  q.zero();
  
  double L = basic_system->getInitialLength();
  double jsx = 1.0/L;

  //
  // Gauss loop
  //
  for (int i = 0; i < numSections; i++) {
    GaussPoint& point = points[i];

    MatrixND<4,6> ka;
    ka.zero();

    double xi6 = 6.0*point.point;

    /*  */ const MatrixND<nsr,nsr> ks = 
      point.material->getTangent<nsr,scheme>(state);
               
    // Perform numerical integration
    // kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    double wti = points[i].weight*jsx;

    for (int k = 0; k < 4; k++) {
      // N
      ka(k,0) += ks(k,0)*wti;
      // Mz
      double tmp = ks(k,1)*wti;
      ka(k,1) += (xi6-4.0)*tmp;
      ka(k,2) += (xi6-2.0)*tmp;
      // My
      tmp = ks(k,2)*wti;
      ka(k,3) += (xi6-4.0)*tmp;
      ka(k,4) += (xi6-2.0)*tmp;
      // T
      ka(k,5) += ks(k,3)*wti;
    }

    // Beam terms
    for (int k = 0; k < 6; k++) {
      // N
      kb(0,k) += ka(0,k);
      // Mz
      double tmp = ka(1,k);
      kb(1,k) += (xi6-4.0)*tmp;
      kb(2,k) += (xi6-2.0)*tmp;
      // My
      tmp = ka(2,k);
      kb(3,k) += (xi6-4.0)*tmp;
      kb(4,k) += (xi6-2.0)*tmp;
      // T
      kb(5,k) += ka(3,k);
    }

    if (state != State::Init) {
      const VectorND<nsr> s = point.material->getResultant<nsr,scheme>();
      // q.addMatrixTransposeVector(1.0, *B, s, wts(i));
      q[0] += s[0]*wt[i];
      q[1] += (xi6-4.0)*s[1]*wt[i];
      q[2] += (xi6-2.0)*s[1]*wt[i];
      q[3] += (xi6-4.0)*s[2]*wt[i];
      q[4] += (xi6-2.0)*s[2]*wt[i];
      q[5] += s[3]*wt[i];
    }  
  }

  q += q0;
  
  return kb;
}


const Matrix &
EulerFrame3d::getTangentStiff()
{
  MatrixND<6,6> kb = this->getBasicTangent(State::Pres, 0);
  VectorND<6>   q  = this->getBasicForce();

  return basic_system->getGlobalStiffMatrix(Matrix(kb), Vector(q));
}


const Matrix &
EulerFrame3d::getInitialStiff()
{
  return basic_system->getInitialGlobalStiffMatrix(this->getBasicTangent(State::Init, 0));
}

int
EulerFrame3d::sendSelf(int tag, Channel &c)
{
  return -1;
}

int
EulerFrame3d::recvSelf(int tag, Channel &c, FEM_ObjectBroker &b)
{
  return -1;
}

void
EulerFrame3d::Print(OPS_Stream &s, int flag)
{
  int numSections = points.size();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << OPS_PRINT_JSON_ELEM_INDENT << "{";
      s << "\"name\": " << this->getTag() << ", ";
      s << "\"type\": \"EulerFrame3d\"" << ", ";
      s << "\"nodes\": [" << connectedExternalNodes(0) << ", " 
                          << connectedExternalNodes(1) << "]";
      s << ", ";

      // Mass
      double mass;
      if (getIntegral(Field::Density, State::Init, mass) == 0)
        s << "\"mass\": " << mass;
      else
        s << "\"massperlength\": " << density;
      s << ", ";

      // Transform
      s << "\"transform\": " << basic_system->getTag();
      s << ", ";

      //
      s << "\"sections\": [";
      for (int i = 0; i < numSections - 1; i++)
              s << points[i].material->getTag() << ", ";
      s << points[numSections - 1].material->getTag() << "], ";
      s << "\"integration\": ";
      stencil->Print(s, flag);
      //
      s << "}";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nEulerFrame3d, element id:  " << this->getTag() << "\n";
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tCoordTransf: " << basic_system->getTag() << "\n";
        s << "\tmass density:  " << density << ", mass_flag: " << mass_flag << "\n";

        double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
        double L = basic_system->getInitialLength();
        double jsx = 1.0 / L;

        N = q[0];
        Mz1 = q[1];
        Mz2 = q[2];
        Vy = (Mz1 + Mz2)*jsx;
        My1 = q[3];
        My2 = q[4];
        Vz = -(My1 + My2)*jsx;
        T = q[5];

        s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
                << -N + p0[0] << ' ' << Mz1 << ' ' << Vy + p0[1] << ' ' << My1 << ' ' << Vz + p0[3] << ' ' << -T << "\n";
        s << "\tEnd 2 Forces (P Mz Vy My Vz T): "
                << N << ' ' << Mz2 << ' ' << -Vy + p0[2] << ' ' << My2 << ' ' << -Vz + p0[4] << ' ' << T << "\n";
        s << "Number of sections: " << numSections << "\n";
        stencil->Print(s, flag);

        for (int i = 0; i < numSections; i++)
          points[i].material->Print(s, flag);

        //  if (density != 0)
        //    opserr << "Mass: \n" << this->getMass();
  }
}


const Matrix &
EulerFrame3d::getMass()
{
    if (!mass_initialized) {
      if (this->getIntegral(Field::Density, State::Init, total_mass) != 0)
        ;
      if (this->getIntegral(Field::PolarInertia, State::Init, twist_mass) != 0)
        ;
      mass_initialized = true;
    }

    if (total_mass == 0.0) {

        thread_local MatrixND<12,12> M{0.0};
        thread_local Matrix Wrapper{M};
        return Wrapper;

    } else if (mass_flag == 0)  {

        thread_local MatrixND<12,12> M{0.0};
        thread_local Matrix Wrapper{M};
        // lumped mass matrix
        double m = 0.5*total_mass;
        M(0,0) = m;
        M(1,1) = m;
        M(2,2) = m;
        M(6,6) = m;
        M(7,7) = m;
        M(8,8) = m;
        return Wrapper;

    } else {
        // consistent (cubic, prismatic) mass matrix

        double L  = basic_system->getInitialLength();
        double m  = total_mass/420.0;
        double mx = twist_mass;
        thread_local MatrixND<12,12> ml{0};

        ml(0,0) = ml(6,6) = m*140.0;
        ml(0,6) = ml(6,0) = m*70.0;

        ml(3,3) = ml(9,9) = mx/3.0; // Twisting
        ml(3,9) = ml(9,3) = mx/6.0;

        ml( 2, 2) = ml( 8, 8) =  m*156.0;
        ml( 2, 8) = ml( 8, 2) =  m*54.0;
        ml( 4, 4) = ml(10,10) =  m*4.0*L*L;
        ml( 4,10) = ml(10, 4) = -m*3.0*L*L;
        ml( 2, 4) = ml( 4, 2) = -m*22.0*L;
        ml( 8,10) = ml(10, 8) = -ml(2,4);
        ml( 2,10) = ml(10, 2) =  m*13.0*L;
        ml( 4, 8) = ml( 8, 4) = -ml(2,10);

        ml( 1, 1) = ml( 7, 7) =  m*156.0;
        ml( 1, 7) = ml( 7, 1) =  m*54.0;
        ml( 5, 5) = ml(11,11) =  m*4.0*L*L;
        ml( 5,11) = ml(11, 5) = -m*3.0*L*L;
        ml( 1, 5) = ml( 5, 1) =  m*22.0*L;
        ml( 7,11) = ml(11, 7) = -ml(1,5);
        ml( 1,11) = ml(11, 1) = -m*13.0*L;
        ml( 5, 7) = ml( 7, 5) = -ml(1,11);

        // transform local mass matrix to global system
        return basic_system->getGlobalMatrixFromLocal(ml);
    }
}




Response*
EulerFrame3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","EulerFrame3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    //
    // Compare argv[0] for known response types 
    //

    // Global force
    if (strcmp(argv[0],"forces") == 0 || 
        strcmp(argv[0],"force") == 0  ||
        strcmp(argv[0],"globalForce") == 0 ||
        strcmp(argv[0],"globalForces") == 0) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");


      theResponse = new ElementResponse(this, 1, Vector(12));

    // Local force
    }  else if (strcmp(argv[0],"localForce") == 0 || 
                strcmp(argv[0],"localForces") == 0) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Vy_2");
      output.tag("ResponseType","Vz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");

      theResponse = new ElementResponse(this, 2, Vector(12));

    // chord rotation -
    }  else if (strcmp(argv[0],"chordRotation") == 0 || 
                strcmp(argv[0],"chordDeformation") == 0 || 
                strcmp(argv[0],"basicDeformation") == 0) {

      output.tag("ResponseType","eps");
      output.tag("ResponseType","thetaZ_1");
      output.tag("ResponseType","thetaZ_2");
      output.tag("ResponseType","thetaY_1");
      output.tag("ResponseType","thetaY_2");
      output.tag("ResponseType","thetaX");

      theResponse = new ElementResponse(this, 3, Vector(6));

    // plastic rotation -
    } else if (strcmp(argv[0],"plasticRotation") == 0 || 
               strcmp(argv[0],"plasticDeformation") == 0) {

      output.tag("ResponseType","epsP");
      output.tag("ResponseType","thetaZP_1");
      output.tag("ResponseType","thetaZP_2");
      output.tag("ResponseType","thetaYP_1");
      output.tag("ResponseType","thetaYP_2");
      output.tag("ResponseType","thetaXP");

      theResponse = new ElementResponse(this, 4, Vector(6));
  

    } else if (strcmp(argv[0],"RayleighForces") == 0 ||
               strcmp(argv[0],"rayleighForces") == 0) {
  
      theResponse =  new ElementResponse(this, 12, Vector(12));
  
    }   
    else if (strcmp(argv[0],"integrationPoints") == 0)
      theResponse = new ElementResponse(this, 10, Vector(points.size()));

    else if (strcmp(argv[0],"integrationWeights") == 0)
      theResponse = new ElementResponse(this, 11, Vector(points.size()));

    else if (strcmp(argv[0],"sectionTags") == 0)
      theResponse = new ElementResponse(this, 110, ID(points.size()));
    
    // section response -
    else if (strstr(argv[0],"sectionX") != 0) {
      if (argc > 2 && (this->setNodes() == 0)) {

        float sectionLoc = atof(argv[1]);
        double L = basic_system->getInitialLength();

        sectionLoc /= L;
        
        float minDistance = fabs(xi[0]-sectionLoc);
        int sectionNum = 0;
        for (int i = 1; i < numSections; i++) {
          if (fabs(xi[i]-sectionLoc) < minDistance) {
            minDistance = fabs(xi[i]-sectionLoc);
            sectionNum = i;
          }
        }
        
        output.tag("GaussPointOutput");
        output.attr("number",sectionNum+1);
        output.attr("eta", points[sectionNum].point * L);
        
        theResponse = points[sectionNum].material->setResponse(&argv[2], argc-2, output);
      }
    }
    
    else if (strcmp(argv[0],"section") ==0) { 
      // Make sure setNodes() succeeds to ensure
      // xi is initialized and L can be determined
      if (argc > 1 && (this->setNodes() == 0)) {
        
        int sectionNum = atoi(argv[1]);
        double L = basic_system->getInitialLength();

        if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

          output.tag("GaussPointOutput");
          output.attr("number",sectionNum);
          output.attr("eta",xi[sectionNum-1]*L);
          
          theResponse =  points[sectionNum-1].material->setResponse(&argv[2], argc-2, output);
          
          output.endTag();

        } else if (sectionNum == 0) { 
          // argv[1] was not an int, we want all sections, 
        
          CompositeResponse *theCResponse = new CompositeResponse();
          int numResponse = 0;
          
          for (int i=0; i<numSections; i++) {
            
            output.tag("GaussPointOutput");
            output.attr("number",i+1);
            output.attr("eta", xi[i]*L);
            
            Response *theSectionResponse = points[i].material->setResponse(&argv[1], argc-1, output);
            
            output.endTag();          
            
            if (theSectionResponse != 0) {
              numResponse = theCResponse->addResponse(theSectionResponse);
            }
          }
          
          if (numResponse == 0) // no valid responses found
            delete theCResponse;
          else
            theResponse = theCResponse;
        }
      }
    }
        // by SAJalali
    else if (strcmp(argv[0], "energy") == 0) {
      return new ElementResponse(this, 13, 0.0);
    }

    if (theResponse == nullptr)
      theResponse = basic_system->setResponse(argv, argc, output);

  output.endTag();
  return theResponse;
}

int 
EulerFrame3d::getResponse(int responseID, Information &info)
{

  if (responseID == 1)
    return info.setVector(this->getResistingForce());

  else if (responseID == 12)
    return info.setVector(this->getRayleighDampingForces());
    
  else if (responseID == 2) {
    double V, M1, M2, T;
    double L = basic_system->getInitialLength();
    static Vector P(12);

    // Axial
    double N = q[0];
    P(6) =  N;
    P(0) = -N+p0[0];
    
    // Torsion
    T = q[5];
    P(9) =  T;
    P(3) = -T;
    
    // Moments about z and shears along y
    M1 = q[1];
    M2 = q[2];
    P(5)  = M1;
    P(11) = M2;
    V = (M1 + M2)/L;
    P(1) =  V+p0[1];
    P(7) = -V+p0[2];
    
    // Moments about y and shears along z
    M1 = q[3];
    M2 = q[4];
    P(4)  = M1;
    P(10) = M2;
    V = (M1 + M2)/L;
    P(2) = -V+p0[3];
    P(8) =  V+p0[4];

    return info.setVector(P);
  }

  // Chord rotation
  else if (responseID == 3)
    return info.setVector(basic_system->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(6);
    static Vector ve(6);
    auto kb = this->getBasicTangent(State::Init, 0);
    kb.solve(q, ve);
    vp = basic_system->getBasicTrialDisp();
    vp -= ve;
    return info.setVector(vp);
  }

  else if (responseID == 10) {
    // ensure we have L, xi[] and wt[]
    if (this->setState(State::Init) != 0)
      return -1;

    double L = basic_system->getInitialLength();
    Vector locs(points.size());
    for (int i = 0; i < numSections; i++)
      locs[i] = wt[i]*L;
    return info.setVector(locs);
  }

  else if (responseID == 11) {
    if (this->setState(State::Init) == 0) {
      double L = basic_system->getInitialLength();
      Vector weights(numSections);
      for (int i = 0; i < numSections; i++)
        weights(i) = wt[i]*L;
      return info.setVector(weights);
    }
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = points[i].material->getTag();
    return info.setID(tags);
  }

  //by SAJalali
  else if (responseID == 13) {
    if (this->setState(State::Init) == 0) {
      double L = basic_system->getInitialLength();
      double energy = 0;
      for (int i = 0; i < numSections; i++)
          energy += points[i].material->getEnergy()*xi[i] * L;

      return info.setDouble(energy);
    }
  }

  else
    return -1;

  return -1;
}


int
EulerFrame3d::setParameter(const char **argv, int argc, Parameter &param)
{

  int status = this->FiniteElement<2,3,6>::setParameter(argv, argc, param);
  if (status != -1)
    return status;

  if (argc < 1)
    return -1;

  if (strstr(argv[0],"sectionX") != 0) {
      if (argc < 3 || (this->setState(State::Init) != 0))
        return -1;
      
      float sectionLoc = atof(argv[1]);
      double L = basic_system->getInitialLength();
      
      sectionLoc /= L;

      float minDistance = fabs(xi[0]-sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i] - sectionLoc) < minDistance) {
          minDistance = fabs(xi[i]-sectionLoc);
          sectionNum = i;
        }
      }
      return points[sectionNum].material->setParameter(&argv[2], argc-2, param);
  }

  // If the parameter belongs to a section or lower
  if (strstr(argv[0],"section") != 0) {
    
    if (argc < 3)
      return -1;
    
    // Get section number
    int sectionNum = atoi(argv[1]) - 1;

    if (sectionNum >= 0 && sectionNum < numSections)
      return points[sectionNum].material->setParameter(&argv[2], argc-2, param);
    else
      return -1;
  }
  
  if (strstr(argv[0], "integration") != 0) {
    
    if (argc < 2)
      return -1;

    return stencil->setParameter(&argv[1], argc-1, param);
  }

  // Default, send to every object
  int ok = 0;
  int result = -1;

  for (int i = 0; i < numSections; i++) {
    ok = points[i].material->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  
  ok = stencil->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  return result;
}



const Matrix &
EulerFrame3d::getInitialStiffSensitivity(int gradNumber)
{
  static Matrix K(12,12);
  K.Zero();
  return K;
}

VectorND<6>
EulerFrame3d::getBasicForceGrad(int gradNumber)
{
  //
  // Return the *conditional* derivative:
  //   dqdh|_v
  //
  //
  // See Haukaas and Scott (2006)
  //
  VectorND<6> dqdh;
  dqdh.zero();
  
  // Perform numerical integration
  for (int i = 0; i < numSections; i++) {

    int order = points[i].material->getOrder();
    const ID &code = points[i].material->getType();
    
    double xi6 = 6.0*xi[i];
    double wti = wt[i];
    
    // Get section stress resultant gradient
    const Vector &dsdh = points[i].material->getStressResultantSensitivity(gradNumber, true);
    
    for (int j = 0; j < order; j++) {
      double sensi = dsdh[j]*wti;
      switch(code(j)) {
      case SECTION_RESPONSE_P:
        dqdh(0) += sensi; 
        break;
      case SECTION_RESPONSE_MZ:
        dqdh(1) += (xi6-4.0)*sensi; 
        dqdh(2) += (xi6-2.0)*sensi; 
        break;
      case SECTION_RESPONSE_MY:
        dqdh(3) += (xi6-4.0)*sensi; 
        dqdh(4) += (xi6-2.0)*sensi; 
        break;
      case SECTION_RESPONSE_T:
        dqdh(5) += sensi; 
        break;
      default:
        break;
      }
    }
  }

  if (basic_system->isShapeSensitivity()) { 
 
    MatrixND<6,6> &Kb = getBasicTangent(State::Pres, 0);
    double L   = basic_system->getInitialLength();
    
    const Vector &A_u = basic_system->getBasicTrialDisp();
    double dLdh = basic_system->getLengthGrad();
    double d1overLdh = -dLdh/(L*L);

    // a^T k_s dadh v
    dqdh.addMatrixVector(1.0, Kb, A_u, d1overLdh);
  }

  return dqdh;
}

const Vector &
EulerFrame3d::getResistingForceSensitivity(int gradNumber)
{
  //
  // dAdh^T q + A^T (dqdh + k dAdh u)
  //
  static Vector P(12);
  P.Zero();

  VectorND<6> dqdh = this->getBasicForceGrad(gradNumber);

  double jsx = 1.0/basic_system->getInitialLength();

  // TODO: No distributed loads
  static Vector dp0dh(6);   

  if (basic_system->isShapeSensitivity()) {
    // k dAdh u
    const Vector &dAdh_u = basic_system->getBasicDisplFixedGrad();
    dqdh.addMatrixVector(1.0, kb, dAdh_u, jsx);

    // dAdh^T q
    P = basic_system->getGlobalResistingForceShapeSensitivity(q, dp0dh, gradNumber);
  }


  // A^T (dqdh + k dAdh u)
  P += basic_system->getGlobalResistingForce(dqdh, dp0dh);
  
  return P;
}


int
EulerFrame3d::commitSensitivity(int gradNumber, int numGrads)
{
  // Get basic deformation and sensitivities
  const Vector &v = basic_system->getBasicTrialDisp();
  
  static Vector dvdh(6);
  dvdh = basic_system->getBasicDisplTotalGrad(gradNumber);
  
  double L = basic_system->getInitialLength();
  double jsx = 1.0/L;

  double d1oLdh = basic_system->getd1overLdh();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = points[i].material->getOrder();
    const ID &code = points[i].material->getType();
    
    Vector e(order);
    
    double xi6 = 6.0*xi[i];
    
    for (int j = 0; j < nsr; j++) {
      switch(scheme[j]) {
      case SECTION_RESPONSE_P:
        e[j] = jsx*dvdh(0) + d1oLdh*v(0); 
        break;
      case SECTION_RESPONSE_MZ:
        e[j] = jsx*((xi6-4.0)*dvdh(1) + (xi6-2.0)*dvdh(2))
             + d1oLdh*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)); 
        break;
      case SECTION_RESPONSE_MY:
        e[j] = jsx*((xi6-4.0)*dvdh(3) + (xi6-2.0)*dvdh(4))
             + d1oLdh*((xi6-4.0)*v(3) + (xi6-2.0)*v(4)); 
        break;
      case SECTION_RESPONSE_T:
        e[j] = jsx*dvdh(5)
             + d1oLdh*v(5); 
        break;
      default:
        e[j] = 0.0; 
        break;
      }
    }

    // Set the section deformations
    points[i].material->commitSensitivity(e,gradNumber,numGrads);
  }
  
  return 0;
}

