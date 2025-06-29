//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
//
// Written in Matlab: Thanh Do
// Created: 07/16
//
#include "FariaPlasticDamage3d.h"
#include <Channel.h>
#include <MatrixND.h>
#include <cmath>
#include <Projector.hh>

#define MND
using namespace OpenSees;


FariaPlasticDamage3d::FariaPlasticDamage3d(int tag,
                                           double E, 
                                           double nu, 
                                           double Ft,
                                           double Fc, 
                                           double beta, 
                                           double Ap,
                                           double An, 
                                           double Bn,
                                           double density)
: NDMaterial(tag,ND_TAG_PlasticDamageConcrete3d),
  E(E), nu(nu),
  ft(Ft), Fc(Fc), beta(beta), Ap(Ap), An(An), Bn(Bn),
  density(density),
  retTangent(C),
  retStress(sig),
  retStrain(eps),
  retInitialTangent(Ce)
{
  this->revertToStart();
  this->commitState();
}

FariaPlasticDamage3d::FariaPlasticDamage3d()
  : NDMaterial(0, ND_TAG_PlasticDamageConcrete3d)
{

}

FariaPlasticDamage3d::~FariaPlasticDamage3d()
{

}


// Stress invariants: octahedral normal and shear
static inline void
StrsInvar(const VectorND<6> &sig, double &sigoct, double &tauoct)
{
  // normal stress
  sigoct = (sig(0) + sig(1) + sig(2))/3.;
  // shear stress
  double J2 = (std::pow((sig(0) - sig(1)),2) 
            +  std::pow((sig(0) - sig(2)),2)
            +  std::pow((sig(1) - sig(2)),2))/6.
            +  std::pow(sig(3),2) + std::pow(sig(4),2) + std::pow(sig(5),2);
  tauoct = std::sqrt(2./3.*J2);  
}

//
// Stress decomposition function: algebraic approach
//
#if 0
#include <concrete/StrsDec.cpp>
#else
template <typename T> static inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
static inline void
StrsDecA(const Vector &sig, VectorND<6> &sigpos, //VectorND<6> &signeg, 
                            MatrixND<6,6> *Qpos) //,   MatrixND<6,6> *Qneg)
{
#if 1
  for (int i=0; i<6; i++) {
    // positive and negative stress tensors
    sigpos(i) = (sig(i) + std::fabs(sig(i))) / 2.;
  }

  // projection tensors
  if (Qpos == nullptr)
    return;

  for (int i=0; i<6; i++) {
    (*Qpos)(i,i) = (1. + sgn(sig(i))) / 2.;
  }
  return;
#else
  for (int i=0; i<6; i++) {
    if (sig(i) > 1e-8) {
      sigpos(i) = sig(i);
      signeg(i) = 0.;
    } else if (sig(i)  < -1.e-8) {
      sigpos(i) = 0.;
      signeg(i) = sig(i);
    } else {
      sigpos(i) = sig(i)/2.;
      signeg(i) = sig(i)/2.;
    }
  }
  if (Qpos == nullptr || Qneg == nullptr)
    return;
  for (int i=0; i<6; i++) {
    if (sig(i) > 1e-8) {
      (*Qpos)(i,i) = 1;
      (*Qneg)(i,i) = 0;
    } else if (sig(i)  < -1.e-8) {
      (*Qpos)(i,i) = 0;
      (*Qneg)(i,i) = 1;
    } else {
      (*Qpos)(i,i) = 0.5;
      (*Qneg)(i,i) = 0.5;
    }
  }
#endif
}
#endif


// Compute dn = 1 - (rn0/rn)*(1-An) - An*exp(Bn*(1 - rn/rn0))
double
negative_damage(double rn0, double rn, double An, double Bn) {
    // 1) Compute alpha = rn0/rn
    double alpha = rn0 / rn;

    // 2) Compute x = Bn * (1 - rn/rn0) without cancellation
    double inv_alpha = 1.0 / alpha;
    double x = Bn * (1.0 - inv_alpha);

    double exm1 = std::expm1(x);

    // 4) Combine terms. Use fma to reduce rounding:
    //    term1 = (1 - alpha)*(1 - An)
    double term1 = std::fma(-(1.0 - An), alpha, (1.0 - An));  
    //    which is = (1 - An) - alpha*(1 - An)
    //
    double term2 = An * exm1;

    // 5) dn = term1 - term2
    double dn = std::fma(-An, exm1, term1);
    return dn;
}

// Compute dp = 1 - (rp0/rp)*exp(Ap*(1 - rp/rp0))
// with improved precision when rp ≈ rp0
double
positive_damage(double rp0, double rp, double Ap) {
    // 1) alpha = rp0 / rp
    double alpha = rp0 / rp;

    // 2) x = Ap * (1 - rp/rp0) = Ap * (1 - 1/alpha)
    double inv_alpha = 1.0 / alpha;
    double x = Ap * (1.0 - inv_alpha);

    // 3) exm1 = exp(x) - 1 computed accurately
    double exm1 = std::expm1(x);

    // 4)
    double term1 = 1.0 - alpha;

    // 5) dp = term1 - alpha * exm1
    //    use fma to compute -(alpha * exm1) + term1 in one rounding
    double dp = std::fma(-alpha, exm1, term1);

    return dp;
}

static double
negative_surface(const VectorND<6> &signeg, double rn, double k)
{
  double sigoct, tauoct;
  StrsInvar(signeg, sigoct, tauoct);                  //  find octahedral stresses
  return std::sqrt(std::sqrt(3.0)*(k*sigoct + tauoct)) - rn;    //  negative equivalent stress
}

int
FariaPlasticDamage3d::setTrialStrain(const Vector &strain)
{

  double F2c = 1.16*Fc; // f2c: biaxial compressive strength
  double k = std::sqrt(2.0)*(F2c - Fc)/(2.*F2c - Fc);
  // initial damage threshold
  double rp0 = ft/sqrt(E);
  double rn0 = sqrt((-k + sqrt(2.0))*Fc/std::sqrt(3.0));

  constexpr static double toln = 1e-10,
                          tolp = 1e-10;
  constexpr static double dn_max = 1.0 - 1e-10,
                          dp_max = 1.0 - 1e-7;

  // retrieve history variables
  eps_p = eps_pCommit;
  VectorND<6> sigeP = sigeCommit;
  rp = rpCommit;
  rn = rnCommit;
  dp = dpCommit;
  dn = dnCommit;

   // elastic compliance
  MatrixND<6,6> Se;
  Ce.invert(Se);

  // current strain
  eps = strain;

  // incremental strain
  VectorND<6> Depse_tr = eps - eps_p; // ?[cmp]
  VectorND<6> Deps = eps - epsCommit;


  //
  // PLASTIC part
  //

  // elastic trial stress
  VectorND<6> sige_tr  =  sigeCommit;
  sige_tr += Ce*Deps;

  // decomposition of trial stress tensor
  VectorND<6> sigpos{}, signeg = sige_tr;
  StrsDecA(sige_tr, sigpos, nullptr);
  signeg -= sigpos;

  MatrixND<6,6> Cbar{};
  if (negative_surface(signeg, rn, k) <= (toln*rn0)) {
    // elastic state, accept trial response
    sige = sige_tr;                                                    
    Cbar = Ce;                                              
  }
  else {
    // Correction
  
    //  norm of trial effective stress
    double nrm = sqrt(pow(sige_tr(0),2) + pow(sige_tr(1),2) + pow(sige_tr(2),2) 
               + 2*pow(sige_tr(3),2) + 2*pow(sige_tr(4),2) + 2*pow(sige_tr(5),2));

    // normalized trial effective stress
    VectorND<6> L_tr = sige_tr/nrm;

    double L_trDotDeps = L_tr.dot(Deps);

    // plastic strain increment
    VectorND<6> Deps_p;
    Deps_p.addVector(0.0, Depse_tr,  beta*E/nrm*L_trDotDeps);

    double lam  = 1.0 - beta*E/nrm * L_trDotDeps;

    sige.addVector(0.0, sige_tr, lam);                 //  corrected effective stress

    // check damage
    signeg = sige;
    StrsDecA(sige, sigpos, nullptr);                    //  decompose the effective stress  
    signeg -= sigpos;

    if ((negative_surface(signeg, rn, k) <= toln*rn0) || (L_trDotDeps <= 0.0)) {
      //  no damage or sige and eps in different direction
      sige = sige_tr;
      Cbar = Ce;
    }

    else {
      //  update plastic strain
      eps_p += Deps_p;

      VectorND<6> L_tr_temp {};
      // tangent in effective space, Cbar
      for (int i=0; i<3; i++) 
        L_tr_temp(i) = L_tr(i);
      for (int i=3; i<6; i++)
        L_tr_temp(i) = 2*L_tr(i);

      double Dlam_Dnrm = 2.0*beta*E/pow(nrm,3)*sige_tr.dot(Deps);

      VectorND<6> Dlam_Deps = L_tr;
      Dlam_Deps *= -beta*E/nrm;       // Dlam_Deps = -beta*E/nrm*L_tr;
      // Dlam_Deps = Dlam_Dnrm * Ce * Dnrm_Dsig + Ce*Dlam_Dsig + Dlam_Deps; 
      Dlam_Deps.addMatrixVector(1.0, Ce, L_tr_temp,         Dlam_Dnrm);
      Dlam_Deps.addMatrixVector(1.0, Ce,      Deps, -beta*E/(nrm*nrm));

      Cbar.addMatrix(Ce, lam);
      Cbar.addTensorProduct(sige_tr, Dlam_Deps, 1.0);
    }
  }


  //
  // DAMAGE part
  //

  // decompose into positive and negative effective stress tensor
  MatrixND<6,6> Qpos{}, Qneg = OpenSees::IImix;
  signeg = sige;
  StrsDecA(sige, sigpos, &Qpos);    // decompose the effective stress
  signeg -= sigpos;                 // signeg = sige - sigpos
  Qneg   -= Qpos;

  // positive damage
  double taup,
         Ddp_Drp = 0.;
  {
  // calculate equivalent stresses
    VectorND<6> tmp = Se*sigpos; // {}
    // Ce.solve(sigpos, tmp);
    taup = std::sqrt(tmp.dot(sigpos));      // positive equivalent stress


    if ((taup - rp) <= (tolp*rp0)) {               // no positive damage
      Ddp_Drp = 0;
    }
    else {                                         // positive damage evolves
      rp = taup;                                   // update rp = max(taup, rp)
      // dp = 1. - rp0/rp * std::exp(Ap*(1. - rp/rp0));
      dp = positive_damage(rp0, rp, Ap);          // update dp
      Ddp_Drp =  (Ap*rp + rp0)/(rp*rp) * std::exp(Ap*(1. - rp/rp0));        
      if (dp < 0) {
        dp = 0;
        Ddp_Drp = 0;
      }       
      // dp = dp*dp_max;                             // cap the damage variable 
      // Ddp_Drp = Ddp_Drp*dp_max;
      if (dp > dp_max) {
        dp = dp_max; 
        Ddp_Drp = 0;
      }
    }
  }

  // negative damage
  double Ddn_Drn = 0;
  double gn = negative_surface(signeg, rn, k); // negative equivalent stress
  {
    if (gn <= toln*rn0) {                    // no negative damage
      Ddn_Drn = 0;
    }
    else {                                  // negative damage evolves
      // rn = taun;                         // update rn
      rn += gn;
      // dn = 1. - rn0/rn*(1-An) - An*std::exp(Bn*(1. - rn/rn0));
      dn = negative_damage(rn0, rn, An, Bn); // update dn
      Ddn_Drn = rn0/(rn*rn)*(1.-An) + An*Bn/rn0*std::exp(Bn*(1. - rn/rn0));
      if (dn < 0) {
        dn = 0;
        Ddn_Drn = 0;
      }
      // cap the damage variable
      if (dn > dn_max) {
        dn = dn_max;
        Ddn_Drn = 0;
      }
    }
  }

  // stress update
  sig.addVector(0.0, sigpos, 1.0-dp);
  sig.addVector(1.0, signeg, 1.0-dn);

  //
  // TANGENT
  //
  {
    C.zero();

    // Negative part
    {
      MatrixND<6,6> Dsigneg_Deps = Qneg*Cbar;

      VectorND<6> Dtaun_Dsigneg{};
      if (gn <= toln) {
        Dtaun_Dsigneg.zero();
      }
      else {
        double sigoct, tauoct;
        StrsInvar(signeg, sigoct, tauoct);             // find octahedral stresses
        double taun = sqrt((sqrt(3.)*(k*sigoct + tauoct)));   // negative equivalent stress

        // norm of deviatoric stress
        VectorND<6> s = IIdevMix*signeg;            // deviatoric stress

        double nrms = sqrt( pow(s(0),2) + pow(s(1),2) + pow(s(2),2) +  
                          2*pow(s(3),2) + 2*pow(s(4),2) + 2*pow(s(5),2));
        VectorND<6> n = s;
        if (abs(nrms) <= 1e-8) //toln) 
          n.zero();
        else {
          n/=nrms;
          double Dtaun_Dtauoct = pow(3,0.25) / 2./sqrt(k*sigoct + tauoct);
          Dtaun_Dsigneg.addVector(1.0,    n, Dtaun_Dtauoct/std::sqrt(3.0)); // += Dtaun_Dtauoct * Dtauoct_Dsigneg
        }

        double Dtaun_Dsigoct = pow(3,0.25) * k/2./sqrt(k*sigoct + tauoct);
        Dtaun_Dsigneg.addVector(0.0, ivol, 1./3.0*Dtaun_Dsigoct); // = Dtaun_Dsigoct * Dsigoct_Dsigneg
      }

      VectorND<6> Ddn_Deps = Dsigneg_Deps ^ Dtaun_Dsigneg;
      Ddn_Deps *= Ddn_Drn;
      C.addMatrix(Dsigneg_Deps , (1-dn));
      C.addTensorProduct(signeg, Ddn_Deps, -1.0);
    }

    // Positive part
    {
      MatrixND<6,6> Dsigpos_Deps = Qpos*Cbar;
      VectorND<6> Dtaup_Dsigpos{};
      if (taup <= tolp) {
        Dtaup_Dsigpos.zero();
      }
      else  {
        Dtaup_Dsigpos = Se*sigpos;
        // Ce.solve(sigpos, Dtaup_Dsigpos);
        Dtaup_Dsigpos/=taup;
      }

      // Ddp_Deps = Ddp_Drp * Dsigpos_Deps' * Dtaup_Dsigpos;
      VectorND<6> Ddp_Deps = Dsigpos_Deps ^ Dtaup_Dsigpos;
      Ddp_Deps *= Ddp_Drp; 

      C.addMatrix(Dsigpos_Deps , (1-dp));
      C.addTensorProduct(sigpos, Ddp_Deps, -1.0);
    }
  }
  return 0;
}


int
FariaPlasticDamage3d::setTrialStrain(Vector const&v1, Vector const&v2)
{
  return this->setTrialStrain(v1);
}


int
FariaPlasticDamage3d::setTrialStrainIncr(const Vector &strain)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

int
FariaPlasticDamage3d::setTrialStrainIncr(const Vector &strain, const Vector &rate)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

const Matrix&
FariaPlasticDamage3d::getTangent()
{
  return retTangent;
}

const Matrix&
FariaPlasticDamage3d::getInitialTangent()
{
  return retInitialTangent;
}

const Vector&
FariaPlasticDamage3d::getStress()
{
  return retStress;
}

const Vector&
FariaPlasticDamage3d::getStrain()
{
  return retStrain;
}

int
FariaPlasticDamage3d::commitState()
{
  rpCommit    = rp;
  rnCommit    = rn;
  dpCommit    = dp;
  dnCommit    = dn;

  epsCommit   = eps;
  sigCommit   = sig;
  sigeCommit  = sige;
  eps_pCommit = eps_p;

  Ccommit = C;
  return 0;
}

int
FariaPlasticDamage3d::revertToLastCommit()
{
  rp = rpCommit;
  rn = rnCommit;
  dp = dpCommit;
  dn = dnCommit;

  eps = epsCommit;
  sig = sigCommit;
  sige = sigeCommit;
  eps_p = eps_pCommit;
  C  = Ccommit;
  return 0;
}

int
FariaPlasticDamage3d::revertToStart()
{
  eps.zero();
  sig.zero();
  sige.zero();
  eps_p.zero();

  double G  = E/2./(1. +    nu);     // Shear modulus
  double K  = E/3./(1. - 2.*nu);     // Bulk  modulus
  Ce.zero();
  Ce.addMatrix(IIvol, K);
  Ce.addMatrix(IIdevCon, 2.*G);
  // Ce.addMatrix(IIdevMix, 2*G);

  C = Ce;

  double F2c = 1.16*Fc;
  double k = sqrt(2.0)*(F2c - Fc)/(2.*F2c - Fc);

  // initial damage threshold
  double rp0 = ft/sqrt(E);
  double rn0 = sqrt((-k+sqrt(2.0))*Fc/sqrt(3.0));
      
  rp = rp0;
  rn = rn0;
  dp = 0.;
  dn = 0.;
  return 0;
}

NDMaterial*
FariaPlasticDamage3d::getCopy(const char *type)
{
  if (strcmp(type,"ThreeDimensional") == 0 || strcmp(type,"3D") == 0) {
    return this->getCopy();

  } else {
    return NDMaterial::getCopy(type);
  }
}

NDMaterial*
FariaPlasticDamage3d::getCopy()
{
  FariaPlasticDamage3d *theCopy =
    new FariaPlasticDamage3d (this->getTag(), E, nu, ft, Fc, beta, Ap, An, Bn, density);

  return theCopy;
}

const char*
FariaPlasticDamage3d::getType() const
{
  return "ThreeDimensional";
}

int
FariaPlasticDamage3d::getOrder() const
{
  return 6;
}


int 
FariaPlasticDamage3d::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(1+8+4 + 5*6);

  data(0) = this->getTag();

  data(1) = E;
  data(2) = nu;
  data(3) = ft;
  data(4) = Fc;
  data(5) = beta;
  data(6) = Ap;
  data(7) = An;
  data(8) = Bn;

  data(9)  = rpCommit;
  data(10) = rnCommit;
  data(11) = dpCommit;
  data(12) = dnCommit;  

  for (int i = 0; i < 6; i++) {
    data(13+i) = epsCommit(i);
    data(13+6+i) = sigCommit(i);
    data(13+12+i) = sigeCommit(i);
    data(13+18+i) = eps_pCommit(i);
    // data(13+24+i) = sigePCommit(i);
  }

  int res = 0;
  int dbTag = this->getDbTag();

  res = theChannel.sendVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FariaPlasticDamage3d::sendSelf -- could not send Vector\n";
    return res;
  }
  Matrix Ccm(Ccommit);
  res = theChannel.sendMatrix(dbTag, commitTag, Ccm);
  if (res < 0) {
    opserr << "FariaPlasticDamage3d::sendSelf -- could not send Ccommit matrix\n";
    return res;
  }  
  
  return res;
}

int 
FariaPlasticDamage3d::recvSelf(int commitTag, Channel &theChannel, 
				  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dbTag = this->getDbTag();

  static Vector data(43);  
  res = theChannel.recvVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FariaPlasticDamage3d::recvSelf -- could not receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));

  E = data(1);
  nu = data(2);
  ft = data(3);
  Fc = data(4);
  beta = data(5);
  Ap = data(6);
  An = data(7);
  Bn = data(8);

  rpCommit = data(9);
  rnCommit = data(10);
  dpCommit = data(11);
  dnCommit = data(12);

  for (int i = 0; i < 6; i++) {
    epsCommit(i) = data(13+i);
    sigCommit(i) = data(13+6+i);
    sigeCommit(i) = data(13+12+i);
    eps_pCommit(i) = data(13+18+i);
  }
  
  Matrix Ccm(Ccommit);
  res = theChannel.recvMatrix(dbTag, commitTag, Ccm);
  if (res < 0) {
    opserr << "FariaPlasticDamage3d::recvSelf -- could not receive Ccommit matrix\n";
    return res;
  }  

  this->revertToStart();
  this->commitState();
  
  return res;
}

void 
FariaPlasticDamage3d::Print(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"type\": \"" << this->getClassType() << "\", ";
    s << "\"tag\": " << this->getTag() << ", ";
    s << "\"E\": " << E << ", ";
    s << "\"nu\": " << nu << ", ";
    s << "\"Ft\": " << ft << ", ";
    s << "\"Fc\": " << Fc << ", ";
    s << "\"beta\": " << beta << ", ";
    s << "\"Ap\": " << Ap << ", ";
    s << "\"An\": " << An << ", ";
    s << "\"Bn\": " << Bn;
    s << "}";
  }

  else {
    opserr << this->getType() << ": " << this->getTag() << "\n";
    opserr << "stress: " << Vector(retStress) << "\n";
    opserr << "strain: " << Vector(retStrain) << "\n";
    opserr << "tangent: " << Matrix(retTangent) << "\n";
  }
}       
