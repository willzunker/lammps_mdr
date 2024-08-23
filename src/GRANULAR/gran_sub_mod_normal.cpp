/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "gran_sub_mod_normal.h"
#include "error.h"
#include "granular_model.h"
#include "math_const.h"
#include "atom.h"
#include "csv_writer.h"
#include "update.h"

#include <cmath>
#include <iostream>
#include <iomanip> 
#include <sstream>

using namespace LAMMPS_NS;
using namespace Granular_NS;

using MathConst::MY_2PI;
using MathConst::MY_PI;

static constexpr double PI27SQ = 266.47931882941264802866;      // 27*PI**2
static constexpr double THREEROOT3 = 5.19615242270663202362;    // 3*sqrt(3)
static constexpr double SIXROOT6 = 14.69693845669906728801;     // 6*sqrt(6)
static constexpr double INVROOT6 = 0.40824829046386307274;      // 1/sqrt(6)
static constexpr double FOURTHIRDS = (4.0 / 3.0);               // 4/3
static constexpr double JKRPREFIX = 1.2277228507842888;         // cbrt(3*PI**2/16)

/* ----------------------------------------------------------------------
   Default normal model
------------------------------------------------------------------------- */

GranSubModNormal::GranSubModNormal(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp)
{
  material_properties = 0;
  cohesive_flag = 0;
}

/* ---------------------------------------------------------------------- */

bool GranSubModNormal::touch()
{
  bool touchflag = (gm->rsq < gm->radsum * gm->radsum);
  return touchflag;
}

/* ---------------------------------------------------------------------- */

double GranSubModNormal::pulloff_distance(double /*radi*/, double /*radj*/)
{
  // called outside of compute(), do not assume correct geometry defined in contact
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double GranSubModNormal::calculate_contact_radius()
{
  return sqrt(gm->dR);
}

/* ---------------------------------------------------------------------- */

void GranSubModNormal::set_fncrit()
{
  Fncrit = fabs(gm->Fntot);
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModNormalNone::GranSubModNormalNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalNone::calculate_forces()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   Hookean normal force
------------------------------------------------------------------------- */

GranSubModNormalHooke::GranSubModNormalHooke(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 2;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHooke::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hooke normal model");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalHooke::calculate_forces()
{
  return k * gm->delta;
}

/* ----------------------------------------------------------------------
   Hertzian normal force
------------------------------------------------------------------------- */

GranSubModNormalHertz::GranSubModNormalHertz(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 2;
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertz::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hertz normal model");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalHertz::calculate_forces()
{
  return k * gm->contact_radius * gm->delta;
}

/* ----------------------------------------------------------------------
   Hertzian normal force with material properties
------------------------------------------------------------------------- */

GranSubModNormalHertzMaterial::GranSubModNormalHertzMaterial(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormalHertz(gm, lmp)
{
  material_properties = 1;
  num_coeffs = 3;
  contact_radius_flag = 1;
  mixed_coefficients = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertzMaterial::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  if (!mixed_coefficients) {
    if (gm->contact_type == PAIR) {
      k = FOURTHIRDS * mix_stiffnessE(Emod, Emod, poiss, poiss);
    } else {
      k = FOURTHIRDS * mix_stiffnessE_wall(Emod, poiss);
    }
  }

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hertz material normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertzMaterial::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);

  k = FOURTHIRDS * coeffs[0];
  mixed_coefficients = 1;

  coeffs_to_local();
}

/* ----------------------------------------------------------------------
   DMT normal force
------------------------------------------------------------------------- */

GranSubModNormalDMT::GranSubModNormalDMT(GranularModel *gm, LAMMPS *lmp) : GranSubModNormal(gm, lmp)
{
  material_properties = 1;
  cohesive_flag = 1;
  num_coeffs = 4;
  contact_radius_flag = 1;
  mixed_coefficients = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];

  if (!mixed_coefficients) {
    if (gm->contact_type == PAIR) {
      k = FOURTHIRDS * mix_stiffnessE(Emod, Emod, poiss, poiss);
    } else {
      k = FOURTHIRDS * mix_stiffnessE_wall(Emod, poiss);
    }
  }

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal DMT normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);

  k = FOURTHIRDS * coeffs[0];
  mixed_coefficients = 1;

  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalDMT::calculate_forces()
{
  Fne = k * gm->contact_radius * gm->delta;
  F_pulloff = 4.0 * MY_PI * cohesion * gm->Reff;
  Fne -= F_pulloff;
  return Fne;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}

/* ----------------------------------------------------------------------
   JKR normal force
------------------------------------------------------------------------- */

GranSubModNormalJKR::GranSubModNormalJKR(GranularModel *gm, LAMMPS *lmp) : GranSubModNormal(gm, lmp)
{
  material_properties = 1;
  cohesive_flag = 1;
  beyond_contact = 1;
  num_coeffs = 4;
  contact_radius_flag = 1;
  mixed_coefficients = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];

  if (!mixed_coefficients) {
    if (gm->contact_type == PAIR) {
      Emix = mix_stiffnessE(Emod, Emod, poiss, poiss);
    } else {
      Emix = mix_stiffnessE_wall(Emod, poiss);
    }
  }

  k = FOURTHIRDS * Emix;

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal JKR normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);

  Emix = coeffs[0];
  mixed_coefficients = 1;

  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

bool GranSubModNormalJKR::touch()
{
  double delta_pulloff, dist_pulloff;
  bool touchflag;

  if (gm->touch) {
    // delta_pulloff defined as positive so center-to-center separation is > radsum
    delta_pulloff = JKRPREFIX * cbrt(gm->Reff * cohesion * cohesion / (Emix * Emix));
    dist_pulloff = gm->radsum + delta_pulloff;
    touchflag = gm->rsq < (dist_pulloff * dist_pulloff);
  } else {
    touchflag = gm->rsq < (gm->radsum * gm->radsum);
  }

  return touchflag;
}

/* ----------------------------------------------------------------------
  called outside of compute(), do not assume geometry defined in contact
------------------------------------------------------------------------- */

double GranSubModNormalJKR::pulloff_distance(double radi, double radj)
{
  double Reff_tmp;

  Reff_tmp = radi * radj / (radi + radj);    // May not be defined
  if (Reff_tmp <= 0) return 0;
  // Defined as positive so center-to-center separation is > radsum
  return JKRPREFIX * cbrt(Reff_tmp * cohesion * cohesion / (Emix * Emix));
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalJKR::calculate_contact_radius()
{
  double R2, dR2, t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3;

  R2 = gm->Reff * gm->Reff;
  dR2 = gm->dR * gm->dR;
  t0 = cohesion * cohesion * R2 * R2 * Emix;
  t1 = PI27SQ * t0;
  t2 = 8.0 * gm->dR * dR2 * Emix * Emix * Emix;
  t3 = 4.0 * dR2 * Emix;

  // in case sqrt(0) < 0 due to precision issues
  sqrt1 = MAX(0, t0 * (t1 + 2.0 * t2));
  t4 = cbrt(t1 + t2 + THREEROOT3 * MY_PI * sqrt(sqrt1));
  t5 = t3 / t4 + t4 / Emix;
  sqrt2 = MAX(0, 2.0 * gm->dR + t5);
  t6 = sqrt(sqrt2);
  sqrt3 = MAX(0, 4.0 * gm->dR - t5 + SIXROOT6 * cohesion * MY_PI * R2 / (Emix * t6));

  return INVROOT6 * (t6 + sqrt(sqrt3));
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalJKR::calculate_forces()
{
  double a2;
  a2 = gm->contact_radius * gm->contact_radius;
  Fne = k * gm->contact_radius * a2 / gm->Reff -
      MY_2PI * a2 * sqrt(4.0 * cohesion * Emix / (MY_PI * gm->contact_radius));
  F_pulloff = 3.0 * MY_PI * cohesion * gm->Reff;

  return Fne;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}

/* ----------------------------------------------------------------------
   MDR contact model

   Contributing authors:
   William Zunker (MIT), Sachith Dunatunga (MIT),
   Dan Bolintineanu (SNL), Joel Clemmer (SNL)
------------------------------------------------------------------------- */

GranSubModNormalMDR::GranSubModNormalMDR(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 6; // Young's Modulus, Poisson's ratio, yield stress, effective surface energy, psi_b, coefficent of restitution
  contact_radius_flag = 1;
  size_history = 23; 

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  //transfer_history_factor[0] = +1;
  for (int i = 0; i < size_history; i++) { 
    transfer_history_factor[i] = +1;
  }
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalMDR::coeffs_to_local()
{
  E = coeffs[0];      // Young's modulus
  nu = coeffs[1];     // Poisson's ratio
  Y = coeffs[2];      // yield stress
  gamma = coeffs[3];  // effective surface energy
  psi_b = coeffs[4];  // bulk response trigger based on ratio of remaining free area: A_{free}/A_{total}
  CoR = coeffs[5];    // coefficent of restitution

  if (E <= 0.0) error->all(FLERR, "Illegal MDR normal model, Young's modulus must be greater than 0");
  if (nu < 0.0 || nu > 0.5) error->all(FLERR, "Illegal MDR normal model, Poisson's ratio must be between 0 and 0.5");
  if (Y < 0.0) error->all(FLERR, "Illegal MDR normal model, yield stress must be greater than or equal to 0");
  if (gamma < 0.0) error->all(FLERR, "Illegal MDR normal model, effective surface energy must be greater than or equal to 0");
  if (psi_b < 0.0 || psi_b > 1.0) error->all(FLERR, "Illegal MDR normal model, psi_b must be between 0 and 1.0");
  if (CoR < 0.0 || CoR > 1.0) error->all(FLERR, "Illegal MDR normal model, coefficent of restitution must be between 0 and 1.0");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalMDR::calculate_forces()

{
  // To understand the structure of the overall code it is important to consider 
  // the following:
  //
  // The MDR contact model was developed by imagining individual particles being 
  // squished between a number of rigid flats (references below). To allow  
  // for many interacting particles, we extend the idea of isolated particles surrounded 
  // by rigid flats. In particular, we imagine placing rigid flats at the overlap 
  // midpoints between particles. The force is calculated seperately on both sides
  // of the contact assuming interaction with a rigid flat. The two forces are then 
  // averaged on either side of the contact to determine the final force. If the 
  // contact is between a particle and wall then only one force evaluation is required.
  //  
  // Zunker and Kamrin, 2024, Part I: https://doi.org/10.1016/j.jmps.2023.105492
  // Zunker and Kamrin, 2024, Part II: https://doi.org/10.1016/j.jmps.2023.105493

  const int i_true = gm->i;               // true i particle index
  const int j_true = gm->j;               // true j particle index
  const double radi_true = gm->radi;      // true i particle initial radius
  const double radj_true = gm->radj;      // true j particle initial radius
    
  double F;                               // average force 
  double F0;                              // force on contact side 0
  double F1;                              // force on contact side 1
  double delta = gm->delta;               // apparent overlap
  if (gm->contact_type == PAIR) delta = gm->delta/2.0; // half displacement to imagine interaction with rigid flat 

  // initialize indexing in history array of different constact history variables 
  const int delta_offset_0 = 0;           // apparent overlap 
  const int delta_offset_1 = 1;           
  const int deltao_offset_0 = 2;          // displacement 
  const int deltao_offset_1 = 3;  
  const int delta_MDR_offset_0 = 4;       // MDR apparent overlap
  const int delta_MDR_offset_1 = 5;
  const int delta_BULK_offset_0 = 6;      // bulk displacement
  const int delta_BULK_offset_1 = 7;
  const int deltamax_MDR_offset_0 = 8;    // maximum MDR apparent overlap 
  const int deltamax_MDR_offset_1 = 9;
  const int Yflag_offset_0 = 10;          // yield flag
  const int Yflag_offset_1 = 11;
  const int deltaY_offset_0 = 12;         // yield displacement
  const int deltaY_offset_1 = 13;
  const int cA_offset_0 = 14;             // contact area intercept
  const int cA_offset_1 = 15;
  const int aAdh_offset_0 = 16;           // adhesive contact radius
  const int aAdh_offset_1 = 17;
  const int Ac_offset_0 = 18;             // contact area
  const int Ac_offset_1 = 19;
  const int eps_bar_offset_0 = 20;        // volume-averaged infinitesimal strain tensor
  const int eps_bar_offset_1 = 21;
  const int penalty_offset_ = 22;         // contact penalty 

  // initialize particle history variables 
  int tmp1, tmp2;
  int index_Ro = atom->find_custom("Ro",tmp1,tmp2);                       // initial radius
  int index_Vcaps = atom->find_custom("Vcaps",tmp1,tmp2);                 // spherical cap volume from intersection of apparent radius particle and contact planes
  int index_Vgeo = atom->find_custom("Vgeo",tmp1,tmp2);                   // geometric particle volume of apparent particle after removing spherical cap volume
  int index_Velas = atom->find_custom("Velas",tmp1,tmp2);                 // particle volume from linear elasticity  
  int index_eps_bar = atom->find_custom("eps_bar",tmp1,tmp2);             // volume-averaged infinitesimal strain tensor
  int index_dRnumerator = atom->find_custom("dRnumerator",tmp1,tmp2);     // summation of numerator terms in calculation of dR
  int index_dRdenominator = atom->find_custom("dRdenominator",tmp1,tmp2); // summation of denominator terms in calculation of dR
  int index_Acon0 = atom->find_custom("Acon0",tmp1,tmp2);                 // total area involved in contacts: Acon^{n} 
  int index_Acon1 = atom->find_custom("Acon1",tmp1,tmp2);                 // total area involved in contacts: Acon^{n+1}
  int index_Atot = atom->find_custom("Atot",tmp1,tmp2);                   // total particle area 
  int index_Atot_sum = atom->find_custom("Atot_sum",tmp1,tmp2);           // running sum of contact area minus cap area
  int index_ddelta_bar = atom->find_custom("ddelta_bar",tmp1,tmp2);       // change in mean surface displacement
  int index_psi = atom->find_custom("psi",tmp1,tmp2);                     // ratio of free surface area to total surface area
  int index_sigmaxx = atom->find_custom("sigmaxx",tmp1,tmp2);             // xx-component of the stress tensor, not necessary for force calculation
  int index_sigmayy = atom->find_custom("sigmayy",tmp1,tmp2);             // yy-component of the stress tensor, not necessary for force calculation  
  int index_sigmazz = atom->find_custom("sigmazz",tmp1,tmp2);             // zz-component of the stress tensor, not necessary for force calculation 
  int index_contacts = atom->find_custom("contacts",tmp1,tmp2);                     // total contacts on particle 
  int index_adhesive_length = atom->find_custom("adhesive_length",tmp1,tmp2);       // total contacts on particle 
  double * Rinitial = atom->dvector[index_Ro];
  double * Vgeo = atom->dvector[index_Vgeo];
  double * Velas = atom->dvector[index_Velas];
  double * Vcaps = atom->dvector[index_Vcaps];
  double * eps_bar = atom->dvector[index_eps_bar];
  double * dRnumerator = atom->dvector[index_dRnumerator];
  double * dRdenominator = atom->dvector[index_dRdenominator];
  double * Acon0 = atom->dvector[index_Acon0];
  double * Acon1 = atom->dvector[index_Acon1]; 
  double * Atot = atom->dvector[index_Atot];
  double * Atot_sum = atom->dvector[index_Atot_sum];
  double * ddelta_bar = atom->dvector[index_ddelta_bar];
  double * psi = atom->dvector[index_psi];
  double * sigmaxx = atom->dvector[index_sigmaxx];
  double * sigmayy = atom->dvector[index_sigmayy];
  double * sigmazz = atom->dvector[index_sigmazz];
  double * contacts = atom->dvector[index_contacts];
  double * adhesive_length = atom->dvector[index_adhesive_length];
  
  double * history = & gm->history[history_index]; // load in all history variables  

  for (int contactSide = 0; contactSide < 2; contactSide++) { 

    double * delta_offset; 
    double * deltao_offset;
    double * delta_MDR_offset;   
    double * delta_BULK_offset; 
    double * deltamax_MDR_offset; 
    double * Yflag_offset; 
    double * deltaY_offset; 
    double * cA_offset;
    double * aAdh_offset; 
    double * Ac_offset; 
    double * eps_bar_offset; 
    double * penalty_offset;

    if (contactSide == 0) {
      if (gm->contact_type == PAIR) {
        gm->i = std::max(i_true,j_true);
        gm->j = std::min(i_true,j_true);
        if (gm->i == i_true) {
          gm->radi = radi_true;
          gm->radj = radj_true;
        } else {
          gm->radi = radj_true;
          gm->radj = radi_true;
        }
      }
      delta_offset = & history[delta_offset_0];
      deltao_offset = & history[deltao_offset_0];
      delta_MDR_offset = & history[delta_MDR_offset_0];
      delta_BULK_offset = & history[delta_BULK_offset_0];
      deltamax_MDR_offset = & history[deltamax_MDR_offset_0];
      Yflag_offset = & history[Yflag_offset_0];
      deltaY_offset = & history[deltaY_offset_0];
      cA_offset = & history[cA_offset_0];
      aAdh_offset = & history[aAdh_offset_0];
      Ac_offset = & history[Ac_offset_0];
      eps_bar_offset = & history[eps_bar_offset_0];
    } else {
      if (gm->contact_type != PAIR) break; // contact with particle-wall requires only one evaluation
      gm->i = std::min(i_true,j_true);
      gm->j = std::max(i_true,j_true);
      if (gm->i == i_true) {
        gm->radi = radi_true;
        gm->radj = radj_true;
      } else {
        gm->radi = radj_true;
        gm->radj = radi_true;
      }
      delta_offset = & history[delta_offset_1];
      deltao_offset = & history[deltao_offset_1];
      delta_MDR_offset = & history[delta_MDR_offset_1];
      delta_BULK_offset = & history[delta_BULK_offset_1];
      deltamax_MDR_offset = & history[deltamax_MDR_offset_1];
      Yflag_offset = & history[Yflag_offset_1];
      deltaY_offset = & history[deltaY_offset_1];
      cA_offset = & history[cA_offset_1];
      aAdh_offset = & history[aAdh_offset_1];
      Ac_offset = & history[Ac_offset_1];
      eps_bar_offset = & history[eps_bar_offset_1];
    }

    // temporary i and j indices
    const int i = gm->i;
    const int j = gm->j;

    // material and geometric property definitions
    // E, nu, Y gamma , psi_b, and CoR are already defined.
    const double G = E/(2.0*(1.0+nu));          // shear modulus
    const double kappa = E/(3.0*(1.0-2.0*nu));  // bulk modulus
    const double Eeff = E/(1.0-pow(nu,2.0));    // composite plane strain modulus

    const double Ro = Rinitial[i];              // initial radius
    const double R = gm->radi;                  // apparent radius
 
    // kinematics 
    const double ddelta = delta - *delta_offset;
    *delta_offset = delta;

    const double deltao = delta - (R - Ro);
    const double ddeltao = deltao - *deltao_offset;
    *deltao_offset = deltao;

    double ddelta_MDR;
    double ddelta_BULK;
    if ( psi[i] < psi_b ) { // if true, bulk response has triggered, split displacement increment between the MDR and BULK components 
      ddelta_MDR = std::min(ddelta-ddelta_bar[i], delta-*delta_MDR_offset);
      ddelta_BULK = ddelta_bar[i];
    } else { // if false, no bulk response, full displacement increment goes to the MDR component
      ddelta_BULK = 0.0;
      ddelta_MDR = ddelta;
    }
    const double delta_MDR = *delta_MDR_offset + ddelta_MDR; // MDR displacement
    *delta_MDR_offset = delta_MDR; // Update old MDR displacement
    const double delta_BULK = std::max(0.0,*delta_BULK_offset+ddelta_BULK); // bulk displacement
    *delta_BULK_offset = delta_BULK; // update old bulk displacement

    if (delta_MDR > *deltamax_MDR_offset) *deltamax_MDR_offset = delta_MDR;
    const double deltamax_MDR = *deltamax_MDR_offset;

    const double pY = Y*(1.75*exp(-4.4*deltamax_MDR/R) + 1.0); // Set value of average pressure along yield surface
    if ( *Yflag_offset == 0.0 && delta_MDR >= deltamax_MDR ) {
    const double phertz = 4*Eeff*sqrt(delta_MDR)/(3*M_PI*sqrt(R));
      if ( phertz > pY ) {
        *Yflag_offset = 1.0;
        *deltaY_offset = delta_MDR;
        *cA_offset = M_PI*(pow(*deltaY_offset,2.0) - *deltaY_offset*R);
      }
    } 

    // MDR force calculation
    double F_MDR;
    double A;                     // height of elliptical indenter
    double B;                     // width of elliptical indenter
    double deltae1D;              // transformed elastic displacement
    double deltaR;                // displacement correction 
    double amax;                  // maximum experienced contact radius
    const double cA = *cA_offset; // contact area intercept

    if ( *Yflag_offset == 0.0 ) { // elastic contact
      A = 4.0*R;              
      B = 2.0*R;              
      deltae1D = delta_MDR;    
      (deltae1D > 0) ?  amax = sqrt(deltae1D*R) : amax = 0.0;  
    } else { // plastic contact
      amax = sqrt((2.0*deltamax_MDR*R - pow(deltamax_MDR,2.0)) + cA/M_PI);              
      A = 4.0*pY/Eeff*amax;                                                             
      B = 2.0*amax;                                                                     
      const double deltae1Dmax = A/2.0; // maximum transformed elastic displacement 
      const double Fmax = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltae1Dmax/A) - (1.0 - 2.0*deltae1Dmax/A)*sqrt(4.0*deltae1Dmax/A - 4.0*pow(deltae1Dmax,2.0)/pow(A,2.0))); // force caused by full submersion of elliptical indenter to depth of A/2
      const double zR = R - (deltamax_MDR - deltae1Dmax); // depth of particle center
      deltaR = (Fmax*(2*pow(amax,2.0)*(-1 + nu) - (-1 + 2*nu)*zR*(-zR + sqrt(pow(amax,2.0) + pow(zR,2.0)))))/((M_PI*pow(amax,2.0))*2*G*sqrt(pow(amax,2.0) + pow(zR,2.0))); 
      deltae1D = (delta_MDR - deltamax_MDR + deltae1Dmax + deltaR)/(1 + deltaR/deltae1Dmax);  // transformed elastic displacement 
    }

    double a_na;
    (deltae1D >= 0.0) ? a_na = B*sqrt(A - deltae1D)*sqrt(deltae1D)/A : a_na = 0.0;
    double aAdh = *aAdh_offset; 
    if ( gamma > 0.0  ) { // adhesive contact
    double g_aAdh;
    
      if (delta_MDR == deltamax_MDR || a_na >= aAdh ) { // case 1: no tensile springs, purely compressive contact
        (deltae1D <= 0.0) ? F_MDR = 0.0 : F_MDR = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltae1D/A) - (1.0 - 2.0*deltae1D/A)*sqrt(4.0*deltae1D/A - 4.0*pow(deltae1D,2.0)/pow(A,2.0))); 
        if ( std::isnan(F_MDR) ) std::cout << "F_MDR is NaN, case 1: no tensile springs" << std::endl;
        *aAdh_offset = a_na;
      } else {
        const double lmax = sqrt(2.0*M_PI*aAdh*gamma/Eeff); 
        g_aAdh = A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh,2.0));
        const double acrit = (-((pow(B,2)*gamma*M_PI)/(pow(A,2)*Eeff)) + (pow(2,0.3333333333333333)*pow(B,4)*pow(gamma,2)*pow(M_PI,1.6666666666666667))/
                            (pow(A,2)*pow(Eeff,2)*pow((27*pow(A,4)*pow(B,4)*gamma)/Eeff - (2*pow(B,6)*pow(gamma,3)*pow(M_PI,2))/pow(Eeff,3) + (3*sqrt(3)*sqrt(27*pow(A,8)*pow(B,8)*pow(Eeff,2)*pow(gamma,2) -
                              4*pow(A,4)*pow(B,10)*pow(gamma,4)*pow(M_PI,2)))/pow(Eeff,2),0.3333333333333333)) + (pow(M_PI/2.,0.3333333333333333)*pow((27*pow(A,4)*pow(B,4)*gamma)/Eeff -
                            (2*pow(B,6)*pow(gamma,3)*pow(M_PI,2))/pow(Eeff,3) + (3*sqrt(3)*sqrt(27*pow(A,8)*pow(B,8)*pow(Eeff,2)*pow(gamma,2) - 4*pow(A,4)*pow(B,10)*pow(gamma,4)*pow(M_PI,2)))/
                            pow(Eeff,2),0.3333333333333333))/pow(A,2))/6;

        if ( deltae1D + lmax - g_aAdh >= 0) { // case 2: tensile springs, but not exceeding critical length --> deltae + lmax - g(aAdhes) >= 0
          const double deltaeAdh = g_aAdh; 
          const double F_na = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltaeAdh/A) - (1.0 - 2.0*deltaeAdh/A)*sqrt(4.0*deltaeAdh/A - 4.0*pow(deltaeAdh,2.0)/pow(A,2.0)));
          const double F_Adhes = 2.0*Eeff*(deltae1D - deltaeAdh)*aAdh;
          F_MDR = F_na + F_Adhes; 
          if ( std::isnan(F_MDR) ) std::cout << "F_MDR is NaN, case 2: tensile springs, but not exceeding critical length" << std::endl;

        } else { // case 3: tensile springs exceed critical length --> deltae + lmax - g(aAdhes) = 0
          // newton-raphson to find aAdh
          const double maxIterations = 100;
          const double error = 1e-16;
          double aAdh_tmp = aAdh;
          double fa; 
          double fa2;
          double dfda;
          for (int lv1 = 0; lv1 < maxIterations; ++lv1) {
            fa = deltae1D + sqrt(2.0*M_PI*aAdh_tmp*gamma/Eeff) - ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh_tmp,2.0)) );
            if (abs(fa) < error) {
              break;
            } 
            dfda = -((aAdh_tmp*A)/(B*sqrt(-pow(aAdh_tmp,2.0) + pow(B,2.0)/4.0))) + (gamma*sqrt(M_PI/2.0))/(Eeff*sqrt((aAdh_tmp*gamma)/Eeff));
            aAdh_tmp = aAdh_tmp - fa/dfda;
            fa2 = deltae1D + sqrt(2.0*M_PI*aAdh_tmp*gamma/Eeff) - ( A/2 - A/B*sqrt(pow(B,2.0)/4 - pow(aAdh_tmp,2.0)) );
            if (abs(fa-fa2) < error) {
              break;
            } 
            if (lv1 == maxIterations-1){
              aAdh_tmp = 0.0;
            }
          }
          aAdh = aAdh_tmp;

          if ( aAdh < acrit ) {
            F_MDR = 0.0;
          } else {
            g_aAdh = A/2.0 - A/B*sqrt(pow(B,2.0)/4.0 - pow(aAdh,2.0));                  
            const double deltaeAdh = g_aAdh; 
            const double F_na = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltaeAdh/A) - (1.0 - 2.0*deltaeAdh/A)*sqrt(4.0*deltaeAdh/A - 4.0*pow(deltaeAdh,2.0)/pow(A,2.0)));
            const double F_Adhes = 2.0*Eeff*(deltae1D - deltaeAdh)*aAdh;
            F_MDR = F_na + F_Adhes; 
            if ( std::isnan(F_MDR) ) std::cout << "F_MDR is NaN, case 3: tensile springs exceed critical length" << std::endl;
          }
          *aAdh_offset = aAdh;
        }
      }
    } else { // non-adhesive contact
      *aAdh_offset = a_na;
      (deltae1D <= 0.0) ? F_MDR = 0.0 : F_MDR = Eeff*(A*B/4.0)*(acos(1.0 - 2.0*deltae1D/A) - (1.0 - 2.0*deltae1D/A)*sqrt(4.0*deltae1D/A - 4.0*pow(deltae1D,2.0)/pow(A,2.0))); 
      if ( std::isnan(F_MDR) ) std::cout << "F_MDR is NaN, non-adhesive case" << std::endl;
    }
    
    contacts[i] += 1;
    adhesive_length[i] += aAdh;

    // contact penalty scheme
    penalty_offset = & history[penalty_offset_];
    double pij = *penalty_offset;
    const double wij = std::max(1.0-pij,0.0);

    // area related calculations 
    double Ac; 
    (*Yflag_offset == 0.0) ? Ac = M_PI*delta*R : Ac = M_PI*((2.0*delta*R - pow(delta,2.0)) + cA/M_PI);
    if (Ac < 0.0 ) Ac = 0.0;
    Atot_sum[i] += wij*(Ac - 2.0*M_PI*R*(deltamax_MDR + delta_BULK));
    Acon1[i] += wij*Ac;

    // bulk force calculation
    double F_BULK;
    (delta_BULK <= 0.0) ? F_BULK = 0.0 : F_BULK = (1.0/Vgeo[i])*Acon0[i]*delta_BULK*kappa*Ac;

    // total force calculation
    (contactSide == 0) ? F0 = F_MDR + F_BULK : F1 = F_MDR + F_BULK;

    // mean surface dipslacement calculation
     *Ac_offset = wij*Ac;

    // radius update scheme quantity calculation
    Vcaps[i] += (M_PI/3.0)*pow(delta,2.0)*(3.0*R - delta);
    
    const double Fntmp = wij*(F_MDR + F_BULK);
    const double fx = Fntmp*gm->nx[0];
    const double fy = Fntmp*gm->nx[1];
    const double fz = Fntmp*gm->nx[2];
    const double bx = -(Ro - deltao)*gm->nx[0];
    const double by = -(Ro - deltao)*gm->nx[1];
    const double bz = -(Ro - deltao)*gm->nx[2];
    const double eps_bar_contact = (1.0/(3*kappa*Velas[i]))*(fx*bx + fy*by + fz*bz);
    eps_bar[i] += eps_bar_contact;
    
    if(delta_MDR == deltamax_MDR && *Yflag_offset > 0.0 && F_MDR > 0.0){
      const double Vo = (4.0/3.0)*M_PI*pow(Ro,3.0);
      dRnumerator[i] += -Vo*(eps_bar_contact - *eps_bar_offset) - wij*M_PI*ddeltao*( 2.0*deltao*Ro - pow(deltao,2.0) + pow(R,2.0) - pow(Ro,2.0) );

      dRdenominator[i] += wij*2.0*M_PI*R*(deltao + R - Ro);
    }
    *eps_bar_offset = eps_bar_contact;

    sigmaxx[i] += (1.0/Velas[i])*(fx*bx);
    sigmayy[i] += (1.0/Velas[i])*(fy*by);
    sigmazz[i] += (1.0/Velas[i])*(fz*bz);
  }

  // rigid flat placement scheme 
  gm->i = i_true;
  gm->j = j_true;
  gm->radi = radi_true;
  gm->radj = radj_true;

  double * penalty_offset = & history[penalty_offset_];
  const double pij = *penalty_offset;
  const double wij = std::max(1.0-pij,0.0);
  *penalty_offset = 0.0;

  // wall force magnifier
  double * deltao_offset = & history[deltao_offset_0];
  const double wallForceMagnifer = 1.0;

  // assign final force
  (gm->contact_type != PAIR) ? F = wij*F0*wallForceMagnifer : F = wij*(F0 + F1)/2; 

  // calculate damping force
  if (F > 0.0) {
    double Eeff;
    double Reff;
    if (gm->contact_type == PAIR) {
      Eeff = E/(2.0*(1.0-pow(nu,2.0)));
      Reff = pow((1/gm->radi + 1/gm->radj),-1);
    } else {
      Eeff = E/(1.0-pow(nu,2.0));
      Reff = gm->radi;
    }
    const double kn = Eeff*Reff;
    const double beta = -log(CoR)/sqrt(pow(log(CoR),2.0) + M_PI*M_PI);
    const double damp_prefactor = beta*sqrt(gm->meff*kn);
    const double F_DAMP = -damp_prefactor*(gm->vnnr);

    F += wij*F_DAMP;
  }

  return F;
}