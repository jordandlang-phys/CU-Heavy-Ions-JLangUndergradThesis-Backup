//#if !defined(MYLIB_ML_PROJECT_HEADER_H)
//#define MYLIB_ML_PROJECT_HEADER_H

//  File Name : jet_ml_constants.h
namespace Jet_ML_Constants {

/*
    The default values below are intended to make managing projects easier.
    These can be changed here, but can also be overwritten when calling any function in Jet_ML_Generator_PYTHIA.cpp
 */

// Default Physical Values
const float math_pi = 3.14159265359;
const float math_e  = 2.71828182846;
const float m_pion  = 0.1396; // [GeV] Pion mass

// PYTHIA Default Values
const float d_beam_power            = 2760.;    // [GeV] Simulated beam power in the dector. Default is 2.76 TeV for ALICE.
const float d_detector_eta          = 0.9;      // Maximum rapidity of the detector. Default is 0.9 for ALICE.
const bool  d_use_voronoi           = false;    // Use Voronoi for jet area determination
const int   d_accept_jet_index      = 0;        // If =0, accepts events with any jet between min/max. If >0, only accepts events with top N jets between min/max

// Thermal Default Values
const int   d_thermal_mean          = 1800;     // Mean number of thermal particles
const int   d_thermal_sigma         = 200;      // Standard deviation for number of thermal particles
const float d_thermal_pt_max        = 100;      // Max pT of a thermal particle (adjust based on pT distribution)
const float d_MH_par_1              = 64547;    // Modified Hagedorn function, parameter 1
const float d_MH_par_2              = 3.076;    // Modified Hagedorn function, parameter 2
const float d_MH_par_3              = 1.126;    // Modified Hagedorn function, parameter 3
const float d_MH_par_4              = -8.491;   // Modified Hagedorn function, parameter 4

// FastJet Default Values
const float d_jet_eta_max           = 0.5;      // Largest jet rapidity to consider.
const float d_pt_hat_min            = 7.5;      // [GeV] ptHatMin value.
const float d_pt_hat_max            = 0.;       // [GeV] ptHatMax value. If 0. then no upper limit.
const float d_fastjet_pt_min        = 8.0;      // [GeV] Minimum pT considered by FastJet for a jet.
const float d_fastjet_radius        = 0.4;      // Jet radius used by FastJet.
const float d_fastjet_match_radius  = 0.3;      // Matches combined jets to PYTHIA jets within this radius.
const float d_jet_pt_raw_max        = 200.;     // [GeV] Maximum raw jet pt allowed

// ML Prep Default Values
const int   d_softest_jet_index     = 0;        // If =0, tries to match with all PYTHIA jets. If >0, matches only to top N PYTHIA jets.

// Nice Plot Colors
const int plot_black    = kGray+3;
const int plot_red      = kPink;
const int plot_green    = kTeal-6;
const int plot_blue     = kAzure;
const int plot_violet   = kViolet-1;

// Marker Styles
const float mark_circ_open[2] = {24, 1.0};
const float mark_circ_fill[2] = {20, 1.0};
const float mark_squa_open[2] = {25, 0.9};
const float mark_squa_fill[2] = {21, 0.9};
const float mark_diam_open[2] = {27, 1.5};
const float mark_diam_fill[2] = {33, 1.5};
const float mark_star_open[2] = {42, 1.4};
const float mark_star_fill[2] = {43, 1.4};
const float mark_plus_open[2] = {28, 0.9};
const float mark_plus_fill[2] = {34, 0.9};
const float mark_triu_open[2] = {26, 1.1};
const float mark_triu_fill[2] = {22, 1.1};
const float mark_trid_open[2] = {32, 1.1};
const float mark_trid_fill[2] = {23, 1.1};

} // End Jet_ML_Constants namespace
