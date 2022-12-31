//#if !defined(MYLIB_ML_PROJECT_HEADER_H)
//#define MYLIB_ML_PROJECT_HEADER_H

//  File Name : jet_ml_constants.h
namespace Jet_ML_Constants {

// Math Default Values
const double math_pi        = 3.14159265359;
const double math_e         = 2.71828182846;

// PYTHIA Default Values
const float  beam_power     = 2760.;    // [GeV] Simulated beam power in the dector. Default is 2.76 TeV for ALICE.
const float  detector_eta   = 0.9;      // Maximum rapidity of the detector. Default is 0.9 for ALICE.
const bool   use_voronoi    = false;    // Use Voronoi for jet clustering

// Thermal Default Values
const int    gaus_mean      = 1800;     // Mean number of thermal particles
const int    gaus_sigma     = 200;      // Standard deviation for number of thermal particles
const double gaus_norm      = 1 / sqrt(2 * math_pi * pow(gaus_sigma,2));    // Normalization factor
const double m_pion         = 0.1396;   // [GeV] Pion mass

// Nice Plot Colors
const int plot_black        = kGray+3;
const int plot_red          = kPink;
const int plot_green        = kTeal-6;
const int plot_blue         = kAzure;
const int plot_violet       = kViolet-1;

// Marker Styles
const double mark_circ_open[2] = {24, 1.0};
const double mark_circ_fill[2] = {20, 1.0};
const double mark_squa_open[2] = {25, 0.9};
const double mark_squa_fill[2] = {21, 0.9};
const double mark_diam_open[2] = {27, 1.5};
const double mark_diam_fill[2] = {33, 1.5};
const double mark_star_open[2] = {42, 1.4};
const double mark_star_fill[2] = {43, 1.4};
const double mark_plus_open[2] = {28, 0.9};
const double mark_plus_fill[2] = {34, 0.9};
const double mark_triu_open[2] = {26, 1.1};
const double mark_triu_fill[2] = {22, 1.1};
const double mark_trid_open[2] = {32, 1.1};
const double mark_trid_fill[2] = {23, 1.1};

} // End Jet_ML_Constants namespace
