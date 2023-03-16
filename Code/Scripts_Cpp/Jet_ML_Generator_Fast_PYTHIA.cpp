#include "Pythia8/Pythia.h"
using namespace Pythia8;
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/Selector.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom3.h"
#include <cmath>
#include <iostream>
#include <filesystem>

#include "jet_ml_constants.h" // Contains default values for many parameters.
using namespace Jet_ML_Constants;

const bool print_out        = true;
const bool debug            = false;
const int  print_every_x    = 1000;

// Sets functions within the Jet_Generator namespace
namespace Jet_Generator {


// Fit functions for thermal generator:
double Func_Gaussian(double* x, double* par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    double fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}
double Func_Power(double* x, double* par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[2]);
    double fitval = par[0] * pow(arg, par[1]);
    return fitval;
}
double Func_Flat(double* x, double* par) {
    double fitval = par[0] * 1;
    return fitval;
}
double Func_ModifiedHagedorn(double* x, double* par) {
    double arg1 = par[0] * pow(x[0], 2) / pow( (pow(x[0], 2) + pow(par[1], 2)), 0.5);
    double arg2 = pow( 1 + (x[0] /par[2]), par[3] );
    double arg3 = 1 / x[0];
    double fitval =  arg1 * arg2 * arg3;
    return fitval;
}

// Random generator functions for thermal generator:
int    Generate_N(TF1* thermal_n_func, float gaus_mean) {
    int n = round( thermal_n_func->GetRandom(0, 2 * gaus_mean) );
    return n;
}
double Generate_Pt(TF1* thermal_pt_func, float thermal_pt_max) {
    double pt = thermal_pt_func->GetRandom(0, thermal_pt_max);
    return pt;
}
double Generate_Eta(TF1* thermal_eta_func, float detector_eta) {
    double eta = thermal_eta_func->GetRandom(-1 * detector_eta, detector_eta);
    return eta;
}
double Generate_Phi(TF1* thermal_phi_func) {
    double phi = thermal_phi_func->GetRandom(-1 * math_pi, math_pi);
    return phi;
}

/*
    Full_Jet_Generator()
    Streamlines the jet generation process by combining PYTHIA and thermal particle generation and jet clustering into a single function
 */
void Jet_Generator_Optimized(
    const char  output_base_name[100],                  // Base/stem of output file name (no file extensions or prefixes).
    const char  output_directory[400],                  // Path to desired output directory.
    const int   event_count,                            // Number of collision events to generate. The number of jets will be greater than this.
    const float pythia_ptbias,                          // Bias applied to the particle pT distribution.
    const float jet_pt_min,                             // Minimum jet pT to accept. At least 1 jet per event will be greater than this.
    const float jet_pt_max,                             // Maximum jet pT to accept. At least 1 jet per event will be less than this.
    const float jet_eta_max = d_jet_eta_max,            // Largest jet rapidity to consider.
    const float pythia_pthatmin = d_pt_hat_min,         // [GeV] PYTHIA ptHatMin value.
    const float pythia_pthatmax = d_pt_hat_max,         // [GeV] PYTHIA ptHatMax value. If 0. then no upper limit.
    const float fastjet_radius = d_fastjet_radius,      // Jet radius used by FastJet.
    const float fastjet_pt_min = d_fastjet_pt_min,      // [GeV] Minimum pT considered by FastJet for a jet.
    const float fastjet_match_radius = d_fastjet_match_radius,  // Matches combined jets to PYTHIA jets within this radius.
    const int   fastjet_accept_index = d_accept_jet_index,      // If =0, accepts events with any jet in min/max. If >0, only accepts if one of top N jets in min/max.
    const int   ml_soft_jet_index = d_softest_jet_index,        // If =0, tries to match with all PYTHIA jets. If >0, matches only to top N PYTHIA jets.
    const int   thermal_mean = d_thermal_mean,          // Mean number of thermal particles
    const int   thermal_sigma = d_thermal_sigma,        // Standard deviation for number of thermal particles
    const float pythia_beampower = d_beam_power,        // [GeV] Simulated beam power in the dector. Default defined by jet_ml_constants.h file.
    const float detector_eta = d_detector_eta,          // Maximum rapidity of the detector. Default defined by jet_ml_constants.h file.
    const float jet_pt_raw_max = d_jet_pt_raw_max,      // Max raw jet pT to accept
    const bool  export_thermal = false,                 // If true, exports thermal particles to root file
    const float thermal_pt_max = d_thermal_pt_max,      // Max thermal pT (typically twice the thermal mean)
    const float thermal_MH_par_1 = d_MH_par_1,          // Modified Hagedorn function parameter 1
    const float thermal_MH_par_2 = d_MH_par_2,          // Modified Hagedorn function parameter 2
    const float thermal_MH_par_3 = d_MH_par_3,          // Modified Hagedorn function parameter 3
    const float thermal_MH_par_4 = d_MH_par_4           // Modified Hagedorn function parameter 4
    ) {
    
    std::cout << "\nStarting Full_Jet_Generator()...\n" << std::endl;
    
    // Builds output file and TTrees
    char output_file_params[200];
    snprintf(output_file_params, 200, "%s_B%.0f_%.0f_%.0f_N%i", output_base_name, pythia_ptbias, jet_pt_min, jet_pt_max, event_count);
    char output_file_path[500];
    snprintf(output_file_path, 500, "%s/Full_%s.root", output_directory, output_file_params);
    TFile* output_file      = new TFile(output_file_path, "UPDATE");
    TTree* gen_param_tree   = new TTree("Generator_Parameters",     "Flat tree of parameters used for jet generation");
    TTree* pyth_par_tree    = new TTree("PYTHIA_Particle_Tree",     "Tree of particles from PYTHIA p+p collisions");
    TTree* pyth_jet_tree    = new TTree("PYTHIA_Jet_Tree",          "Tree of jet clusters from PYTHIA, by event");
    TTree* comb_par_tree    = new TTree("Combined_Particle_Tree",   "Tree of particles combined from PYTHIA and thermal model generation");
    TTree* comb_jet_tree    = new TTree("Combined_Jet_Tree",        "Tree of jet clusters from combined events, by event");
    TTree* ther_par_tree;
    
    std::cout << "Output file and trees successfully initialized." << std::endl;
    
    // Generator Parameter Tree (gen_param_tree)
    // This TTree stores the parameters used for generation in the file to preserve this data for reference
    int   gen_event_count;
    float gen_pythia_ptbias;
    float gen_jet_pt_min;
    float gen_jet_pt_max;
    float gen_jet_eta_max;
    float gen_pythia_pthatmin;
    float gen_pythia_pthatmax;
    int   gen_thermal_mean;
    int   gen_thermal_sigma;
    float gen_thermal_norm;
    float gen_fastjet_radius;
    float gen_fastjet_pt_min;
    float gen_pythia_beampower;
    float gen_detector_eta;
    float gen_thermal_MH_par_1;
    float gen_thermal_MH_par_2;
    float gen_thermal_MH_par_3;
    float gen_thermal_MH_par_4;
    float gen_thermal_pt_max;
    int   gen_fastjet_accept_index;
    int   gen_ml_soft_jet_index;
    float gen_max_pt_raw;
    
    gen_param_tree->Branch("event_count",           &gen_event_count);
    gen_param_tree->Branch("jet_pt_min",            &gen_jet_pt_min);
    gen_param_tree->Branch("jet_pt_max",            &gen_jet_pt_max);
    gen_param_tree->Branch("jet_eta_max",           &gen_jet_eta_max);
    gen_param_tree->Branch("pythia_ptbias",         &gen_pythia_ptbias);
    gen_param_tree->Branch("pythia_beampower",      &gen_pythia_beampower);
    gen_param_tree->Branch("pythia_pthatmin",       &gen_pythia_pthatmin);
    gen_param_tree->Branch("pythia_pthatmin",       &gen_pythia_pthatmax);
    gen_param_tree->Branch("detector_eta",          &gen_detector_eta);
    gen_param_tree->Branch("fastjet_radius",        &gen_fastjet_radius);
    gen_param_tree->Branch("fastjet_pt_min",        &gen_fastjet_pt_min);
    gen_param_tree->Branch("fastjet_accept_index",  &gen_fastjet_accept_index);
    gen_param_tree->Branch("thermal_mean",          &gen_thermal_mean);
    gen_param_tree->Branch("thermal_sigma",         &gen_thermal_sigma);
    gen_param_tree->Branch("thermal_norm",          &gen_thermal_norm);
    gen_param_tree->Branch("thermal_pt_max",        &gen_thermal_pt_max);
    gen_param_tree->Branch("thermal_MH_par_1",      &gen_thermal_MH_par_1);
    gen_param_tree->Branch("thermal_MH_par_2",      &gen_thermal_MH_par_2);
    gen_param_tree->Branch("thermal_MH_par_3",      &gen_thermal_MH_par_3);
    gen_param_tree->Branch("thermal_MH_par_4",      &gen_thermal_MH_par_4);
    gen_param_tree->Branch("ml_soft_jet_index",     &gen_ml_soft_jet_index);
    gen_param_tree->Branch("jet_pt_raw_max",            &gen_max_pt_raw);
    
    float thermal_norm = 1 / sqrt(2 * math_pi * pow(thermal_sigma,2)); // Sets normalization factor for gaussian
    std::cout << "Thermal Norm: " << thermal_norm << std::endl;
    
    gen_event_count         = event_count;
    gen_pythia_ptbias       = pythia_ptbias;
    gen_jet_pt_min          = jet_pt_min;
    gen_jet_pt_max          = jet_pt_max;
    gen_jet_eta_max         = jet_eta_max;
    gen_pythia_pthatmin     = pythia_pthatmin;
    gen_pythia_pthatmax     = pythia_pthatmax;
    gen_thermal_mean        = thermal_mean;
    gen_thermal_sigma       = thermal_sigma;
    gen_thermal_norm        = thermal_norm;
    gen_fastjet_radius      = fastjet_radius;
    gen_fastjet_pt_min      = fastjet_pt_min;
    gen_pythia_beampower    = pythia_beampower;
    gen_detector_eta        = detector_eta;
    gen_thermal_MH_par_1    = thermal_MH_par_1;
    gen_thermal_MH_par_2    = thermal_MH_par_2;
    gen_thermal_MH_par_3    = thermal_MH_par_3;
    gen_thermal_MH_par_4    = thermal_MH_par_4;
    gen_thermal_pt_max      = thermal_pt_max;
    gen_fastjet_accept_index= fastjet_accept_index;
    gen_ml_soft_jet_index   = ml_soft_jet_index;
    gen_max_pt_raw          = jet_pt_raw_max;
    
    gen_param_tree->Fill();
    gen_param_tree->Write("", TObject::kOverwrite);
    delete gen_param_tree;
    
    std::cout << "Generator_Parameters tree built and written to." << std::endl;
    
    // Maximum number of particles and jets
    const int p_par_max = 1000;
    const int t_par_max = int( thermal_mean + 5 * thermal_sigma );
    const int c_par_max = p_par_max + t_par_max;
    const int p_jet_max = 50;
    const int c_jet_max = 100;
    // Number of constituent particles has to be hard coded (or switch to using vectors)
    
    // PYTHIA Particle Tree (pyth_par_tree)
    // This TTree stores PYTHIA particle information
    // DO NOT USE DOUBLES! Must hard-code numbers for maximum compatibility
    int     p_particle_n;           // particle index
    float   p_particle_pt[p_par_max];     // transverse momentum
    float   p_particle_phi[p_par_max];    // azimuthal angle
    float   p_particle_eta[p_par_max];    // pseudorapidity
    float   p_particle_m[p_par_max];      // particle mass
    int     p_particle_pid[p_par_max];    // numerical code for particle type
    
    pyth_par_tree->Branch("particle_n",     &p_particle_n);
    pyth_par_tree->Branch("particle_pt",    p_particle_pt,    "particle_pt[particle_n]/F");
    pyth_par_tree->Branch("particle_eta",   p_particle_eta,   "particle_eta[particle_n]/F");
    pyth_par_tree->Branch("particle_phi",   p_particle_phi,   "particle_phi[particle_n]/F");
    pyth_par_tree->Branch("particle_m",     p_particle_m,     "particle_m[particle_n]/F");
    pyth_par_tree->Branch("particle_pid",   p_particle_pid,   "particle_pid[particle_n]/I");
    
    std::cout << "PYTHIA_Particle_Tree tree built." << std::endl;
    
    // PYTHIA Jet Tree (pyth_jet_tree)
    // This TTree stores PYTHIA jet information
    int     p_jet_n;
    float   p_jet_pt[p_jet_max];
    float   p_jet_y[p_jet_max];
    float   p_jet_phi[p_jet_max];
    float   p_jet_mass[p_jet_max];
    float   p_jet_area[p_jet_max];
    float   p_jet_area_err[p_jet_max];
    int     p_jet_const_n[p_jet_max];
    float   p_jet_const_pt[p_jet_max][400];
    float   p_jet_const_eta[p_jet_max][400];
    float   p_jet_const_phi[p_jet_max][400];
    
    pyth_jet_tree->Branch("jet_n",          &p_jet_n);
    pyth_jet_tree->Branch("jet_pt",         p_jet_pt,         "jet_pt[jet_n]/F");
    pyth_jet_tree->Branch("jet_y",          p_jet_y,          "jet_y[jet_n]/F");
    pyth_jet_tree->Branch("jet_phi",        p_jet_phi,        "jet_phi[jet_n]/F");
    pyth_jet_tree->Branch("jet_mass",       p_jet_mass,       "jet_mass[jet_n]/F");
    pyth_jet_tree->Branch("jet_area",       p_jet_area,       "jet_area[jet_n]/F");
    pyth_jet_tree->Branch("jet_area_err",   p_jet_area_err,   "jet_area_err[jet_n]/F");
    pyth_jet_tree->Branch("jet_const_n",    p_jet_const_n,    "jet_const_n[jet_n]/I");
    pyth_jet_tree->Branch("jet_const_pt",   p_jet_const_pt,   "jet_const_pt[jet_n][400]/F");
    pyth_jet_tree->Branch("jet_const_eta",  p_jet_const_eta,  "jet_const_eta[jet_n][400]/F");
    pyth_jet_tree->Branch("jet_const_phi",  p_jet_const_phi,  "jet_const_phi[jet_n][400]/F");
    // NOTE: The hard coded size of the second array for the last 3 variables will create a LOT of 0's in the trees - this will not be included in ML
    
    std::cout << "PYTHIA_Jet_Tree tree built." << std::endl;
    
    // Thermal Particle Tree (ther_par_tree)
    // This TTree stores both thermal background particle information
    int     t_particle_n;           // particle index
    float   t_particle_pt[t_par_max];     // transverse momentum
    float   t_particle_phi[t_par_max];    // azimuthal angle
    float   t_particle_eta[t_par_max];    // pseudorapidity
    float   t_particle_m[t_par_max];      // particle mass
    
    if ( export_thermal ) {
        ther_par_tree = new TTree("Thermal_Particle_Tree", "Tree of particles from thermal toy model");
        ther_par_tree->Branch("particle_n",     &t_particle_n);
        ther_par_tree->Branch("particle_pt",    t_particle_pt,    "particle_pt[particle_n]/F");
        ther_par_tree->Branch("particle_eta",   t_particle_eta,   "particle_eta[particle_n]/F");
        ther_par_tree->Branch("particle_phi",   t_particle_phi,   "particle_phi[particle_n]/F");
        ther_par_tree->Branch("particle_m",     t_particle_m,     "particle_m[particle_n]/F");
        std::cout << "Thermal_Jet_Tree tree built." << std::endl;
    }
    
    // Combined Particle Tree (comb_par_tree)
    // This TTree stores both PYTHIA and thermal background particle information
    int   c_particle_n;             // number of total particles in event
    float c_particle_pt[c_par_max];      // transverse momentum
    float c_particle_phi[c_par_max];     // azimuthal angle
    float c_particle_eta[c_par_max];     // pseudorapidity
    float c_particle_m[c_par_max];       // particle mass
    int   c_particle_source[c_par_max];  // source of particles (0 = PYTHIA, 1 = Thermal)
    
    comb_par_tree->Branch("particle_n",     &c_particle_n);
    comb_par_tree->Branch("particle_pt",    c_particle_pt,     "particle_pt[particle_n]/F");
    comb_par_tree->Branch("particle_eta",   c_particle_eta,    "particle_eta[particle_n]/F");
    comb_par_tree->Branch("particle_phi",   c_particle_phi,    "particle_phi[particle_n]/F");
    comb_par_tree->Branch("particle_m",     c_particle_m,      "particle_m[particle_n]/F");
    comb_par_tree->Branch("particle_source",c_particle_source, "particle_source[particle_n]/I");
    
    std::cout << "Combined_Particle_Tree tree built." << std::endl;
    
    // Combined Jet Tree (comb_jet_tree)
    // This TTree stores combined event jet information
    int     c_jet_n;
    float   c_jet_pt[c_jet_max];
    float   c_jet_y[c_jet_max];
    float   c_jet_phi[c_jet_max];
    float   c_jet_mass[c_jet_max];
    float   c_jet_area[c_jet_max];
    float   c_jet_area_err[c_jet_max];
    int     c_jet_const_n[c_jet_max];
    float   c_jet_const_pt[c_jet_max][400];
    float   c_jet_const_eta[c_jet_max][400];
    float   c_jet_const_phi[c_jet_max][400];
    
    comb_jet_tree->Branch("jet_n",          &c_jet_n);
    comb_jet_tree->Branch("jet_pt",         c_jet_pt,         "jet_pt[jet_n]/F");
    comb_jet_tree->Branch("jet_y",          c_jet_y,          "jet_y[jet_n]/F");
    comb_jet_tree->Branch("jet_phi",        c_jet_phi,        "jet_phi[jet_n]/F");
    comb_jet_tree->Branch("jet_mass",       c_jet_mass,       "jet_mass[jet_n]/F");
    comb_jet_tree->Branch("jet_area",       c_jet_area,       "jet_area[jet_n]/F");
    comb_jet_tree->Branch("jet_area_err",   c_jet_area_err,   "jet_area_err[jet_n]/F");
    comb_jet_tree->Branch("jet_const_n",    c_jet_const_n,    "jet_const_n[jet_n]/I");
    comb_jet_tree->Branch("jet_const_pt",   c_jet_const_pt,   "jet_const_pt[jet_n][400]/F");
    comb_jet_tree->Branch("jet_const_eta",  c_jet_const_eta,  "jet_const_eta[jet_n][400]/F");
    comb_jet_tree->Branch("jet_const_phi",  c_jet_const_phi,  "jet_const_phi[jet_n][400]/F");
    // NOTE: The hard coded size of the second array for the last 3 variables will create a LOT of 0's in the trees - this will not be included in ML
    
    std::cout << "Combined_Jet_Tree tree built." << std::endl;
    
    // PYTHIA Setup
    // LHC process and output selection. Initialization.
    Pythia pythia;
    Settings& pythia_settings = pythia.settings;
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");
    pythia_settings.parm("Beams:eCM", pythia_beampower);  // [GeV]
    pythia_settings.parm("PhaseSpace:pTHatMin", pythia_pthatmin);
    pythia_settings.parm("PhaseSpace:pTHatMax", pythia_pthatmax);
    pythia.readString("HardQCD:all = on");  // Turns on hard scattering
    if ( pythia_ptbias != 0. ) {
        pythia_settings.parm("PhaseSpace:bias2SelectionPow", pythia_ptbias);
        pythia.readString("PhaseSpace:bias2Selection = on");
        pythia.readString("PhaseSpace:bias2SelectionRef = 100.");
    }
    pythia.init();
    
    std::cout << "PYTHIA settings initialized." << std::endl;
    
    // Thermal Background Setup
    // This initializes the distribution functions for thermal particles
    TF1* thermal_n_func   = new TF1("thermal_n_func",   Func_Gaussian, 0, t_par_max, 3);
    TF1* thermal_pt_func  = new TF1("thermal_pt_func",  Func_ModifiedHagedorn, 0, thermal_pt_max, 4);
    TF1* thermal_eta_func = new TF1("thermal_eta_func", Func_Flat, -detector_eta, detector_eta, 1);
    TF1* thermal_phi_func = new TF1("thermal_phi_func", Func_Flat, -math_pi, math_pi, 1);
    thermal_n_func   ->SetParameters(thermal_norm, thermal_mean, thermal_sigma);
    thermal_pt_func  ->SetParameters(thermal_MH_par_1, thermal_MH_par_2, thermal_MH_par_3, thermal_MH_par_4);
    thermal_eta_func ->SetParameter(0, 1.);
    thermal_phi_func ->SetParameter(0, 1.);
    
    std::cout << "Thermal model settings initialized." << std::endl;
    
    // FastJet Setup
    // These settings are used to cluster both PYTHIA events and combined events
    fastjet::JetDefinition jet_definition(fastjet::antikt_algorithm, fastjet_radius);
    fastjet::AreaDefinition area_definition;
    if (!d_use_voronoi) {
        double ghost_eta_max = detector_eta;
        double ghost_area    = 0.05;
        int    active_area_repeats = 5;
        fastjet::GhostedAreaSpec ghost_spec(ghost_eta_max, active_area_repeats, ghost_area);
        area_definition = fastjet::AreaDefinition(fastjet::active_area, ghost_spec);
    }
    else {
        double effective_Rfact = 1.0;
        area_definition = fastjet::VoronoiAreaSpec(effective_Rfact);
    }
    
    std::cout << "FastJet settings initialized." << std::endl;
    
//    float const_1_pt_arr[10000];
//    float const_2_pt_arr[10000];
//    float const_3_pt_arr[10000];
//    int   const_arr_count = 0;
    
    // Generate Events
    // Loop to generate PYTHIA events. Skip if error.
    int const_arr_size = sizeof(float) * 400;
    int e = 0;
    while ( e < event_count ) {
        // ----- PYTHIA -----
        if (!pythia.next()) continue;

        p_particle_n = 0;

        for ( int p = 0; p < pythia.event.size(); p++ ) {
            if (!pythia.event[p].isFinal()) continue;                       // Ends if final particle
            if (!pythia.event[p].isCharged()) continue;                     // Skips neutral particles
            if (pythia.event[p].pT() < 0.15) continue;                      // Skips particles with pT < 0.15
            if (fabs(pythia.event[p].eta()) > detector_eta) continue;       // Skips particles with rapidity outside the detector rapidity
            p_particle_pt[p_particle_n]  = pythia.event[p].pT();
            p_particle_eta[p_particle_n] = pythia.event[p].eta();
            p_particle_phi[p_particle_n] = pythia.event[p].phi();
            p_particle_pid[p_particle_n] = pythia.event[p].id();
            p_particle_m[p_particle_n]   = pythia.event[p].m();             // Outputs mass as negative??
            p_particle_n++;
        }
        
        if ( p_particle_n == 0 ) continue;
                
        // ---------- FASTJET ----------
        TLorentzVector pyth_vector;
        std::vector<fastjet::PseudoJet> pyth_particles;
        for ( int p = 0 ; p < p_particle_n ; p++) {
            pyth_vector.SetPtEtaPhiM( p_particle_pt[p], p_particle_eta[p], p_particle_phi[p], p_particle_m[p] );
            pyth_particles.push_back( fastjet::PseudoJet( pyth_vector.Px(), pyth_vector.Py(), pyth_vector.Pz(), pyth_vector.E() ) );
        }
        fastjet::ClusterSequenceArea p_jet_clusters( pyth_particles, jet_definition, area_definition );
        std::vector<fastjet::PseudoJet> p_jet_list = sorted_by_pt( p_jet_clusters.inclusive_jets(fastjet_pt_min) );
        
//        std::cout << "Event: " << e << " -----" << std::endl;
//        for (int z = 0 ; z < p_jet_list.size() ; z++ ) {
//            std::cout << "Jet: " << z << ", pT: " << p_jet_list[z].pt() << std::endl;
//        }
        
        // Checks jets to see if the event should be accepted
        p_jet_n = p_jet_list.size();
        if ( p_jet_n == 0 ) continue;
        if ( p_jet_n >= 100 ) p_jet_n = 100;    // Prevents more jets from being recorded than space in jet_n vector - this should never happen
        int max_accept_index = p_jet_n;
        if ( fastjet_accept_index != 0 ) max_accept_index = fastjet_accept_index;
        bool accept_event = false;
        float jet_pt_temp;
        float jet_y_temp;
        for ( int j = 0 ; j < max_accept_index ; j++ ) {
            jet_pt_temp = p_jet_list[j].pt();
            jet_y_temp = p_jet_list[j].rap();
            if ( (jet_pt_temp >= jet_pt_min) && (jet_pt_temp <= jet_pt_max) &&
                (jet_y_temp >= -jet_eta_max) && (jet_y_temp <= jet_eta_max) ) {
                accept_event = true;
                break;
            }
            else continue;
        }
        
        if ( accept_event ) {
            // Add PYTHIA Particles
            pyth_par_tree->Fill();
            
//            // Used to check particle pT
//            if ( const_arr_count < 10000 ) {
//                float p_pt_arr[200];
//                for ( int i = 0 ; i < p_particle_n ; i++ ) p_pt_arr[i] = p_particle_pt[i];
//                std::sort(std::begin(p_pt_arr), std::begin(p_pt_arr) + p_particle_n, std::greater<>());
//                const_1_pt_arr[const_arr_count] = p_pt_arr[0];
//                const_2_pt_arr[const_arr_count] = p_pt_arr[1];
//                const_3_pt_arr[const_arr_count] = p_pt_arr[2];
////                std::cout << p_pt_arr[0] << ", " << p_pt_arr[1] << ", " << p_pt_arr[2] << std::endl;
//                const_arr_count++;
//            }
            
            // Add PYTHIA Jets
            for ( int j = 0 ; j < p_jet_n ; j++ ) {
                // Add jet data
                p_jet_pt[j]         = p_jet_list[j].pt();
                p_jet_y[j]          = p_jet_list[j].rap();
                p_jet_phi[j]        = p_jet_list[j].phi();
                p_jet_mass[j]       = p_jet_list[j].m();
                p_jet_area[j]       = p_jet_list[j].area();
                p_jet_area_err[j]   = p_jet_list[j].area_error();
                std::vector<fastjet::PseudoJet> p_jet_constituents = sorted_by_pt(p_jet_list[j].constituents());
                p_jet_const_n[j]    = p_jet_constituents.size();
                // Zero out constituent array and add constituent data
                memset(p_jet_const_pt[j], 0, const_arr_size);
                memset(p_jet_const_eta[j], 0, const_arr_size);
                memset(p_jet_const_phi[j], 0, const_arr_size);
                for ( int p = 0 ; p < p_jet_const_n[j] ; p++) {
                    p_jet_const_pt[j][p]    = p_jet_constituents[p].pt();
                    p_jet_const_eta[j][p]   = p_jet_constituents[p].eta();
                    p_jet_const_phi[j][p]   = p_jet_constituents[p].phi();
                }
            }
            pyth_jet_tree->Fill();
            // Add Combined Event Particles
            c_particle_n = 0;
            for ( int p = 0; p < pythia.event.size(); p++ ) {
                if (!pythia.event[p].isFinal()) continue;                   // Ends if final particle
                if (!pythia.event[p].isCharged()) continue;                 // Skips neutral particles
                if (pythia.event[p].pT() < 0.15) continue;                  // Skips particles with pT < 0.15
                if (fabs(pythia.event[p].eta()) > detector_eta) continue;   // Skips particles with rapidity outside the detector rapidity
                c_particle_pt[c_particle_n]     = pythia.event[p].pT();
                c_particle_eta[c_particle_n]    = pythia.event[p].eta();
                c_particle_phi[c_particle_n]    = pythia.event[p].phi();
                c_particle_m[c_particle_n]      = pythia.event[p].m();       // Outputs mass as negative
                c_particle_source[c_particle_n] = 1;
                c_particle_n++;
            }
            int thermal_n = Generate_N(thermal_n_func, thermal_mean);
            if ( export_thermal ) t_particle_n = thermal_n;
            for ( int p = 0 ; p < thermal_n ; p++) {
                c_particle_pt[c_particle_n]     = Generate_Pt(thermal_pt_func, thermal_pt_max);
                c_particle_eta[c_particle_n]    = Generate_Eta(thermal_eta_func, detector_eta);
                c_particle_phi[c_particle_n]    = Generate_Phi(thermal_phi_func);
                c_particle_m[c_particle_n]      = m_pion;
                c_particle_source[c_particle_n] = 0;
                if ( export_thermal ) {
                    t_particle_pt[p]    = c_particle_pt[c_particle_n];
                    t_particle_eta[p]   = c_particle_eta[c_particle_n];
                    t_particle_phi[p]   = c_particle_phi[c_particle_n];
                    t_particle_m[p]     = c_particle_m[c_particle_n];
                }
                c_particle_n++;
            }
            comb_par_tree->Fill();
            if ( export_thermal ) ther_par_tree->Fill();;
            
            // Cluster Combined Events
            TLorentzVector comb_vector;
            std::vector<fastjet::PseudoJet> comb_particles;
            for ( int p = 0 ; p < c_particle_n ; p++) {
                comb_vector.SetPtEtaPhiM( c_particle_pt[p], c_particle_eta[p], c_particle_phi[p], c_particle_m[p] );
                comb_particles.push_back( fastjet::PseudoJet( comb_vector.Px(), comb_vector.Py(), comb_vector.Pz(), comb_vector.E() ) );
            }
            fastjet::ClusterSequenceArea c_jet_clusters(comb_particles, jet_definition, area_definition);
            std::vector<fastjet::PseudoJet> c_jet_list = sorted_by_pt(c_jet_clusters.inclusive_jets(fastjet_pt_min));
            c_jet_n = c_jet_list.size();
            if ( c_jet_n >= 100 ) c_jet_n = 100; // Limits to top 100 jets (this should never happen)
            for ( int j = 0 ; j < c_jet_n ; j++ ) {
                // Add jet data
                c_jet_pt[j]         = c_jet_list[j].pt();
                c_jet_y[j]          = c_jet_list[j].rap();
                c_jet_phi[j]        = c_jet_list[j].phi();
                c_jet_mass[j]       = c_jet_list[j].m();
                c_jet_area[j]       = c_jet_list[j].area();
                c_jet_area_err[j]   = c_jet_list[j].area_error();
                std::vector<fastjet::PseudoJet> c_jet_constituents = sorted_by_pt(c_jet_list[j].constituents());
//                std::cout << "Event: " << e << " Jet: " << j << " with " << c_jet_constituents.size() << " constituents -----" << std::endl;
//                for (int z = 0 ; z < c_jet_constituents.size() ; z++ ) {
//                    std::cout << "Const: " << z << ", pT: " << c_jet_constituents[z].pt() << std::endl;
//                }
                c_jet_const_n[j]    = c_jet_constituents.size();
                // Zero out constituent array and add constituent data
                memset(c_jet_const_pt[j], 0, const_arr_size);
                memset(c_jet_const_eta[j], 0, const_arr_size);
                memset(c_jet_const_phi[j], 0, const_arr_size);
                for ( int p = 0 ; p < c_jet_const_n[j] ; p++) {
                    c_jet_const_pt[j][p]   = c_jet_constituents[p].pt();
                    c_jet_const_eta[j][p]  = c_jet_constituents[p].eta();
                    c_jet_const_phi[j][p]  = c_jet_constituents[p].phi();
                }
            }
            comb_jet_tree->Fill();
            
            // Prints values for jets from events that have been accepted
            if ( ( e % print_every_x ) == 0 ) std::cout << "Event " << e << " has " << c_jet_n << " total jets with " << p_jet_n << " truth jets." << std::endl;
            if ( print_out && ( e % print_every_x ) == 0 ) {
                std::cout << "=== NEW EVENT # " << e << " ===" << std::endl;
                std::cout << "--- Summary of PYTHIA Jets --- " << std::endl;
                for ( int j = 0; j < p_jet_n; j++ ) {
                    std::cout << " -> jet " << j << ", pT / y / phi = " << p_jet_pt[j] << " / " << p_jet_y[j] << " / " << p_jet_phi[j] << std::endl;
                }
                std::cout << "--- Summary of Combined Jets --- " << std::endl;
                for ( int j = 0; j < c_jet_n; j++ ) {
                    std::cout << " -> Jet " << j << ", pT / y / phi = " << c_jet_pt[j] << " / " << c_jet_y[j] << " / " << c_jet_phi[j] << std::endl;
                }
            }
            e++;
        }
        else continue;
    }
    std::cout << "Event generation complete!" << std::endl;
    if ( debug ) pythia.stat(); // Prints PYTHIA stats at the end
    
//    std::cout << "\n SUMMARY OF CONSTITUENTS" << std::endl;
//    std::cout << "Const. 1 Mean: " << TMath::Mean(10000, const_1_pt_arr) << std::endl;
//    std::cout << "Const. 2 Mean: " << TMath::Mean(10000, const_2_pt_arr) << std::endl;
//    std::cout << "Const. 3 Mean: " << TMath::Mean(10000, const_3_pt_arr) << "\n" << std::endl;
    
    // Write trees to output file
    pyth_par_tree ->Write("", TObject::kOverwrite);
    pyth_jet_tree ->Write("", TObject::kOverwrite);
    comb_par_tree ->Write("", TObject::kOverwrite);
    comb_jet_tree ->Write("", TObject::kOverwrite);
    if ( export_thermal ) ther_par_tree ->Write("", TObject::kOverwrite);
    
    std::cout << "Output trees written to." << std::endl;
    
    // Delete particle trees (not needed for ML Prep)
    delete pyth_par_tree;
    delete comb_par_tree;
    
    // Start ML Prep
    std::cout << "Starting Machine Learning Prep..." << std::endl;
    
    char mlprep_file_path[500];
    snprintf(mlprep_file_path, 500, "%s/Full_%s_ML_Prep.root", output_directory, output_file_params);
    TFile* mlprep_file = new TFile(mlprep_file_path, "UPDATE");
    mlprep_file->cd();
    char mlprep_tree_name[200];
    snprintf(mlprep_tree_name, 200, "ML_%s", output_file_params);
    char mlprep_tree_description[200];
    snprintf(mlprep_tree_description, 200, "Generated with bias of pT^%.0f, %.1f - %.1f GeV, %i events, jet match radius of %.2f", pythia_ptbias, jet_pt_min, jet_pt_max, gen_event_count, fastjet_match_radius);
    TTree* mlprep_tree = new TTree(mlprep_tree_name, mlprep_tree_description);
    
    // Jet Machine Learning Tree (Jet_ML_...)
    // This Flat TTree stores jets for machine learning
    int    jet_index;
    float  jet_pt_raw;
    float  jet_pt_corr;
    float  jet_pt_true;
    float  jet_mass;
    float  jet_area;
    float  jet_area_err;
    int    jet_const_n;
    float  const_pt_mean;
    float  const_pt_median;
    float  const_1_pt;
    float  const_2_pt;
    float  const_3_pt;
    float  const_4_pt;
    float  const_5_pt;
    float  const_6_pt;
    float  const_7_pt;
    float  const_8_pt;
    float  const_9_pt;
    float  const_10_pt;
    float  jet_y;
    float  jet_phi;
    float  jet_rho;
    
    mlprep_tree->Branch("jet_index",            &jet_index);
    mlprep_tree->Branch("jet_pt_raw",           &jet_pt_raw,        "jet_pt_raw/F");
    mlprep_tree->Branch("jet_pt_corr",          &jet_pt_corr,       "jet_pt_corr/F");
    mlprep_tree->Branch("jet_pt_true",          &jet_pt_true,       "jet_pt_true/F");
    mlprep_tree->Branch("jet_mass",             &jet_mass,          "jet_mass/F");
    mlprep_tree->Branch("jet_area",             &jet_area,          "jet_area/F");
    mlprep_tree->Branch("jet_area_err",         &jet_area_err,      "jet_area_err/F");
    mlprep_tree->Branch("jet_const_n",          &jet_const_n,       "jet_const_n/I");
    mlprep_tree->Branch("const_pt_mean",        &const_pt_mean,     "const_pt_mean/F");
    mlprep_tree->Branch("const_pt_median",      &const_pt_median,   "const_pt_median/F");
    mlprep_tree->Branch("const_1_pt",           &const_1_pt,        "const_1_pt/F");
    mlprep_tree->Branch("const_2_pt",           &const_2_pt,        "const_2_pt/F");
    mlprep_tree->Branch("const_3_pt",           &const_3_pt,        "const_3_pt/F");
    mlprep_tree->Branch("const_4_pt",           &const_4_pt,        "const_4_pt/F");
    mlprep_tree->Branch("const_5_pt",           &const_5_pt,        "const_5_pt/F");
    mlprep_tree->Branch("const_6_pt",           &const_6_pt,        "const_6_pt/F");
    mlprep_tree->Branch("const_7_pt",           &const_7_pt,        "const_7_pt/F");
    mlprep_tree->Branch("const_8_pt",           &const_8_pt,        "const_8_pt/F");
    mlprep_tree->Branch("const_9_pt",           &const_9_pt,        "const_9_pt/F");
    mlprep_tree->Branch("const_10_pt",          &const_10_pt,       "const_10_pt/F");
    mlprep_tree->Branch("jet_y",                &jet_y,             "jet_y/F");
    mlprep_tree->Branch("jet_phi",              &jet_phi,           "jet_phi/F");
    mlprep_tree->Branch("jet_rho",              &jet_rho,           "jet_rho/F");
    
    std::cout << "Jet_ML_ tree built." << std::endl;
    
    // Iterate through events to match and extract jets
    int   truth_matched_counter = 0;
    int   truth_unmatched_counter = 0;
    int   truth_outside_pt_counter = 0;
    int   truth_outside_y_counter = 0;
    int   raw_outside_pt_counter = 0;
    float fastjet_match_rad_sq = pow(fastjet_match_radius, 2);
    float background_pt_median[event_count];
    std::cout << "Starting jet matching..." << std::endl;
    for ( e = 0 ; e < event_count ; e++ ) {
        if ( print_out && ( e % print_every_x ) == 0 ) std::cout << "Getting jets from event " << e << std::endl;
        comb_jet_tree->GetEntry(e);
        pyth_jet_tree->GetEntry(e);
        
        // Finds the median pt of the background jets
        float  background_jet_pt[100];
        float  background_jet_area[100];
        float  background_rho_arr[100];
        float  background_area_median;
        int    background_jet_n = 0;
        for ( int j = 2 ; j < c_jet_n ; j++ ) { // Skips top 2 hardest jets in each event
            background_jet_pt[j-2]   = c_jet_pt[j];
            background_jet_area[j-2] = c_jet_area[j];
            background_rho_arr[j-2]  = c_jet_pt[j] / c_jet_area[j];
            background_jet_n++;
        }
        background_area_median  = TMath::Median(background_jet_n, background_jet_area);
        background_pt_median[e] = TMath::Median(background_jet_n, background_jet_pt) * background_area_median;
        float background_rho    = TMath::Median(background_jet_n, background_rho_arr);
        char rho_output[100];
        snprintf(rho_output, 100, "Event %i / %i bg_jets / rho: %.3f", e, background_jet_n, background_rho);
        if ( debug ) std::cout << rho_output << std::endl;
        
        // Iterate through pythia jets to match jet pT_true
        int softest_jet_index = ml_soft_jet_index;
        if ( ml_soft_jet_index == 0 ) softest_jet_index = c_jet_n;
        for ( int k = 0 ; k < p_jet_n ; k++ ) {
            int c_match = -1;
            int c_match_tracker[100] = {0};
            
            // Iterates through combined jets from hardest to softest to try matching
            for ( int j = 0 ; j < softest_jet_index ; j++ ) {
                if ( c_match_tracker[j] == 1 ) continue; // Skips if the combined jet has already been paired
                float jet_y_sq = pow((p_jet_y[k] - c_jet_y[j]), 2);
                float jet_phi_sq = pow((p_jet_phi[k] - c_jet_phi[j]), 2);
                if ( jet_phi_sq + jet_y_sq > fastjet_match_rad_sq ) { // Triggers the 2pi offset check only if needed (hopefully speeds this up)
                    float jet_phi_sq_B;
                    if ( p_jet_phi[k] > c_jet_phi[j] ) {
                        jet_phi_sq_B = pow((p_jet_phi[k] - 2*math_pi - c_jet_phi[j]), 2);
                    }
                    else {
                        jet_phi_sq_B = pow((c_jet_phi[k] - 2*math_pi - p_jet_phi[j]), 2);
                    }
                    if ( jet_phi_sq > jet_phi_sq_B ) jet_phi_sq = jet_phi_sq_B;
                    
                    if ( jet_phi_sq + jet_y_sq > fastjet_match_rad_sq ) continue; // Skips if the distance between PYTHIA and combined jets is larger than match radius
                }
                c_match = j;
                c_match_tracker[j] = 1;
                if ( print_out && ( e % print_every_x ) == 0 ) std::cout << "Combined jet " << j << " matches truth jet " << k << std::endl;
                if ( debug ) std::cout << "Match found!" << std::endl;
                break;
            }
            
            // Skips if jet can't be matched to a PYTHIA jet OR if jet raw is greater than 200 GeV
            if ( c_match < 0 ) {
                truth_unmatched_counter++;
                if ( debug) std::cout << "Event: " << e << " Truth Jet: " << k << " has no match found." << std::endl;
                if ( debug) std::cout << "    p_phi: " << p_jet_phi[k] << ", p_y: " << p_jet_y[k] << std::endl;
                continue;
            }
            // Skips if jet pT_raw is greater than the pax accepted pT_raw
            if ( c_jet_pt[c_match] > jet_pt_raw_max ) {
                raw_outside_pt_counter++;
                continue;
            }
            // Skips if jet pT is outside range
            if ( p_jet_pt[k] < jet_pt_min || p_jet_pt[k] > jet_pt_max ) {
                truth_outside_pt_counter++;
                if ( debug) std::cout << "Event: " << e << " Truth Jet: " << k << " is outside pt min/max. pT = " << p_jet_pt[k] << std::endl;
                continue;
            }
            // Skips if jet rapidity is outside range
            if ( p_jet_y[k] < -jet_eta_max || p_jet_y[k] > jet_eta_max ) {
                truth_outside_y_counter++;
                if ( debug) std::cout << "Event: " << e << " Truth Jet: " << k << " is outside rapidity. y = " << p_jet_y[k] << std::endl;
                continue;
            }
            truth_matched_counter++;
            
            // Iterate through constituent particles to collect their pT for mean and median
            float  jet_const_pt_arr[400];
            for ( int p = 0 ; p < c_jet_const_n[c_match] ; p++ ) {
                jet_const_pt_arr[p] = c_jet_const_pt[c_match][p];
            }
            std::sort(jet_const_pt_arr, jet_const_pt_arr + jet_const_n, std::greater<>());
            
            jet_index       = c_match;
            jet_pt_raw      = c_jet_pt[c_match];
            jet_pt_corr     = c_jet_pt[c_match] - (background_rho * c_jet_area[c_match]);
            jet_mass        = c_jet_mass[c_match];
            jet_area        = c_jet_area[c_match];
            jet_area_err    = c_jet_area_err[c_match];
            jet_const_n     = c_jet_const_n[c_match];
            const_pt_mean   = TMath::Mean(jet_const_n, jet_const_pt_arr);
            const_pt_median = TMath::Median(jet_const_n, jet_const_pt_arr);
            const_1_pt      = jet_const_pt_arr[0];
            const_2_pt      = jet_const_pt_arr[1];
            const_3_pt      = jet_const_pt_arr[2];
            const_4_pt      = jet_const_pt_arr[3];
            const_5_pt      = jet_const_pt_arr[4];
            const_6_pt      = jet_const_pt_arr[5];
            const_7_pt      = jet_const_pt_arr[6];
            const_8_pt      = jet_const_pt_arr[7];
            const_9_pt      = jet_const_pt_arr[8];
            const_10_pt     = jet_const_pt_arr[9];
            jet_y           = c_jet_y[c_match];
            jet_phi         = c_jet_phi[c_match];
            jet_rho         = background_rho;
            jet_pt_true     = p_jet_pt[k];
            
            mlprep_tree->Fill();
            
            if ( print_out && jet_pt_true != 0 && (e % print_every_x) == 0 ) std::cout << e << "-" << c_match << ": Truth_PYTHIA Jet: " << jet_pt_true << " ----- " << std::endl;
        }
    }
    
    char truth_match_results[200];
    snprintf(truth_match_results, 200, "Matched %i jets. %i jets unmatched. %i jets outside pt. %i jets outside y. %i jets with pt_raw above %.0f GeV. Unmatched ratio is: %.4f", truth_matched_counter, truth_unmatched_counter, truth_outside_pt_counter, truth_outside_y_counter, raw_outside_pt_counter, jet_pt_raw_max, float(truth_unmatched_counter)/(truth_matched_counter + truth_unmatched_counter) );
    std::cout << truth_match_results << std::endl;
    
    mlprep_tree->Write("", TObject::kOverwrite);
    
    std::cout << "ML tree written to." << std::endl;
    
    delete mlprep_tree;
    delete pyth_jet_tree;
    delete comb_jet_tree;
    output_file->Close();
    mlprep_file->Close();
    
    std::cout << "File saved and closed. Function complete!" << std::endl;
}

} // Ends the Jet_Generator namespace



// ----- DEMO CODE -----



/*
   Demo()
   Note, this is included for example purposes.
 */
int Demo() {

    // Defines working directories for data and plots
    char dir_master[200];
    snprintf(dir_master, 200, "../../Files/Demo");
    
    char dir_data[200];
    snprintf(dir_data, 200, "%s/Data", dir_master);
    
    char dir_plots[200];
    snprintf(dir_plots, 200, "%s/Plots", dir_master);
    
    // Creates directories if they do not already exist
    std::__fs::filesystem::create_directories(dir_master);
    std::__fs::filesystem::create_directories(dir_data);
    std::__fs::filesystem::create_directories(dir_plots);
    
    // Example using MINIMUM number of parameters and default values . Note beam_power and detector_eta have defai
    Jet_Generator::Jet_Generator_Optimized(
        "Train",        // output_base_name (pt_bias, pt_min, pt_max, and event_cout will be added to this name in the functions)
        dir_data,       // output_directory
        1000,           // event_count (Generates 100,000 events but will produce MORE total jets)
        0.,             // pythia_ptbias (Here NO BIAS is applied)
        10.,            // jet_pt_min (Events will have at least one jet above
        90.             // jet_pt_max
        );
    
    // Example defining ALL parameters.
    // NOTE: Must define ALL parameters up to the last one used!
    // (e.g. If defining fastjet_radius then jet_eta_max, pt_hat_min, pt_hat_max, must also be defined - it's not like Python!)
    Jet_Generator::Jet_Generator_Optimized(
        "Train",    // output_base_name (Typically just "Train" or "Test")
        dir_data,   // output_directory
        100000,     // event_count [int]
        8.,         // pythia_ptbias (Generates events with bias of pT^N for event energies distribution)
        10.,        // jet_pt_min [GeV]
        90.,        // jet_pt_max [GeV]
        0.5,        // jet_eta_max (Reasonable value is detector_eta - fastjet_radius)
        7.5,        // pythia_pthatmin (Reasonable value is 0.75 * jet_pt_min)
        0.,         // pythia_pthatmax (If =0. then no upper limit)
        0.4,        // fastjet_radius
        8.,         // fastjet_pt_min (Minimum pT considered by FastJet for a jet)
        0.3,        // fastjet_match_radius (Reasonable value is 0.75 * fastjet_radius)
        0,          // fastjet_accept_index [int] (If =0, accepts events with any jet in min/max. If >0, only accepts events with one of top N jets in min/max.)
        0,          // ml_soft_jet_index [int] (If =0, tries to match with all PYTHIA jets. If >0, matches only to top N PYTHIA jets)
        1800,       // thermal_mean
        200,        // thermal_sigma
        2760.,      // pythia_beampower [GeV] (Simulated beam power in the detector)
        0.9,        // detector_eta (Maximum rapidity of the detector)
        64547,      // Modified Hagedorn function parameter 1
        3.076,      // Modified Hagedorn function parameter 2
        1.126,      // Modified Hagedorn function parameter 3
        -8.491,     // Modified Hagedorn function parameter 4
        100         // Max thermal pT (typically twice the thermal mean)
        );
    
    return 0;
}
