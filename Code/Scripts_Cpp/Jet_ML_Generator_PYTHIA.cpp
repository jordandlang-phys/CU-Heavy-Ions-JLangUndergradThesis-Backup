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
#include "TRandom.h"
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



// ----- PYTHIA GENERATOR -----



/*
    Pythia_Generator()
    Generates events in PYTHIA, clusters into jets, and filters based on input parameters.
    Returns a character string pointer of the output file path.
 */
char* Pythia_Generator(
    char  output_base_name[100],            // Base/stem of output file name (no file extensions or prefixes).
    char  output_directory[400],            // Path to desired output directory.
    int   event_count,                      // Number of collision events to generate. The number of jets will be greater than this.
    float pt_bias_power,                    // Bias applied to the particle pT distribution. If =0. then no bias applied. If =-1 then distribution is flattened
    float jet_pt_min,                       // Minimum jet pT to accept. At least 1 jet per event will be greater than this.
    float jet_pt_max,                       // Maximum jet pT to accept. At least 1 jet per event will be less than this.
    float jet_eta_max = d_jet_eta_max,      // Maximum rapidity to consider.
    float pt_hat_min = d_pt_hat_min,        // [GeV] ptHatMin value.
    float pt_hat_max = d_pt_hat_max,        // [GeV] ptHatMax value. If 0. then no upper limit.
    int   accept_jet_index = d_accept_jet_index,    // If =0, accepts events with any jet between min/max. If >0, only accepts events with one of top N jets between min/max.
    float fastjet_radius = d_fastjet_radius,        // Jet radius used by FastJet.
    float fastjet_pt_min = d_fastjet_pt_min,        // [GeV] Minimum pT considered by FastJet for a
    float beam_power = d_beam_power,        // [GeV] Simulated beam power in the dector. Default defined by jet_ml_constants.h file.
    float detector_eta = d_detector_eta     // Maximum rapidity of the detector. Default defined by jet_ml_constants.h file.
    ) {
    
    char output_file_path[500];
    snprintf(output_file_path, 500, "%s/Pyth_%s_B%.0f_%.0f_%.0f_N%i.root", output_directory, output_base_name, pt_bias_power, jet_pt_min, jet_pt_max, event_count);
    TFile* output_file      = new TFile(output_file_path, "UPDATE");
    TTree* pythia_tree      = new TTree("Pythia_Tree","Tree of particle jet events from PYTHIA p+p collisions");
    TTree* jet_tree         = new TTree("FastJet_Tree","Tree of jet clusters by event from PYTHIA");
    
    if (print_out) std::cout << "Output file and trees successfully generated." << std::endl;
    
    // Generator Parameter Tree
    // This TTree stores the parameters used for generation in the file to preserve this data for reference
    TTree* gen_param_tree   = new TTree("Generator_Parameters", "Flat tree of parameters used for file generation");
    
    int   gen_event_count;
    float gen_pt_bias_power;
    float gen_jet_pt_min;
    float gen_jet_pt_max;
    float gen_jet_eta_max;
    float gen_pt_hat_min;
    float gen_pt_hat_max;
    float gen_fastjet_radius;
    float gen_fastjet_pt_min;
    float gen_beam_power;
    float gen_detector_eta;
    
    gen_param_tree->Branch("event_count",       &gen_event_count);
    gen_param_tree->Branch("pt_bias_power",     &gen_pt_bias_power);
    gen_param_tree->Branch("jet_pt_min",        &gen_jet_pt_min);
    gen_param_tree->Branch("jet_pt_max",        &gen_jet_pt_max);
    gen_param_tree->Branch("jet_eta_max",       &gen_jet_eta_max);
    gen_param_tree->Branch("pt_hat_min",        &gen_pt_hat_min);
    gen_param_tree->Branch("pt_hat_max",        &gen_pt_hat_max);
    gen_param_tree->Branch("fastjet_radius",    &gen_fastjet_radius);
    gen_param_tree->Branch("fastjet_pt_min",    &gen_fastjet_pt_min);
    gen_param_tree->Branch("beam_power",        &gen_beam_power);
    gen_param_tree->Branch("detector_eta",      &gen_detector_eta);
    
    gen_event_count     = event_count;
    gen_pt_bias_power   = pt_bias_power;
    gen_jet_pt_min      = jet_pt_min;
    gen_jet_pt_max      = jet_pt_max;
    gen_jet_eta_max     = jet_eta_max;
    gen_pt_hat_min      = pt_hat_min;
    gen_pt_hat_max      = pt_hat_max;
    gen_fastjet_radius  = fastjet_radius;
    gen_fastjet_pt_min  = fastjet_pt_min;
    gen_beam_power      = beam_power;
    gen_detector_eta    = detector_eta;
    
    gen_param_tree->Fill();
    gen_param_tree->Write("", TObject::kOverwrite);
    
    // Output Variables
    
    // Initializes arrays for particle info. DO NOT USE DOUBLES!!!
    // Must hard-code numbers for maximum compatibility
    int     particle_n;             // particle index
    float   particle_pt[200];       // transverse momentum
    float   particle_phi[200];      // azimuthal angle
    float   particle_eta[200];      // pseudorapidity
    int     particle_pid[200];      // numerical code for particle type
    int     particle_m[200];        // particle mass
    
    int     jet_n;
    float   jet_pt[100];
    float   jet_y[100];
    float   jet_phi[100];
    float   jet_mass[100];
    float   jet_area[100];
    float   jet_area_err[100];
    int     jet_const_n[100];
    float   jet_const_pt[100][400];
    float   jet_const_eta[100][400];
    float   jet_const_phi[100][400];
    
    // Builds particle and jet branches
    // Note: Can only do variable size for first variable in array
    pythia_tree->Branch("particle_n",   &particle_n);
    pythia_tree->Branch("particle_pt",  particle_pt,    "particle_pt[particle_n]/F");
    pythia_tree->Branch("particle_eta", particle_eta,   "particle_eta[particle_n]/F");
    pythia_tree->Branch("particle_phi", particle_phi,   "particle_phi[particle_n]/F");
    pythia_tree->Branch("particle_pid", particle_pid,   "particle_pid[particle_n]/I");
    pythia_tree->Branch("particle_m",   particle_m,     "particle_m[particle_n]/I");
    
    jet_tree->Branch("jet_n",           &jet_n);
    jet_tree->Branch("jet_pt",          jet_pt,         "jet_pt[jet_n]/F");
    jet_tree->Branch("jet_y",           jet_y,          "jet_y[jet_n]/F");
    jet_tree->Branch("jet_phi",         jet_phi,        "jet_phi[jet_n]/F");
    jet_tree->Branch("jet_mass",        jet_mass,       "jet_mass[jet_n]/F");
    jet_tree->Branch("jet_area",        jet_area,       "jet_area[jet_n]/F");
    jet_tree->Branch("jet_area_err",    jet_area_err,   "jet_area_err[jet_n]/F");
    jet_tree->Branch("jet_const_n",     jet_const_n,    "jet_const_n[jet_n]/I");
    jet_tree->Branch("jet_const_pt",    jet_const_pt,   "jet_const_pt[jet_n][400]/F");
    jet_tree->Branch("jet_const_eta",   jet_const_eta,  "jet_const_eta[jet_n][400]/F");
    jet_tree->Branch("jet_const_phi",   jet_const_phi,  "jet_const_phi[jet_n][400]/F");
    
    // --- PYTHIA Setup ---
    // --- LHC process and output selection. Initialization.
    Pythia pythia;
    Settings& pythia_settings = pythia.settings;
    pythia_settings.parm("Beams:eCM", beam_power);  // [GeV]
    pythia_settings.parm("PhaseSpace:pTHatMin", pt_hat_min);
    pythia_settings.parm("PhaseSpace:pTHatMax", pt_hat_max);
    pythia.readString("HardQCD:all = on");  // Turns on hard scattering
    if ( pt_bias_power != 0. ) {
        pythia_settings.parm("PhaseSpace:bias2SelectionPow", pt_bias_power);
        pythia.readString("PhaseSpace:bias2Selection = on");
        pythia.readString("PhaseSpace:bias2SelectionRef = 100.");
    }
    pythia.init();
    
    if ( print_out ) std::cout << "PYTHIA settings initialized." << std::endl;
    
    // --- FastJet Setup ---
    // --- Cluster each PYTHIA event
    fastjet::JetDefinition jet_definition(fastjet::antikt_algorithm, fastjet_radius);
    fastjet::AreaDefinition area_definition;
    
    if (!d_use_voronoi) {
        double ghost_eta_max = detector_eta;
        double ghost_area    = 0.05;
        int    active_area_repeats = 5;
        fastjet::GhostedAreaSpec ghost_spec(ghost_eta_max, active_area_repeats, ghost_area);
        area_definition = fastjet::AreaDefinition(fastjet::active_area, ghost_spec);
    } else {
        double effective_Rfact = 1.0;
        area_definition = fastjet::VoronoiAreaSpec(effective_Rfact);
    }
    
    if ( print_out ) std::cout << "FastJet settings initialized." << std::endl;
    
    // Loop to generate PYTHIA events. Skip if error.
    int e = 0;
    while ( e < event_count ) {
        // ---- PYTHIA -----
        if (!pythia.next()) continue;

        particle_n = 0;

        for ( int p = 0; p < pythia.event.size(); p++ ) {
            if (!pythia.event[p].isFinal()) continue;                       // Ends if final particle
            if (!pythia.event[p].isCharged()) continue;                     // Skips neutral particles
            if (pythia.event[p].pT() < 0.15) continue;                      // Skips particles with pT < 0.15
            if (fabs(pythia.event[p].eta()) > detector_eta) continue;       // Skips particles with rapidity outside the detector rapidity

            particle_pt[particle_n]  = pythia.event[p].pT();
            particle_eta[particle_n] = pythia.event[p].eta();
            particle_phi[particle_n] = pythia.event[p].phi();
            particle_pid[particle_n] = pythia.event[p].id();
            particle_m[particle_n]   = -1 * pythia.event[p].m();        // Outputs mass as negative?

            particle_n++;
        }
        
        // ---------- FASTJET ----------
        
        TLorentzVector input_vector;
        std::vector<fastjet::PseudoJet> input_particles;
        
        for ( int p = 0 ; p < particle_n ; p++) {
            input_vector.SetPtEtaPhiM( particle_pt[p], particle_eta[p], particle_phi[p], particle_m[p] );
            input_particles.push_back( fastjet::PseudoJet( input_vector.Px(), input_vector.Py(), input_vector.Pz(), input_vector.E() ) );
        }
        
        fastjet::ClusterSequenceArea jet_clusters( input_particles, jet_definition, area_definition );
        std::vector<fastjet::PseudoJet> jet_list = sorted_by_pt( jet_clusters.inclusive_jets(fastjet_pt_min) );
        
        // Checks jets in an event to see if the event should be accepted
        jet_n = jet_list.size();
        if ( jet_n == 0 ) continue;
        if ( jet_n >= 100 ) jet_n = 100;    // Prevents more jets from being recorded than space in jet_n vector
        
        int max_check_index = jet_n;
        if ( accept_jet_index != 0 ) max_check_index = accept_jet_index;
        
        float jet_pt_temp;
        float jet_y_temp;
        
        bool accept_event = false;
        
        for ( int j = 0 ; j < max_check_index ; j++ ) {
            jet_pt_temp = jet_list[j].pt();
            jet_y_temp = jet_list[j].rap();
            
            if ( (jet_pt_temp >= jet_pt_min) && (jet_pt_temp <= jet_pt_max) &&
                (jet_y_temp >= -jet_eta_max) && (jet_y_temp <= jet_eta_max) ) {
                accept_event = true;
                break;
            }
            else continue;
        }
        
        if ( accept_event ) {
            
            // Sets values for each jet in an accepted event
            for ( int j = 0 ; j < jet_n ; j++ ) {
                jet_pt[j]       = jet_list[j].pt();
                jet_y[j]        = jet_list[j].rap();
                jet_phi[j]      = jet_list[j].phi();
                jet_mass[j]     = jet_list[j].m();
                jet_area[j]     = jet_list[j].area();
                jet_area_err[j] = jet_list[j].area_error();
                
                std::vector<fastjet::PseudoJet> jet_constituents = jet_list[j].constituents();
                jet_const_n[j]  = jet_constituents.size();
                
                for ( int p = 0 ; p < jet_constituents.size() ; p++) {
                    jet_const_pt[j][p]  = jet_constituents[p].pt();
                    jet_const_eta[j][p] = jet_constituents[p].eta();
                    jet_const_phi[j][p] = jet_constituents[p].phi();
                }
            }
            
            // Fills trees
            pythia_tree ->Fill();
            jet_tree    ->Fill();
            e++;
            
            // Prints values for jets that have passed
            
            if ( print_out && ( e % print_every_x ) == 0 ) {
                std::cout << "== NEW EVENT # " << e << " ==" << std::endl;
                std::cout << "-- Summary of Particles -- " << std::endl;
                
                for ( int p = 0 ; p < particle_n ; p++ ) {
                    std::cout << " -> particle #" << particle_n << ", id = " << particle_pid[p] << ", pT/eta/phi = " <<  pythia.event[p].pT() << " / " <<  pythia.event[p].y() << " / " <<  pythia.event[p].phi() << std::endl;
                }
                
                std::cout << "== NEW EVENT # " << e << " ==" << std::endl;
                std::cout << "-- Summary of Jets -- " << std::endl;
                
                for ( int j = 0; j < jet_n; j++ ) {
                    std::cout << " -> jet " << j << ", pT / y / phi = " << jet_pt[j] << " / " << jet_y[j] << " / " << jet_phi[j] << std::endl;
                }
            }
        }
        else continue;
    }
    
    if ( print_out ) pythia.stat(); // Prints PYTHIA stats at the end

    pythia_tree ->Write("", TObject::kOverwrite);
    jet_tree    ->Write("", TObject::kOverwrite);
    std::cout << "Output trees written to." << std::endl;

    delete pythia_tree;
    delete jet_tree;
    
    output_file->Close();
    
    std::cout << "Files saved and closed." << std::endl;
    
    return output_file_path;
}



// ----- THERMAL GENERATOR -----



// Definition of fit functions:

double Func_Gaussian(double *x,double *par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    double fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

double Func_Power(double *x,double *par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[2]);
    double fitval = par[0] * pow(arg, par[1]);
    return fitval;
}

double Func_Flat(double *x,double *par) {
    double fitval = par[0] * 1;
    return fitval;
}

double Func_ModifiedHagedorn(double* x, double* par) {
    double arg1 = par[0] * pow(x[0], 2) / pow( (pow(x[0], 2) + pow(par[1], 2)), 0.5);
    double arg2 = pow( 1 + (x[0] /par[2]), par[3] );
    double arg3 = 1 / x[0];
    double fitval =  arg1 * arg2; // * arg3;
    return fitval;
}

// Definition of random generator functions:

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
    Thermal_Generator()
    Generates a defined number of random thermal particle backgrounds.
    Returns a character string pointer of the output file path.
 */
char* Thermal_Generator(
    char  output_base_name[100],            // Base/stem of output file name (no file extensions or prefixes).
    char  output_directory[400],            // Path to desired output directory.
    int   event_count,                      // Number of collision events to generate. The number of jets will be greater than this.
    float pt_bias_power,                    // Bias applied to the particle pT distribution.
    float jet_pt_min,                       // Minimum jet pT to accept. At least 1 jet per event will be greater than this.
    float jet_pt_max,                       // Maximum jet pT to accept. At least 1 jet per event will be less than this.
    int   thermal_mean = d_thermal_mean,    // Mean number of thermal particles.
    int   thermal_sigma = d_thermal_sigma,  // Standard deviation of the mean number of thermal particles for Gaussian distribution.
    float beam_power = d_beam_power,        // [GeV] Simulated beam power in the dector. Default defined by jet_ml_constants.h file.
    float detector_eta = d_detector_eta,    // Maximum rapidity of the detector. Default defined by jet_ml_constants.h file.
    float MH_par_1 = d_MH_par_1,
    float MH_par_2 = d_MH_par_2,
    float MH_par_3 = d_MH_par_3,
    float MH_par_4 = d_MH_par_4,
    float thermal_pt_max = d_thermal_pt_max
    ) {
    
    float thermal_norm = 1 / sqrt(2 * math_pi * pow(thermal_sigma,2)); // Sets normalization factor for gaussian
    std::cout << "Thermal Norm: " << thermal_norm << std::endl;
    
    char output_file_path[500];
    snprintf(output_file_path, 500, "%s/Ther_%s_B%.0f_%.0f_%.0f_N%i.root", output_directory, output_base_name, pt_bias_power, jet_pt_min, jet_pt_max, event_count);
    TFile* output_file = new TFile(output_file_path, "UPDATE");
    
    if (print_out) std::cout << "File Created: " << output_file << std::endl;
    TTree* output_tree = new TTree("Thermal_Tree","Tree of thermal particles as background for events");
    if (print_out) std::cout << "TTree Created: Thermal_Tree" << std::endl;
    
    // Generator Parameter Tree
    // This TTree stores the parameters used for generation in the file to preserve this data for reference
    TTree* gen_param_tree   = new TTree("Generator_Parameters","Flat tree of parameters used for file generation");
    
    int   gen_event_count;
    float gen_pt_bias_power;
    float gen_jet_pt_min;
    float gen_jet_pt_max;
    int   gen_gaus_mean;
    int   gen_gaus_sigma;
    float gen_gaus_norm;
    float gen_beam_power;
    float gen_detector_eta;
    
    gen_param_tree->Branch("event_count",       &gen_event_count);
    gen_param_tree->Branch("pt_bias_power",     &gen_pt_bias_power);
    gen_param_tree->Branch("jet_pt_min",        &gen_jet_pt_min);
    gen_param_tree->Branch("jet_pt_max",        &gen_jet_pt_max);
    gen_param_tree->Branch("gaus_mean",         &gen_gaus_mean);
    gen_param_tree->Branch("gaus_sigma",        &gen_gaus_sigma);
    gen_param_tree->Branch("gaus_norm",         &gen_gaus_norm);
    gen_param_tree->Branch("beam_power",        &gen_beam_power);
    gen_param_tree->Branch("detector_eta",      &gen_detector_eta);
    
    gen_event_count     = event_count;
    gen_pt_bias_power   = pt_bias_power;
    gen_jet_pt_min      = jet_pt_min;
    gen_jet_pt_max      = jet_pt_max;
    gen_gaus_mean       = thermal_mean;
    gen_gaus_sigma      = thermal_sigma;
    gen_gaus_norm       = thermal_norm;
    gen_beam_power      = beam_power;
    gen_detector_eta    = detector_eta;
    
    gen_param_tree->Fill();
    gen_param_tree->Write("", TObject::kOverwrite);
    
    // Output Variables
    
    TF1* thermal_n_func   = new TF1("thermal_n_func",   Func_Gaussian, 0, 2 * thermal_mean, 3);
    TF1* thermal_pt_func  = new TF1("thermal_pt_func",  Func_ModifiedHagedorn, 0, 100, 4);
    TF1* thermal_eta_func = new TF1("thermal_eta_func", Func_Flat, -detector_eta, detector_eta, 1);
    TF1* thermal_phi_func = new TF1("thermal_phi_func", Func_Flat, -math_pi, math_pi, 1);

    thermal_n_func   ->SetParameters(thermal_norm, thermal_mean, thermal_sigma);
    thermal_pt_func  ->SetParameters(MH_par_1, MH_par_2, MH_par_3, MH_par_4);
    thermal_eta_func ->SetParameter(0, 1.);
    thermal_phi_func ->SetParameter(0, 1.);
    
    // Initializes arrays for particle info
    int   particle_n;         // particle index
    float particle_pt[4000];  // transverse momentum
    float particle_phi[4000]; // azimuthal angle
    float particle_eta[4000]; // pseudorapidity
    
    // Builds TTree branches
    output_tree->Branch("particle_n",   &particle_n);
    output_tree->Branch("particle_pt",  particle_pt,     "particle_pt[particle_n]/F");
    output_tree->Branch("particle_eta", particle_eta,    "particle_eta[particle_n]/F");
    output_tree->Branch("particle_phi", particle_phi,    "particle_phi[particle_n]/F");
    if (print_out) std::cout << "TTree branches built." <<   std::endl;
    
    for (int e = 0 ; e < event_count ; e++) {
        // Initializes particle_n to 0 for the event
        particle_n = 0;
        
        // Generates random number of particles for thermal background
        int nParticles = Generate_N(thermal_n_func, thermal_mean);
        
        if ( print_out && (e % print_every_x) == 0 ) {
            std::cout << "----- Printing Event #" << e << " -----" << std::endl;
        }
        else if ( (e % print_every_x) == 0 ) {
            std::cout << "Generating Event #" << e << "..." << std::endl;
        }
        
        // For each particle, assigns pt, eta, phi, then increments n
        for (int p = 0 ; p < nParticles ; p++) {
            particle_pt[particle_n]  = Generate_Pt(thermal_pt_func, thermal_pt_max);
            particle_eta[particle_n] = Generate_Eta(thermal_eta_func, detector_eta);
            particle_phi[particle_n] = Generate_Phi(thermal_phi_func);
            
            if ( print_out && (e % print_every_x) == 0 ) {
                if ( p % 100 == 0 ) {
                    std::cout << " -> particle #" << particle_n << ", pT/eta/phi = " <<  particle_pt[particle_n] << " / " <<  particle_eta[particle_n] << " / " <<  particle_phi[particle_n] << std::endl;
                }
            }
            
            // Increments particle_n for the next particle
            particle_n++;
        }
        output_tree->Fill();
        if ( print_out && (e % print_every_x) == 0 ) std::cout << "Event #" << e << " TTree filled." << std::endl;
    }
    
    output_tree->Write("", TObject::kOverwrite);
    
    delete output_tree;
    
    output_file->Close();
    
    if ( print_out ) std::cout << "File written to and closed." << std::endl;
    
    return output_file_path;
}



// ----- COMBINE EVENTS -----



/*
    Combine_Events()
    Combines generated PYTHIA and thermal particle into combined events.
    Returns a character string pointer of the output file path.
 */
char* Combine_Events(
    char  output_base_name[100],            // Base/stem of output file name (no file extensions or prefixes).
    char  output_directory[400],            // Path to desired output directory.
    char  pythia_file_path[500],            // File path of input PYTHIA events
    char  thermal_file_path[500]            // File path of input thermal background events
    ) {
    
    // Open PYTHIA events data and access thermal TTree
    TFile* pythia_file = new TFile(pythia_file_path, "READ");
    std::cout << "Reading PYTHIA File" << std::endl;
    
    TTree* pythia_tree = (TTree*) pythia_file->Get("Pythia_Tree");
    TTree* pythia_param_tree = (TTree*) pythia_file->Get("Generator_Parameters");
    
    int   pythia_n;
    float pythia_pt[200];
    float pythia_eta[200];
    float pythia_phi[200];
    
    pythia_tree->SetBranchAddress("particle_n",   &pythia_n);
    pythia_tree->SetBranchAddress("particle_pt",  pythia_pt);
    pythia_tree->SetBranchAddress("particle_eta", pythia_eta);
    pythia_tree->SetBranchAddress("particle_phi", pythia_phi);
    
    int   p_gen_event_count;
    float p_gen_pt_bias_power;
    float p_gen_jet_pt_min;
    float p_gen_jet_pt_max;
    float p_gen_jet_eta_max;
    float p_gen_pt_hat_min;
    float p_gen_pt_hat_max;
    float p_gen_fastjet_radius;
    float p_gen_fastjet_pt_min;
    float p_gen_beam_power;
    float p_gen_detector_eta;
    
    pythia_param_tree->SetBranchAddress("event_count",      &p_gen_event_count);
    pythia_param_tree->SetBranchAddress("pt_bias_power",    &p_gen_pt_bias_power);
    pythia_param_tree->SetBranchAddress("jet_pt_min",       &p_gen_jet_pt_min);
    pythia_param_tree->SetBranchAddress("jet_pt_max",       &p_gen_jet_pt_max);
    pythia_param_tree->SetBranchAddress("jet_eta_max",      &p_gen_jet_eta_max);
    pythia_param_tree->SetBranchAddress("pt_hat_min",       &p_gen_pt_hat_min);
    pythia_param_tree->SetBranchAddress("pt_hat_max",       &p_gen_pt_hat_max);
    pythia_param_tree->SetBranchAddress("fastjet_radius",   &p_gen_fastjet_radius);
    pythia_param_tree->SetBranchAddress("fastjet_pt_min",   &p_gen_fastjet_pt_min);
    pythia_param_tree->SetBranchAddress("beam_power",       &p_gen_beam_power);
    pythia_param_tree->SetBranchAddress("detector_eta",     &p_gen_detector_eta);
    
    pythia_param_tree->GetEntry(0);
    
    if (print_out) std::cout << "TTree Accessed: Pythia_Tree" << std::endl;
    
    // Open thermal events data and access thermal TTree
    TFile* thermal_file = new TFile(thermal_file_path, "READ");
    std::cout << "Reading Thermal File" << std::endl;
    TTree* thermal_tree = (TTree*) thermal_file->Get("Thermal_Tree");
    TTree* thermal_param_tree = (TTree*) thermal_file->Get("Generator_Parameters");
    
    int   thermal_n;
    float thermal_pt[4000];
    float thermal_eta[4000];
    float thermal_phi[4000];
    
    thermal_tree->SetBranchAddress("particle_n",   &thermal_n);
    thermal_tree->SetBranchAddress("particle_pt",  thermal_pt);
    thermal_tree->SetBranchAddress("particle_eta", thermal_eta);
    thermal_tree->SetBranchAddress("particle_phi", thermal_phi);
    
    int   t_gen_gaus_mean;
    int   t_gen_gaus_sigma;
    float t_gen_gaus_norm;
    
    thermal_param_tree->SetBranchAddress("gaus_mean",       &t_gen_gaus_mean);
    thermal_param_tree->SetBranchAddress("gaus_sigma",      &t_gen_gaus_sigma);
    thermal_param_tree->SetBranchAddress("gaus_norm",       &t_gen_gaus_norm);
    
    if (print_out) std::cout << "TTree Accessed: Thermal_Tree" << std::endl;
    
    // Create new combined events data file
    char combined_file_path[500];
    pythia_param_tree->GetEntry(0);
    snprintf(combined_file_path, 500, "%s/Comb_%s_B%.0f_%.0f_%.0f_N%i.root", output_directory, output_base_name, p_gen_pt_bias_power, p_gen_jet_pt_min, p_gen_jet_pt_max, p_gen_event_count);
    TFile* combined_file = new TFile(combined_file_path, "UPDATE");
    std::cout << "Combined File Created" << std::endl;
    
    // Create new combined events TTree
    TTree* combined_tree = new TTree("Combined_Tree","Tree containing jet and thermal background particles.");
    TTree* combined_param_tree = new TTree("Generator_Parameters", "Flat tree of parameters (PYTHIA and Thermal) used for event generation");
    
    // Creates flat TTree of generator parameters
    int   gen_event_count;
    float gen_pt_bias_power;
    float gen_jet_pt_min;
    float gen_jet_pt_max;
    float gen_jet_eta_max;
    float gen_pt_hat_min;
    float gen_pt_hat_max;
    float gen_fastjet_radius;
    float gen_fastjet_pt_min;
    float gen_beam_power;
    float gen_detector_eta;
    int   gen_gaus_mean;
    float gen_gaus_sigma;
    float gen_gaus_norm;
    
    combined_param_tree->Branch("event_count",       &gen_event_count);
    combined_param_tree->Branch("pt_bias_power",     &gen_pt_bias_power);
    combined_param_tree->Branch("jet_pt_min",        &gen_jet_pt_min);
    combined_param_tree->Branch("jet_pt_max",        &gen_jet_pt_max);
    combined_param_tree->Branch("jet_eta_max",       &gen_jet_eta_max);
    combined_param_tree->Branch("pt_hat_min",        &gen_pt_hat_min);
    combined_param_tree->Branch("pt_hat_max",        &gen_pt_hat_max);
    combined_param_tree->Branch("fastjet_radius",    &gen_fastjet_radius);
    combined_param_tree->Branch("fastjet_pt_min",    &gen_fastjet_pt_min);
    combined_param_tree->Branch("beam_power",        &gen_beam_power);
    combined_param_tree->Branch("detector_eta",      &gen_detector_eta);
    combined_param_tree->Branch("gaus_mean",         &gen_gaus_mean);
    combined_param_tree->Branch("gaus_sigma",        &gen_gaus_sigma);
    combined_param_tree->Branch("gaus_norm",         &gen_gaus_norm);
    
    pythia_param_tree->GetEntry(0);
    thermal_param_tree->GetEntry(0);
    
    gen_event_count     = p_gen_event_count;
    gen_pt_bias_power   = p_gen_pt_bias_power;
    gen_jet_pt_min      = p_gen_jet_pt_min;
    gen_jet_pt_max      = p_gen_jet_pt_max;
    gen_jet_eta_max     = p_gen_jet_eta_max;
    gen_pt_hat_min      = p_gen_pt_hat_min;
    gen_pt_hat_max      = p_gen_pt_hat_max;
    gen_fastjet_radius  = p_gen_fastjet_radius;
    gen_fastjet_pt_min  = p_gen_fastjet_pt_min;
    gen_gaus_mean       = t_gen_gaus_mean;
    gen_gaus_sigma      = t_gen_gaus_sigma;
    gen_gaus_norm       = t_gen_gaus_norm;
    gen_beam_power      = p_gen_beam_power;
    gen_detector_eta    = p_gen_detector_eta;
    
    combined_param_tree->Fill();
    combined_param_tree->Write("", TObject::kOverwrite);
    
    if (print_out) std::cout << "TTree Created: Combined_Tree" << std::endl;
    
    // Build out combined TTree
    int   combined_n;
    float combined_pt[4200];
    float combined_eta[4200];
    float combined_phi[4200];
    int   combined_jet_class;
    
    combined_tree->Branch("particle_n",         &combined_n);
    combined_tree->Branch("particle_pt",        combined_pt,    "particle_pt[particle_n]/F");
    combined_tree->Branch("particle_eta",       combined_eta,   "particle_eta[particle_n]/F");
    combined_tree->Branch("particle_phi",       combined_phi,   "particle_phi[particle_n]/F");
    combined_tree->Branch("particle_jet_class", &combined_jet_class);
    
    if (print_out) std::cout << "Combined TTree branches built." <<   std::endl;
    
    // Merge particles into a single TTree for analysis
    for (int e = 0 ; e < pythia_tree->GetEntries() ; e++) {
        pythia_tree->GetEntry(e);
        thermal_tree->GetEntry(e);
        
        int pythia_count = pythia_n;
        int thermal_count = thermal_n;
        
        if ( print_out && (e % print_every_x) == 0 ) std::cout << "Pythia Count: " << pythia_count << ", Therm Count: " << thermal_count << std::endl;
        
        combined_n = 0;
                
        // Loop to add jet particles
        for (int p = 0 ; p < pythia_count ; p++) {
            combined_n += 1;
            combined_pt[p] = pythia_pt[p];
            combined_eta[p] = pythia_eta[p];
            combined_phi[p] = pythia_phi[p];
            combined_jet_class = 1;
            if ( debug && (e % print_every_x) == 0 ) std::cout << "source / pT / eta / phi: " << combined_jet_class << combined_pt[p] << combined_eta[p] << combined_phi[p] << std::endl;
        }
        
        // Loop to add thermal particles
        for (int t = 0 ; t < thermal_count ; t++) {
            combined_n += 1;
            combined_pt[pythia_count + t - 1]  = thermal_pt[t];
            combined_eta[pythia_count + t - 1] = thermal_eta[t];
            combined_phi[pythia_count + t - 1] = thermal_phi[t];
            combined_jet_class = 0;
            if ( debug && (e % print_every_x) == 0 ) std::cout << "source / pT / eta / phi: " << combined_jet_class << combined_pt[t] << combined_eta[t] << combined_phi[t] << std::endl;
        }
        
        if ( print_out && (e % print_every_x) == 0 ) std::cout << "Event #" << e << " has " << combined_n << " total particles." << std::endl;
        
        // Fills tree
        combined_tree->Fill();
        if ( print_out && (e % print_every_x) == 0 ) std::cout << "Event #" << e << " TTree filled." << std::endl;
    }
    combined_tree->Write("", TObject::kOverwrite);
    
    delete pythia_tree;
    delete thermal_tree;
    delete combined_tree;
    
    // Closes all files
    combined_file   ->Close();
    pythia_file     ->Close();
    thermal_file    ->Close();
    
    if ( print_out ) std::cout << "File written to and closed." << std::endl;
    
    return combined_file_path;
}



// ----- COMBINED EVENT JET CLUSTERER -----



/*
    Combined_Jet_Clusterer()
    Clusters particles from combined events into jets using anti-kt.
 */
void Combined_Jet_Clusterer(
    char   combined_file_path[500]      // File path of combined event file to use as input.
    ) {
    // Input Variables and Trees
    TFile* input_file = new TFile(combined_file_path, "UPDATE");
    TTree* input_tree = (TTree*) input_file->Get("Combined_Tree");
    
    int   input_n;
    float input_pt[4200];
    float input_eta[4200];
    float input_phi[4200];
    TLorentzVector input_vector;

    input_tree->SetBranchAddress("particle_n",      &input_n);
    input_tree->SetBranchAddress("particle_pt",     input_pt);
    input_tree->SetBranchAddress("particle_eta",    input_eta);
    input_tree->SetBranchAddress("particle_phi",    input_phi);
    
    TTree* input_param_tree = (TTree*) input_file->Get("Generator_Parameters");
    
    float gen_fastjet_radius;
    float gen_fastjet_pt_min;
    float gen_beam_power;
    float gen_detector_eta;
    
    input_param_tree->SetBranchAddress("fastjet_radius",    &gen_fastjet_radius);
    input_param_tree->SetBranchAddress("fastjet_pt_min",    &gen_fastjet_pt_min);
    input_param_tree->SetBranchAddress("beam_power",        &gen_beam_power);
    input_param_tree->SetBranchAddress("detector_eta",      &gen_detector_eta);
    
    input_param_tree->GetEntry(0);
    
    float fastjet_radius = gen_fastjet_radius;
    float fastjet_pt_min = gen_fastjet_pt_min;
    float detector_eta   = gen_detector_eta;
    
    std::cout << fastjet_radius << ", " << fastjet_pt_min << std::endl;
    
    if (print_out) std::cout << "Input file and tree successfully accessed." << std::endl;
    
    // Defines the jet definition and jet area definition
    fastjet::JetDefinition jet_definition(fastjet::antikt_algorithm, gen_fastjet_radius);
    fastjet::AreaDefinition area_definition;
    
    if (!d_use_voronoi) {
        double ghost_etamax = detector_eta;
        double ghost_area    = 0.05;
        int    active_area_repeats = 5;
        fastjet::GhostedAreaSpec ghost_spec(ghost_etamax, active_area_repeats, ghost_area);
        area_definition = fastjet::AreaDefinition(fastjet::active_area, ghost_spec);
    } else {
        double effective_Rfact = 1.0;
        area_definition = fastjet::VoronoiAreaSpec(effective_Rfact);
    }
    
    if ( print_out ) std::cout << "FastJet settings initialized." << std::endl;
        
    // Sets jet data to process
    TTree* output_tree = new TTree("FastJet_Tree","Tree of jet clusters");
    
    int     jet_n;
    float   jet_pt[100];
    float   jet_y[100];
    float   jet_phi[100];
    float   jet_mass[100];
    float   jet_area[100];
    float   jet_area_err[100];
    int     jet_const_n[100];
    float   jet_const_pt[100][400];
    float   jet_const_eta[100][400];
    float   jet_const_phi[100][400];
    
    // Note: Can only do variable size for first variable in array
    output_tree->Branch("jet_n",            &jet_n);
    output_tree->Branch("jet_pt",           jet_pt,         "jet_pt[jet_n]/F");
    output_tree->Branch("jet_y",            jet_y,          "jet_y[jet_n]/F");
    output_tree->Branch("jet_phi",          jet_phi,        "jet_phi[jet_n]/F");
    output_tree->Branch("jet_mass",         jet_mass,       "jet_mass[jet_n]/F");
    output_tree->Branch("jet_area",         jet_area,       "jet_area[jet_n]/F");
    output_tree->Branch("jet_area_err",     jet_area_err,   "jet_area_err[jet_n]/F");
    output_tree->Branch("jet_const_n",      jet_const_n,    "jet_const_n[jet_n]/I");
    output_tree->Branch("jet_const_pt",     jet_const_pt,   "jet_const_pt[jet_n][400]/F");
    output_tree->Branch("jet_const_eta",    jet_const_eta,  "jet_const_eta[jet_n][400]/F");
    output_tree->Branch("jet_const_phi",    jet_const_phi,  "jet_const_phi[jet_n][400]/F");
    
    if ( print_out ) std::cout << "Output tree variables initialized." << std::endl;
    
    // Loop to fill vector with particles
    for ( int e = 0 ; e < input_tree->GetEntries() ; e++ ) {
        
        if ( (e % print_every_x) == 0 ) std::cout << "Processing Event " << e << "." << std::endl;
        
        input_tree->GetEntry(e);
        std::vector<fastjet::PseudoJet> input_particles;
        
        for ( int p = 0 ; p < input_n ; p++) {
            input_vector.SetPtEtaPhiM(input_pt[p], input_eta[p], input_phi[p], m_pion);
            input_particles.push_back(fastjet::PseudoJet( input_vector.Px(), input_vector.Py(), input_vector.Pz(), input_vector.E() ));
        }
        
        fastjet::ClusterSequenceArea jet_clusters(input_particles, jet_definition, area_definition);
        std::vector<fastjet::PseudoJet> jet_list = sorted_by_pt(jet_clusters.inclusive_jets(fastjet_pt_min));
        
        jet_n = jet_list.size();
        if ( jet_n >= 100 ) continue;
        if ( print_out && ( e % print_every_x ) == 0 ) std::cout << "Event #" << e << " has " << jet_n << " jets." << std::endl;
        
        if ( jet_n == 0 ) continue;
        
        for ( int j = 0 ; j < jet_n ; j++ ) {
            jet_pt[j]   = jet_list[j].pt();
            jet_y[j]    = jet_list[j].rap();
            jet_phi[j]  = jet_list[j].phi();
            jet_mass[j] = jet_list[j].m();
            jet_area[j] = jet_list[j].area();
            jet_area_err[j] = jet_list[j].area_error();
            
            std::vector<fastjet::PseudoJet> jet_constituents = jet_list[j].constituents();
            jet_const_n[j]  = jet_constituents.size();
                        
            for ( int p = 0 ; p < jet_const_n[j] ; p++) {
                jet_const_pt[j][p]   = jet_constituents[p].pt();
                jet_const_eta[j][p]  = jet_constituents[p].eta();
                jet_const_phi[j][p]  = jet_constituents[p].phi();
            }
        }
        
        // Prints values for jets
        if ( print_out && ( e % print_every_x ) == 0 ) {
            std::cout << "== NEW EVENT # " << e << " ==" << std::endl;
            std::cout << "-- Summary of Jets -- " << std::endl;
            for ( int j = 0; j < jet_n; j++ ) {
                std::cout << " -> jet " << j << ", pT / y / phi = " << jet_list[j].pt() << " / " << jet_list[j].rap() << " / " << jet_list[j].phi() << std::endl;
            }
        }
        
        output_tree->Fill();
        if ( print_out && ( e % print_every_x ) == 0 ) std::cout << "Event #" << e << " TTree filled." << std::endl;
    }
    output_tree->Write("", TObject::kOverwrite);
    std::cout << "Output file written to." << std::endl;
    
    delete output_tree;
    input_file->Close();
    delete input_file;
        
    std::cout << "Files closed." << std::endl;
}



// ----- MACHINE LEARNING FILE PREP ----



/*
    Check_File_Compatibility()
    Used to check that two files (e.g. a PYTHIA jet file and a combined event jet file) are compatible.
    ONLY checks event count, pT bias power, and jet pT min/max. File path order does not matter.
 */
bool Check_File_Compatibility(
    char  file_path_A[500],     // File path of first file to compare.
    char  file_path_B[500]      // File path of first file to compare.
    ) {
    
    bool files_ok = true;
    
    TFile* file_A = new TFile(file_path_A, "READ");
    TTree* tree_A;
    if (file_A->GetListOfKeys()->Contains("Generator_Parameters")) {
        tree_A = (TTree*) file_A->Get("Generator_Parameters"); }
    else files_ok = false;
        
    TFile* file_B = new TFile(file_path_B, "READ");
    TTree* tree_B;
    if (file_B->GetListOfKeys()->Contains("Generator_Parameters")) {
        tree_B = (TTree*) file_B->Get("Generator_Parameters"); }
    else files_ok = false;
    
    if (files_ok) {
        int   a_event_count;
        float a_pt_bias_power;
        float a_jet_pt_min;
        float a_jet_pt_max;
        
        tree_A->SetBranchAddress("event_count",       &a_event_count);
        tree_A->SetBranchAddress("pt_bias_power",     &a_pt_bias_power);
        tree_A->SetBranchAddress("jet_pt_min",        &a_jet_pt_min);
        tree_A->SetBranchAddress("jet_pt_max",        &a_jet_pt_max);
        
        int   b_event_count;
        float b_pt_bias_power;
        float b_jet_pt_min;
        float b_jet_pt_max;
        
        tree_B->SetBranchAddress("event_count",       &b_event_count);
        tree_B->SetBranchAddress("pt_bias_power",     &b_pt_bias_power);
        tree_B->SetBranchAddress("jet_pt_min",        &b_jet_pt_min);
        tree_B->SetBranchAddress("jet_pt_max",        &b_jet_pt_max);
        
        tree_A->GetEntry(0);
        tree_B->GetEntry(0);
        
        // Checks that files use the same generator data
        if ((a_event_count   != b_event_count)   ||
            (a_pt_bias_power != b_pt_bias_power) ||
            (a_jet_pt_min    != b_jet_pt_min)    ||
            (a_jet_pt_max    != b_jet_pt_max) ) {
            files_ok = false;
            std::cout << "Input files are incompatible! Generator data is not identical between files." << std::endl;
            }
        }
        
    return files_ok;
    }

/*
    Make_Flat_Jet_Distribution()
    --- WORK IN PROGRESS ---
    Should flatten a TTree based on the lowest number of jets in a pT range.
    Need to modify to use random generator to throw out jets
 */
/*
void Make_Flat_Jet_Distribution(
    char  mlprep_file_path[500],        // Full name of ML_Prep root file
    char  input_tree_name[100],
    float bin_width,                    // [GeV]
    float jet_pt_min,
    float jet_pt_max
    ) {
    
    TFile* file = new TFile(mlprep_file_path, "UPDATE");
    TTree* input_tree = (TTree*) file->Get(input_tree_name);
    
    float jet_pt_true;
    input_tree->SetBranchAddress("jet_pt_true", &jet_pt_true);
    
    char output_tree_name[100];
    snprintf(output_tree_name, 100, "%s_Flat", input_tree_name);
    char output_tree_desc[100];
    snprintf(output_tree_desc, 100, "Flattened version of %s", input_tree_name);
    TTree* output_tree = input_tree->CloneTree(0);
    output_tree->SetObject(output_tree_name, output_tree_desc);
    
    const int   pt_bins = int(jet_pt_max - jet_pt_min);             // Assumes 1 GeV wide bins
    const float pt_bin_width = (jet_pt_max - jet_pt_min)/pt_bins;   // Sets bin width
    int   jet_pt_bin_count[pt_bins] = {0};        // Array of number of jets in each bin
    float jet_pt_bin_upper[pt_bins];        // Array of upper bound of each array
    int   bin_skip_array[pt_bins] = {0};    // Array of integers, where every Nth jet in the corresponding array index is skipped to flatten
    int   bin_counter[pt_bins] = {0};       // Keeps track of the number of jets in each bin
    
    // Sets upper bounds for each bin
    for ( int i = 0 ; i < pt_bins ; i++ ) {
        jet_pt_bin_upper[i] = i * pt_bin_width + jet_pt_min;
    }
    
    // Counts number of entries in each bin
    for ( int j = 0 ; j < input_tree->GetEntries() ; j++ ) {
        input_tree->GetEntry(j);
        for ( int i = 0 ; i < pt_bins ; i++ ) {
            if ( jet_pt_true < jet_pt_bin_upper[i] ) {
                jet_pt_bin_count[i]++;
                break;
            }
        }
    }
    
    // Sets the Nth entry to skip for each bin to flatten the array
    int min_bin_count = std::min_element(std::begin(jet_pt_bin_count), std::end(jet_pt_bin_count));
    for ( int i = 0 ; i < pt_bins ; i++ ) {
        bin_skip_array[i] = int(std::nearbyint( jet_pt_bin_count[i] / (jet_pt_bin_array[i] - min_bin_count) ));
    }
    
    // Copies set number of entries for each bin
    for ( int j = 0 ; j < input_tree->GetEntries() ; j++ ) {
        input_tree->GetEntry(j);
        for ( int i = 0 ; i < pt_bins ; i++ ) {
            if ( jet_pt_true < jet_pt_bin_upper[i] ) {
                bin_counter[i]++;
                if ( bin_counter[i] % bin_skip_array[i] == 0 ) break;
                output_tree->Fill();
            }
        }
    }
    
    output_tree->Write("", TObject::kOverwrite);
    file->Write("", TObject::kOverwrite);
    
    delete output_tree;
}
*/

/*
    Jet_ML_Prep()
    Prepares a file for machine learning.
 */
void Jet_ML_Prep(
    char  output_base_name[100],        // Base/stem of output file name (no file extensions or prefixes).
    char  output_directory[400],        // Path to desired output directory.
    char  combined_file_path[500],      // Full name of combined output file
    char  pythia_file_path[500],        // Full name of combined output file
    float jet_pt_min = -1.,             // If >0, overrides min jet pt accepted. If <0 uses min from input file parameters
    float jet_pt_max = -1.,             // If >0, overrides max jet pt accepted. If <0 uses max from input file parameters
    float jet_eta_max = d_jet_eta_max,  // Highest rapidity to accept jets
    int   softest_jet_index = d_softest_jet_index,          // Only accepts the top [softest_jet_index] jets. If softest_jet_index = 0, accepts all jets
    float fastjet_match_radius = d_fastjet_match_radius     // Two jets are matched if they are within this radius squared.
    ) {
    
    if ( !Check_File_Compatibility(combined_file_path, pythia_file_path) ) return;
    
    // Combined file of jets from thermal and PYTHIA data
    TFile* combined_file = new TFile(combined_file_path, "READ");
    TTree* comb_fastjet_tree = (TTree*) combined_file->Get("FastJet_Tree");
    
    std::cout << "Reading combined file." << std::endl;
    
    int     c_jet_n;
    float   c_jet_pt[100];
    float   c_jet_y[100];
    float   c_jet_phi[100];
    float   c_jet_mass[100];
    float   c_jet_area[100];
    float   c_jet_area_err[100];
    int     c_jet_const_n[100];
    float   c_jet_const_pt[100][400];
    float   c_jet_const_eta[100][400];
    float   c_jet_const_phi[100][400];
    
    comb_fastjet_tree->SetBranchAddress("jet_n",            &c_jet_n);
    comb_fastjet_tree->SetBranchAddress("jet_pt",           c_jet_pt);
    comb_fastjet_tree->SetBranchAddress("jet_y",            c_jet_y);
    comb_fastjet_tree->SetBranchAddress("jet_phi",          c_jet_phi);
    comb_fastjet_tree->SetBranchAddress("jet_area",         c_jet_area);
    comb_fastjet_tree->SetBranchAddress("jet_area_err",     c_jet_area_err);
    comb_fastjet_tree->SetBranchAddress("jet_mass",         c_jet_mass);
    comb_fastjet_tree->SetBranchAddress("jet_const_n",      c_jet_const_n);
    comb_fastjet_tree->SetBranchAddress("jet_const_pt",     c_jet_const_pt);
    comb_fastjet_tree->SetBranchAddress("jet_const_eta",    c_jet_const_eta);
    comb_fastjet_tree->SetBranchAddress("jet_const_phi",    c_jet_const_phi);
    
    TTree* comb_param_tree = (TTree*) combined_file->Get("Generator_Parameters");
    
    int   c_gen_event_count;
    float c_gen_pt_bias_power;
    float c_gen_jet_pt_min;
    float c_gen_jet_pt_max;
    float c_gen_jet_eta_max;
    float c_gen_pt_hat_min;
    float c_gen_pt_hat_max;
    float c_gen_fastjet_radius;
    float c_gen_fastjet_pt_min;
    int   c_gen_gaus_mean;
    float c_gen_gaus_sigma;
    float c_gen_gaus_norm;
    float c_gen_beam_power;
    float c_gen_detector_eta;
    
    comb_param_tree->SetBranchAddress("event_count",        &c_gen_event_count);
    comb_param_tree->SetBranchAddress("pt_bias_power",      &c_gen_pt_bias_power);
    comb_param_tree->SetBranchAddress("jet_pt_min",         &c_gen_jet_pt_min);
    comb_param_tree->SetBranchAddress("jet_pt_max",         &c_gen_jet_pt_max);
    comb_param_tree->SetBranchAddress("jet_eta_max",        &c_gen_jet_eta_max);
    comb_param_tree->SetBranchAddress("pt_hat_min",         &c_gen_pt_hat_min);
    comb_param_tree->SetBranchAddress("pt_hat_max",         &c_gen_pt_hat_max);
    comb_param_tree->SetBranchAddress("fastjet_radius",     &c_gen_fastjet_radius);
    comb_param_tree->SetBranchAddress("fastjet_pt_min",     &c_gen_fastjet_pt_min);
    comb_param_tree->SetBranchAddress("gaus_mean",          &c_gen_gaus_mean);
    comb_param_tree->SetBranchAddress("gaus_sigma",         &c_gen_gaus_sigma);
    comb_param_tree->SetBranchAddress("gaus_norm",          &c_gen_gaus_norm);
    comb_param_tree->SetBranchAddress("beam_power",         &c_gen_beam_power);
    comb_param_tree->SetBranchAddress("detector_eta",       &c_gen_detector_eta);
    
    // File of jets from PYTHIA
    TFile* pythia_file = new TFile(pythia_file_path, "READ");
    TTree* pyth_fastjet_tree = (TTree*) pythia_file->Get("FastJet_Tree");
    
    std::cout << "Reading PYTHIA file." << std::endl;
    
    int     p_jet_n;
    float   p_jet_pt[100];
    float   p_jet_y[100];
    float   p_jet_phi[100];
    float   p_jet_mass[100];
    float   p_jet_area[100];
    float   p_jet_area_err[100];
    int     p_jet_const_n[100];
    float   p_jet_const_pt[100][400];
    float   p_jet_const_eta[100][400];
    float   p_jet_const_phi[100][400];
    
    pyth_fastjet_tree->SetBranchAddress("jet_n",         &p_jet_n);
    pyth_fastjet_tree->SetBranchAddress("jet_pt",        p_jet_pt);
    pyth_fastjet_tree->SetBranchAddress("jet_y",         p_jet_y);
    pyth_fastjet_tree->SetBranchAddress("jet_phi",       p_jet_phi);
    pyth_fastjet_tree->SetBranchAddress("jet_area",      p_jet_area);
    pyth_fastjet_tree->SetBranchAddress("jet_area_err",  p_jet_area_err);
    pyth_fastjet_tree->SetBranchAddress("jet_mass",      p_jet_mass);
    pyth_fastjet_tree->SetBranchAddress("jet_const_n",   p_jet_const_n);
    pyth_fastjet_tree->SetBranchAddress("jet_const_pt",  p_jet_const_pt);
    pyth_fastjet_tree->SetBranchAddress("jet_const_eta", p_jet_const_eta);
    pyth_fastjet_tree->SetBranchAddress("jet_const_phi", p_jet_const_phi);
    
    comb_param_tree->GetEntry(0);
    
    char output_file_params[200];
    snprintf(output_file_params, 200, "%s_B%.0f_%.0f_%.0f_N%i", output_base_name, c_gen_pt_bias_power, c_gen_jet_pt_min, c_gen_jet_pt_max, c_gen_event_count);
    
    char output_file_path[200];
    snprintf(output_file_path, 200, "%s/ML_Prep_%s.root", output_directory, output_file_params);
    TFile* output_file = new TFile(output_file_path, "UPDATE");
    
    char output_tree_name[200];
    snprintf(output_tree_name, 200, "Jet_ML_%s", output_file_params);
    char output_tree_description[200];
    snprintf(output_tree_description, 200, "Generated with bias of pT^%.0f, %.1f - %.1f GeV, %i events, jet match radius of %.2f", c_gen_pt_bias_power, c_gen_jet_pt_min, c_gen_jet_pt_max, c_gen_event_count, fastjet_match_radius);
    TTree* output_tree = new TTree(output_tree_name, output_tree_description);
    
    std::cout << "----- Preparing " << output_tree_name << " -----" << std::endl;
    
    // Creates flat TTree of generator parameters
    TTree* output_param_tree = new TTree("ML_Prep_Parameters", "Flat tree of parameters used for event generation and ML preparation");
    
    std::cout << "Making output parameter tree." << std::endl;
    
    int   gen_event_count;
    float gen_pt_bias_power;
    float gen_jet_pt_min;
    float gen_jet_pt_max;
    float gen_jet_eta_max;
    float gen_pt_hat_min;
    float gen_pt_hat_max;
    float gen_fastjet_radius;
    float gen_fastjet_pt_min;
    float gen_fastjet_match_radius;
    float gen_beam_power;
    float gen_detector_eta;
    int   gen_gaus_mean;
    float gen_gaus_sigma;
    float gen_gaus_norm;
    
    output_param_tree->Branch("event_count",          &gen_event_count);
    output_param_tree->Branch("pt_bias_power",        &gen_pt_bias_power);
    output_param_tree->Branch("jet_pt_min",           &gen_jet_pt_min);
    output_param_tree->Branch("jet_pt_max",           &gen_jet_pt_max);
    output_param_tree->Branch("jet_eta_max",          &gen_jet_eta_max);
    output_param_tree->Branch("pt_hat_min",           &gen_pt_hat_min);
    output_param_tree->Branch("pt_hat_max",           &gen_pt_hat_max);
    output_param_tree->Branch("fastjet_radius",       &gen_fastjet_radius);
    output_param_tree->Branch("fastjet_pt_min",       &gen_fastjet_pt_min);
    output_param_tree->Branch("fastjet_match_radius", &gen_fastjet_match_radius);
    output_param_tree->Branch("beam_power",           &gen_beam_power);
    output_param_tree->Branch("detector_eta",         &gen_detector_eta);
    output_param_tree->Branch("gaus_mean",            &gen_gaus_mean);
    output_param_tree->Branch("gaus_sigma",           &gen_gaus_sigma);
    output_param_tree->Branch("gaus_norm",            &gen_gaus_norm);
    
    gen_event_count             = c_gen_event_count;
    gen_pt_bias_power           = c_gen_pt_bias_power;
    gen_jet_pt_min              = c_gen_jet_pt_min;
    gen_jet_pt_max              = c_gen_jet_pt_max;
    gen_jet_eta_max             = c_gen_jet_eta_max;
    gen_pt_hat_min              = c_gen_pt_hat_min;
    gen_pt_hat_max              = c_gen_pt_hat_max;
    gen_fastjet_radius          = c_gen_fastjet_radius;
    gen_fastjet_pt_min          = c_gen_fastjet_pt_min;
    gen_fastjet_match_radius    = fastjet_match_radius;
    gen_gaus_mean               = c_gen_gaus_mean;
    gen_gaus_sigma              = c_gen_gaus_sigma;
    gen_gaus_norm               = c_gen_gaus_norm;
    gen_beam_power              = c_gen_beam_power;
    gen_detector_eta            = c_gen_detector_eta;
    
    // Sets the pt_min and pt_max cuts - if not defined these will be established by the input files
    if ( jet_pt_min < 0 ) jet_pt_min = c_gen_jet_pt_min;
    if ( jet_pt_max < 0 ) jet_pt_max = c_gen_jet_pt_max;
    
    output_param_tree->Fill();
    output_param_tree->Write("", TObject::kOverwrite);
    
    std::cout << "Output parameter tree made." << std::endl;
    std::cout << "Jet pt min / max: " << c_gen_jet_pt_min << ", " << c_gen_jet_pt_max << std::endl;
    
    // Plotted Data Output Tree
    // NOTE: This assumes jets have already been sorted from highest E to lowest E by FastJet!
    
    int    nEvent = pyth_fastjet_tree->GetEntries();
    float  background_pt_median[nEvent];
    
    // ML X value variables (features)
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
    
    output_tree->Branch("jet_index",            &jet_index);
    output_tree->Branch("jet_pt_raw",           &jet_pt_raw,        "jet_pt_raw/F");
    output_tree->Branch("jet_pt_corr",          &jet_pt_corr,       "jet_pt_corr/F");
    output_tree->Branch("jet_pt_true",          &jet_pt_true,       "jet_pt_true/F");
    output_tree->Branch("jet_mass",             &jet_mass,          "jet_mass/F");
    output_tree->Branch("jet_area",             &jet_area,          "jet_area/F");
    output_tree->Branch("jet_area_err",         &jet_area_err,      "jet_area_err/F");
    output_tree->Branch("jet_const_n",          &jet_const_n,       "jet_const_n/I");
    output_tree->Branch("const_pt_mean",        &const_pt_mean,     "const_pt_mean/F");
    output_tree->Branch("const_pt_median",      &const_pt_median,   "const_pt_median/F");
    output_tree->Branch("const_1_pt",           &const_1_pt,        "const_1_pt/F");
    output_tree->Branch("const_2_pt",           &const_2_pt,        "const_2_pt/F");
    output_tree->Branch("const_3_pt",           &const_3_pt,        "const_3_pt/F");
    output_tree->Branch("const_4_pt",           &const_4_pt,        "const_4_pt/F");
    output_tree->Branch("const_5_pt",           &const_5_pt,        "const_5_pt/F");
    output_tree->Branch("const_6_pt",           &const_6_pt,        "const_6_pt/F");
    output_tree->Branch("const_7_pt",           &const_7_pt,        "const_7_pt/F");
    output_tree->Branch("const_8_pt",           &const_8_pt,        "const_8_pt/F");
    output_tree->Branch("const_9_pt",           &const_9_pt,        "const_9_pt/F");
    output_tree->Branch("const_10_pt",          &const_10_pt,       "const_10_pt/F");
    output_tree->Branch("jet_y",                &jet_y,             "jet_y/F");
    output_tree->Branch("jet_phi",              &jet_phi,           "jet_phi/F");
    output_tree->Branch("jet_rho",              &jet_rho,           "jet_rho/F");
    
    int truth_matched_counter = 0;
    int truth_unmatched_counter = 0;
    int truth_outside_pt_counter = 0;
    int truth_outside_y_counter = 0;
    float fastjet_match_rad_sq = pow(fastjet_match_radius, 2);
    
    std::cout << "Starting jet finding on combined tree." << std::endl;
    
    // Iterate through events to extract jets
    for ( int e = 0 ; e < comb_fastjet_tree->GetEntries() ; e++ ) {
        comb_fastjet_tree->GetEntry(e);
        pyth_fastjet_tree->GetEntry(e);
        
        if ( print_out && ( e % print_every_x ) == 0 ) std::cout << "Getting jets from event " << e << std::endl;
        
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
        
        // Loops through jets for machine learning predictors
        if ( print_out && ( e % print_every_x ) == 0 ) std::cout << "Jets for event " << e << ": " << c_jet_n << std::endl;
        
        if ( softest_jet_index == 0 ) softest_jet_index = c_jet_n;
        
        for ( int k = 0 ; k < p_jet_n ; k++ ) {
            // Iterate through pythia jets to match jet pT_true
            int c_match = -1;
            
            int c_match_tracker[100] = {0};
            
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
            // Skips if jet can't be matched to a PYTHIA jet
            if ( c_match < 0 ) {
                truth_unmatched_counter++;
                if ( debug) std::cout << "Event: " << e << " Truth Jet: " << k << " has no match found." << std::endl;
                if ( debug) std::cout << "    p_phi: " << p_jet_phi[k] << ", p_y: " << p_jet_y[k] << std::endl;
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
            
            if ( print_out && jet_pt_true != 0 && (e % print_every_x) == 0 ) std::cout << e << "-" << c_match << ": Truth_PYTHIA Jet: " << jet_pt_true << " ----- " << std::endl;
            
            output_tree->Fill();
        }
    }
    
    char truth_match_results[200];
    snprintf(truth_match_results, 200, "Matched %i jets. %i jets unmatched. %i jets outside pt. %i jets outside y. Unmatched ratio is: %.4f", truth_matched_counter, truth_unmatched_counter, truth_outside_pt_counter, truth_outside_y_counter, float(truth_unmatched_counter)/(truth_matched_counter + truth_unmatched_counter) );
    std::cout << truth_match_results << std::endl;
    
    output_tree->Write("", TObject::kOverwrite);
    output_file->Write("", TObject::kOverwrite);
    
    delete output_tree;
    delete output_param_tree;
    
    std::cout << "Data written to file." << std::endl;
    
    pythia_file->Close();
    combined_file->Close();
    output_file->Close();
    
    delete pythia_file;
    delete combined_file;
    delete output_file;
    
    std::cout << "Files closed." << std::endl;
    
    std::cout << "----- Completed " << output_tree_name << " -----" << std::endl;
}



// ----- JET GENERATOR -----



/*
    Jet_Generator_Full()
    Simplifies the jet generation process by chaining above functions together with a single set of input parameters.
    Default values (defined in the header file) can be used.
 */
void Jet_Generator_Full(
    char  output_base_name[100],                // Base/stem of output file name (no file extensions or prefixes).
    char  output_directory[400],                // Path to desired output directory.
    int   event_count,                          // Number of collision events to generate. The number of jets will be greater than this.
    float pt_bias_power,                        // Bias applied to the particle pT distribution.
    float jet_pt_min,                           // Minimum jet pT to accept. At least 1 jet per event will be greater than this.
    float jet_pt_max,                           // Maximum jet pT to accept. At least 1 jet per event will be less than this.
    float jet_eta_max = d_jet_eta_max,          // Largest jet rapidity to consider.
    float pt_hat_min = d_pt_hat_min,            // [GeV] ptHatMin value.
    float pt_hat_max = d_pt_hat_max,            // [GeV] ptHatMax value. If 0. then no upper limit.
    float fastjet_radius = d_fastjet_radius,    // Jet radius used by FastJet.
    float fastjet_pt_min = d_fastjet_pt_min,    // [GeV] Minimum pT considered by FastJet for a jet.
    float fastjet_match_radius = d_fastjet_match_radius,    // Matches combined jets to PYTHIA jets within this radius.
    int   softest_jet_index = d_softest_jet_index,          // If =0, tries to match with all PYTHIA jets. If >0, matches only to top N PYTHIA jets.
    int   accept_jet_index = d_accept_jet_index,            // If =0, accepts events with any jet between min/max. If >0, only accepts events with one of top N jets between min/max.
    int   thermal_mean = d_thermal_mean,        // Mean number of thermal particles
    int   thermal_sigma = d_thermal_sigma,      // Standard deviation for number of thermal particles
    float beam_power = d_beam_power,            // [GeV] Simulated beam power in the dector. Default defined by jet_ml_constants.h file.
    float detector_eta = d_detector_eta         // Maximum rapidity of the detector. Default defined by jet_ml_constants.h file.
    ) {
    
    std::cout << ">>> Generating PYTHIA Events <<<" << std::endl;
    
    char pythia_file_path[500];
    snprintf(pythia_file_path, 500, "%s", Pythia_Generator(
        output_base_name,
        output_directory,
        event_count,
        pt_bias_power,
        jet_pt_min,
        jet_pt_max,
        jet_eta_max,
        pt_hat_min,
        pt_hat_max,
        accept_jet_index,
        fastjet_radius,
        fastjet_pt_min,
        beam_power,
        detector_eta
        ));

    std::cout << ">>> Generating Thermal Events <<<" << std::endl;
    char thermal_file_path[500];
    snprintf(thermal_file_path, 500, "%s", Thermal_Generator(
        output_base_name,
        output_directory,
        event_count,
        pt_bias_power,
        jet_pt_min,
        jet_pt_max,
        thermal_mean,
        thermal_sigma,
        beam_power,
        detector_eta
        ));
    
    std::cout << ">>> Combining PYTHIA and Thermal Events <<<" << std::endl;
    char combined_file_path[500];
    snprintf(combined_file_path, 500, "%s", Combine_Events(
        output_base_name,
        output_directory,
        pythia_file_path,
        thermal_file_path
        ));
    
    std::cout << ">>> Clustering Combined Jets <<<" << std::endl;
    Combined_Jet_Clusterer(
        combined_file_path
        );
    
    std::cout << ">>> Preparing ML File <<<" << std::endl;
    Jet_ML_Prep(
        output_base_name,
        output_directory,
        combined_file_path,
        pythia_file_path,
        jet_pt_min,
        jet_pt_max,
        jet_eta_max,
        softest_jet_index,
        fastjet_match_radius
        );
        
    std::cout << ">>> Event Generation and Jet Clustering Complete! <<<" << std::endl;
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
    
    // Example using minimum number of parameters and default values . Note beam_power and detector_eta have defai
    Jet_Generator::Jet_Generator_Full(
        "Train",        // output_base_name (pt_bias, pt_min, pt_max, and event_cout will be added to this name in the functions)
        dir_data,       // output_directory
        1000,           // event_count (Generates 100,000 events but will produce MORE total jets)
        0.,             // pt_bias_power (Here NO BIAS is applied)
        10.,            // jet_pt_min (Events will have at least one jet above
        90.             // jet_pt_max
        );
    
    // Example defining all parameters.
    // NOTE: Must define ALL parameters up to the last one used!
    // (e.g. If defining fastjet_radius then jet_eta_max, pt_hat_min, pt_hat_max, must also be defined - it's not like Python!)
    Jet_Generator::Jet_Generator_Full(
        "Train",        // output_base_name (Typically just "Train" or "Test")
        dir_data,       // output_directory
        1000,           // event_count [int]
        0.,             // pt_bias_power (Generates events with bias of pT^N for event energies distribution)
        10.,            // jet_pt_min [GeV]
        90.,            // jet_pt_max [GeV]
        0.5,            // jet_eta_max (Reasonable value is detector_eta - fastjet_radius)
        7.5,            // pt_hat_min (Reasonable value is 0.75 * jet_pt_min)
        0.,             // pt_hat_max (If =0. then no upper limit)
        0.4,            // fastjet_radius
        5.,             // fastjet_pt_min (Minimum pT considered by FastJet for a jet)
        0.3,            // fastjet_match_radius (Reasonable value is less than fastjet_radius)
        0,              // softest_jet_index [int] (If =0, tries to match with all PYTHIA jets. If >0, matches only to top N PYTHIA jets)
        0,              // accept_jet_index [int] (If =0, accepts events with any jet in min/max. If >0, only accepts events with one of top N jets in min/max.)
        1800,           // thermal_mean
        200,            // thermal_sigma
        2760.,          // beam_power [GeV] (Simulated beam power in the detector)
        0.9             // detector_eta (Maximum rapidity of the detector)
        );
    
    return 0;
}
