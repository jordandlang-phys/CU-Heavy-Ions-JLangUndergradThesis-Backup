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
#include "TMath.h"
#include "TRandom.h"
#include <cmath>
#include <iostream>
#include <filesystem>
#include "_header.h"
using namespace Project_Constants;

const bool print_out   = true;
const int  print_every_x = 100;
const bool debug       = false;
const bool use_voronoi = false; // Uses Voronoi for jet clustering



// ----- PYTHIA GENERATOR -----


// Add input for file name base/root, file path, etc.
// Make sure that nothing in functions relies on the header file!
// Remove jet energy (we don't use it anywhere)
void Pythia_Generator(
    char* file_name,
    int func_eventCount,
    float func_beamPower,
    float func_ptHatMin,
    float func_ptBiasPow,
    float func_slimMin,
    float func_slimMax,
    float func_slimRap) {
    
    char file_path[200];
    int temp_int = sprintf(file_path, "%s/Pyth_%s", dir_data, file_name);
    TFile* output_file = new TFile(file_path, "UPDATE");
    TTree* pythia_tree = new TTree("Pythia_Tree","Tree of particle jet events from PYTHIA p+p collisions");
    TTree* jet_tree    = new TTree("FastJet_Tree","Tree of jet clusters by event from PYTHIA");
    
    if (print_out) std::cout << "Output file and trees successfully generated." << std::endl;
    
    // Output Variables
    
    // Initializes arrays for particle info, DO NOT USE DOUBLES!!!
    const int max_parts  = 200;
    const int max_jets   = 100;
    const int max_consts = 200;
    
    int     particle_n;                 // particle index
    float   particle_pt[max_parts];     // transverse momentum
    float   particle_phi[max_parts];    // azimuthal angle
    float   particle_eta[max_parts];    // pseudorapidity
    int     particle_pid[max_parts];    // numerical code for particle type
    int     particle_m[max_parts];      // particle mass
    
    int     jet_n;                      // jet index
    float   jet_pt[max_jets];           // transverse momentum
    float   jet_y[max_jets];            // rapidity (FastJet does not output eta)
    float   jet_phi[max_jets];          // azimuthal angle
    float   jet_E[max_jets];            // jet energy
    float   jet_mass[max_jets];
    float   jet_area[max_jets];
    float   jet_area_err[max_jets];
    int     jet_const_n[max_jets];
    float   jet_const_pt[max_jets][max_consts];
    float   jet_const_eta[max_jets][max_consts];
    float   jet_const_phi[max_jets][max_consts];
    float   jet_const_E[max_jets][max_consts];
    
    // Builds particle and jet branches
    // Note: Can only do variable size for first variable in array
    pythia_tree->Branch("particle_n",   &particle_n);
    pythia_tree->Branch("particle_pt",  particle_pt,  "particle_pt[particle_n]/F");
    pythia_tree->Branch("particle_eta", particle_eta, "particle_eta[particle_n]/F");
    pythia_tree->Branch("particle_phi", particle_phi, "particle_phi[particle_n]/F");
    pythia_tree->Branch("particle_pid", particle_pid, "particle_pid[particle_n]/I");
    pythia_tree->Branch("particle_m",   particle_m,   "particle_m[particle_n]/I");
    
    jet_tree->Branch("jet_n",        &jet_n);
    jet_tree->Branch("jet_pt",       jet_pt,       "jet_pt[jet_n]/F");
    jet_tree->Branch("jet_y",        jet_y,        "jet_y[jet_n]/F");
    jet_tree->Branch("jet_phi",      jet_phi,      "jet_phi[jet_n]/F");
    jet_tree->Branch("jet_E",        jet_E,        "jet_E[jet_n]/F");
    jet_tree->Branch("jet_mass",     jet_mass,     "jet_mass[jet_n]/F");
    jet_tree->Branch("jet_area",     jet_area,     "jet_area[jet_n]/F");
    jet_tree->Branch("jet_area_err", jet_area_err, "jet_area_err[jet_n]/F");
    jet_tree->Branch("jet_const_n",     jet_const_n,    "jet_const_n[jet_n]/I");
    jet_tree->Branch("jet_const_pt",    jet_const_pt,   "jet_const_pt[jet_n][400]/F");
    jet_tree->Branch("jet_const_eta",   jet_const_eta,  "jet_const_eta[jet_n][400]/F");
    jet_tree->Branch("jet_const_phi",   jet_const_phi,  "jet_const_phi[jet_n][400]/F");
    jet_tree->Branch("jet_const_E",     jet_const_E,    "jet_const_E[jet_n][400]/F");
    
    // --- PYTHIA Setup ---
    // --- LHC process and output selection. Initialization.
    Pythia pythia;
    Settings& pythia_settings = pythia.settings;
    pythia_settings.parm("Beams:eCM", func_beamPower);  // Uses units of GeV
    pythia_settings.parm("PhaseSpace:pTHatMin", func_ptHatMin );
    pythia_settings.parm("PhaseSpace:pTHatMax", 0.);
    pythia.readString("HardQCD:all = on");              // Turns on hard scattering
    if ( func_ptBiasPow >= 0. ) {
        pythia_settings.parm("PhaseSpace:bias2SelectionPow", func_ptBiasPow);
        pythia.readString("PhaseSpace:bias2Selection = on");
        pythia.readString("PhaseSpace:bias2SelectionRef = 100.");
    }
    pythia.init();
    // -- to try later
    // multiple ptHatMin
    // combine and reweight ptHatMins
    
    if ( print_out ) std::cout << "PYTHIA settings initialized." << std::endl;
    
    // --- FastJet Setup ---
    // --- Cluster each PYTHIA event
    fastjet::JetDefinition jet_definition(fastjet::antikt_algorithm, fj_jetR);
    fastjet::AreaDefinition area_definition;
    
    if (!use_voronoi) {
        double ghost_etamax = 0.9;
        double ghost_area    = 0.05;
        int    active_area_repeats = 5;
        fastjet::GhostedAreaSpec ghost_spec(ghost_etamax, active_area_repeats, ghost_area);
        area_definition = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
    } else {
        double effective_Rfact = 1.0;
        area_definition = fastjet::VoronoiAreaSpec(effective_Rfact);
    }
    
    if ( print_out ) std::cout << "FastJet settings initialized." << std::endl;
    
    // Loop to generate PYTHIA events. Skip if error.
    int e = 0;
    while ( e < func_eventCount ) {
        // ---- PYTHIA -----
        if (!pythia.next()) continue;

        particle_n = 0;

        for ( int p = 0; p < pythia.event.size(); p++ ) {
            if (!pythia.event[p].isFinal()) continue;           // Ends if final particle
            if (!pythia.event[p].isCharged()) continue;         // Skips neutral particles
            if (pythia.event[p].pT() < 0.15) continue;          // Skips particles with pT < 0.15
            if (fabs(pythia.event[p].eta()) > 0.9) continue;    // Skips particles with eta outside +/-0.9

            particle_pt[particle_n]  = pythia.event[p].pT();
            particle_eta[particle_n] = pythia.event[p].eta();
            particle_phi[particle_n] = pythia.event[p].phi();
            particle_pid[particle_n] = pythia.event[p].id();
            particle_m[particle_n] = -1 * pythia.event[p].m(); // Outputs mass as negative?

            particle_n++;
        }
        
        // ---------- FASTJET ----------
        
        TLorentzVector input_vec;
        std::vector<fastjet::PseudoJet> input_particles;
        
        for ( int p = 0 ; p < particle_n ; p++) {
            input_vec.SetPtEtaPhiM(particle_pt[p], particle_eta[p], particle_phi[p], particle_m[p]);
            input_particles.push_back(fastjet::PseudoJet( input_vec.Px(), input_vec.Py(), input_vec.Pz(), input_vec.E() ));
        }
        
        fastjet::ClusterSequenceArea jet_clusters(input_particles, jet_definition, area_definition);
        std::vector<fastjet::PseudoJet> jet_list = sorted_by_pt(jet_clusters.inclusive_jets(fj_jetptmin));
        // jet_clusters.inclusive_jets(fj_jetptmin)  // Use if NOT recording area
        
        jet_n = jet_list.size();
        if ( jet_n == 0 ) continue;
        if ( jet_n >= max_jets ) jet_n = max_jets;
        
        bool acceptJet = false;
        
        for ( int j = 0 ; j < jet_n ; j++ ) {
            jet_pt[j]   = jet_list[j].pt();
            jet_y[j]    = jet_list[j].rap();
            jet_phi[j]  = jet_list[j].phi();
            jet_E[j]    = jet_list[j].E();
            jet_mass[j] = jet_list[j].m();
            jet_area[j] = jet_list[j].area();
            jet_area_err[j] = jet_list[j].area_error();
            
            if ( (jet_pt[j] >= func_slimMin) && (jet_pt[j] <= func_slimMax) &&
                (jet_y[j] >= -func_slimRap) && (jet_y[j] <= func_slimRap) ) acceptJet = true;
            
            std::vector<fastjet::PseudoJet> jet_constituents = jet_list[j].constituents();
            jet_const_n[j]  = jet_constituents.size();
            
            for ( int p = 0 ; p < jet_constituents.size() ; p++) {
                jet_const_pt[j][p]  = jet_constituents[p].pt();
                jet_const_eta[j][p] = jet_constituents[p].eta();
                jet_const_phi[j][p] = jet_constituents[p].phi();
                jet_const_E[j][p]   = jet_constituents[p].E();
            }
        }
        
        if ( acceptJet ) {
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
}



// ----- THERMAL GENERATOR -----



TF1* n_func;
TF1* pt_func;
TF1* eta_func;
TF1* phi_func;

double Gaussian_Func(double *x,double *par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    double fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

double Power_Func(double *x,double *par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[2]);
    double fitval = par[0] * pow(arg, par[1]);
    return fitval;
}

double Flat_Func(double *x,double *par) {
    double fitval = par[0] * 1;
    return fitval;
}

double ModifiedHagedorn_Func(double* x, double* par) {
    double arg1 = par[0] * pow(x[0], 2) / pow( (pow(x[0], 2) + pow(par[1], 2)), 0.5);
    double arg2 = pow( 1 + (x[0] /par[2]), par[3] );
    double arg3 = 1 / x[0];
    double fitval =  arg1 * arg2; // * arg3;
    return fitval;
}

int GenN() {
    int n = round( n_func->GetRandom(0, 2 * gaus_mean) );
    return n;
}

double GenPt() {
//    double pt = pt_func->GetRandom(0, 100);
    double pt = pt_func->GetRandom(0, 100);
    return pt;
}

double GenEta() {
    double eta = eta_func->GetRandom(-1 * eta_max, eta_max);
    return eta;
}

double GenPhi() {
    double phi = phi_func->GetRandom(-1 * math_pi, math_pi);
    return phi;
}

void Thermal_Generator(char* file_name, int func_eventCount) {
    
    char file_path[200];
    int temp_int = sprintf(file_path, "%s/Ther_%s", dir_data, file_name);
    TFile* output_file = new TFile(file_path, "UPDATE");
    
    if (print_out) std::cout << "File Created: " << output_file << std::endl;
    TTree* output_tree = new TTree("Thermal_Tree","Tree of thermal particles as background for events");
    if (print_out) std::cout << "TTree Created: Thermal_Tree" << std::endl;
    
    n_func = new TF1("n_func", Gaussian_Func, 0, 2 * gaus_mean, 3);
    pt_func = new TF1("pt_func", ModifiedHagedorn_Func, 0, 100, 4);
    eta_func = new TF1("eta_func", Flat_Func, -1 * eta_max, eta_max, 1);
    phi_func = new TF1("phi_func", Flat_Func, -1 * math_pi, math_pi, 1);

    n_func->SetParameters(gaus_norm, gaus_mean, gaus_sigma);
    pt_func->SetParameters(64547, 3.076, 1.126, -8.491);
    eta_func->SetParameter(0, 1.);
    phi_func->SetParameter(0, 1.);
    
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
    
    for (int e = 0 ; e < func_eventCount ; e++) {
        // Initializes particle_n to 0 for the event
        particle_n = 0;
        
        // Generates random number of particles for thermal background
        int nParticles = GenN();
        
        if ( print_out && (e % print_every_x) == 0 ) {
            std::cout << "----- Printing Event #" << e << " -----" << std::endl;
        }
        else if ( (e % print_every_x) == 0 ) {
            std::cout << "Generating Event #" << e << "..." << std::endl;
        }
        
        // For each particle, assigns pt, eta, phi, then increments n
        for (int p = 0 ; p < nParticles ; p++) {
            particle_pt[particle_n] = GenPt();
            particle_eta[particle_n] = GenEta();
            particle_phi[particle_n] = GenPhi();
            
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
}

//for (int e = 0 ; e < tree->GetEntries ; e++ ) { tree->GetEntry(e); for (int j = 0 ; j < jet_n ; j++) { temp_hist->Fill(jet_pt[j]) }}

// ----- COMBINE EVENTS -----



void Combine_Events(char* file_name) {
    
    // Open PYTHIA events data and access thermal TTree
    char pythia_file_path[200];
    sprintf(pythia_file_path, "%s/Pyth_%s", dir_data, file_name);
    TFile* pythia_file = new TFile(pythia_file_path, "READ");
    std::cout << "Reading PYTHIA File" << std::endl;
    
    TTree* pythia_tree = (TTree*) pythia_file->Get("Pythia_Tree");
    
    int   pythia_n;
    float pythia_pt[200];
    float pythia_eta[200];
    float pythia_phi[200];
    
    pythia_tree->SetBranchAddress("particle_n",   &pythia_n);
    pythia_tree->SetBranchAddress("particle_pt",  pythia_pt);
    pythia_tree->SetBranchAddress("particle_eta", pythia_eta);
    pythia_tree->SetBranchAddress("particle_phi", pythia_phi);
    
    if (print_out) std::cout << "TTree Accessed: Pythia_Tree" << std::endl;
    
    // Open thermal events data and access thermal TTree
    char thermal_file_path[200];
    sprintf(thermal_file_path, "%s/Ther_%s", dir_data, file_name);
    TFile* thermal_file = new TFile(thermal_file_path, "READ");
    std::cout << "Reading Thermal File" << std::endl;
    
    TTree* thermal_tree = (TTree*) thermal_file->Get("Thermal_Tree");
    
    int   thermal_n;
    float thermal_pt[4000];
    float thermal_eta[4000];
    float thermal_phi[4000];
    
    thermal_tree->SetBranchAddress("particle_n",   &thermal_n);
    thermal_tree->SetBranchAddress("particle_pt",  thermal_pt);
    thermal_tree->SetBranchAddress("particle_eta", thermal_eta);
    thermal_tree->SetBranchAddress("particle_phi", thermal_phi);
    
    if (print_out) std::cout << "TTree Accessed: Thermal_Tree" << std::endl;
    
    // Create new combined events data file
    char combined_file_path[200];
    sprintf(combined_file_path, "%s/Comb_%s", dir_data, file_name);
    TFile *combined_file = new TFile(combined_file_path, "UPDATE");
    std::cout << "Combined File Created" << std::endl;
    
    // Create new combined events TTree
    TTree *combined_tree = new TTree("Combined_Tree","Tree containing jet and thermal background particles.");
    
    if (print_out) std::cout << "TTree Created: Combined_Tree" << std::endl;
    
    // Build out combined TTree
    int   combined_n;
    float combined_pt[4200];
    float combined_eta[4200];
    float combined_phi[4200];
    int   combined_jet_n;
    int   combined_jet_class;
    
    combined_tree->Branch("particle_n",     &combined_n);
    combined_tree->Branch("particle_pt",    combined_pt,    "particle_pt[particle_n]/F");
    combined_tree->Branch("particle_eta",   combined_eta,   "particle_eta[particle_n]/F");
    combined_tree->Branch("particle_phi",   combined_phi,   "particle_phi[particle_n]/F");
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
        }
        
        // Loop to add thermal particles
        for (int t = 0 ; t < thermal_count ; t++) {
            combined_n += 1;
            combined_pt[pythia_count + t - 1] = thermal_pt[t];
            combined_eta[pythia_count + t - 1] = thermal_eta[t];
            combined_phi[pythia_count + t - 1] = thermal_phi[t];
            combined_jet_class = 0;
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
    combined_file->Close();
    pythia_file->Close();
    thermal_file->Close();
    
    if ( print_out ) std::cout << "File written to and closed." << std::endl;
}



// ----- COMBINED EVENT JET CLUSTERER -----



void Jet_Clusterer(char* file_name) {
    
    char* input_tree_name = "Combined_Tree";
    
    // Defines the jet definition and jet area definition
    fastjet::JetDefinition jet_definition(fastjet::antikt_algorithm, fj_jetR);
    fastjet::AreaDefinition area_definition;
    
    bool use_voronoi = false;
    if (!use_voronoi) {
        double ghost_etamax = 0.9;
        double ghost_area    = 0.05;
        int    active_area_repeats = 5;
        fastjet::GhostedAreaSpec ghost_spec(ghost_etamax, active_area_repeats, ghost_area);
        area_definition = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
    } else {
        double effective_Rfact = 1.0;
        area_definition = fastjet::VoronoiAreaSpec(effective_Rfact);
    }
    
    if ( print_out ) std::cout << "FastJet settings initialized." << std::endl;
        
    // Input Variables and Trees
    int   input_ntotal = 0;
    int   input_n;
    float input_pt[4200];
    float input_eta[4200];
    float input_phi[4200];
    TLorentzVector input_vec;

    char input_file_path[200];
    int temp_int_1 = sprintf(input_file_path, "%s/Comb_%s", dir_data, file_name);
    TFile* input_file = new TFile(input_file_path, "UPDATE");
    TTree* input_tree = (TTree*) input_file->Get(input_tree_name);
    TTree* output_tree = new TTree("FastJet_Tree","Tree of jet clusters");

    input_tree->SetBranchAddress("particle_n", &input_n);
    input_tree->SetBranchAddress("particle_pt", input_pt);
    input_tree->SetBranchAddress("particle_eta", input_eta);
    input_tree->SetBranchAddress("particle_phi", input_phi);
    
    if (print_out) std::cout << "Input file and tree successfully accessed." << std::endl;
    
//    // Sets output file
//    char output_file_path[200];
//    int temp_int_2 = sprintf(output_file_path, "%s/Comb_%s", dir_data, file_name);
//    TFile* output_file = new TFile(output_file_path, "UPDATE");
//    TTree* output_tree = new TTree("FastJet_Tree","Tree of jet clusters");
    
    if (print_out) std::cout << "Output file and tree successfully generated." << std::endl;
    
    // Output Variables
    const int max_jets = 100;
    const int max_parts = 400;
    
    int     jet_n;
    float   jet_pt[max_jets];
    float   jet_y[max_jets];
    float   jet_phi[max_jets];
    float   jet_E[max_jets];
    float   jet_mass[max_jets];
    float   jet_area[max_jets];
    float   jet_area_err[max_jets];
    int     jet_const_n[max_jets];
    float   jet_const_pt[max_jets][max_parts];
    float   jet_const_eta[max_jets][max_parts];
    float   jet_const_phi[max_jets][max_parts];
    float   jet_const_E[max_jets][max_parts];
    
    // Note: Can only do variable size for first variable in array
    output_tree->Branch("jet_n",        &jet_n);
    output_tree->Branch("jet_pt",       jet_pt,       "jet_pt[jet_n]/F");
    output_tree->Branch("jet_y",        jet_y,        "jet_y[jet_n]/F");
    output_tree->Branch("jet_phi",      jet_phi,      "jet_phi[jet_n]/F");
    output_tree->Branch("jet_E",        jet_E,        "jet_E[jet_n]/F");
    output_tree->Branch("jet_mass",     jet_mass,     "jet_mass[jet_n]/F");
    output_tree->Branch("jet_area",     jet_area,     "jet_area[jet_n]/F");
    output_tree->Branch("jet_area_err", jet_area_err, "jet_area_err[jet_n]/F");
    output_tree->Branch("jet_const_n",      jet_const_n,    "jet_const_n[jet_n]/I");
    output_tree->Branch("jet_const_pt",     jet_const_pt,   "jet_const_pt[jet_n][400]/F");
    output_tree->Branch("jet_const_eta",    jet_const_eta,  "jet_const_eta[jet_n][400]/F");
    output_tree->Branch("jet_const_phi",    jet_const_phi,  "jet_const_phi[jet_n][400]/F");
    output_tree->Branch("jet_const_E",      jet_const_E,    "jet_const_E[jet_n][400]/F");
    
    if ( print_out ) std::cout << "Output tree variables initialized." << std::endl;
    
    // Loop to fill vector with particles
    for ( int e = 0 ; e < input_tree->GetEntries() ; e++ ) {
        
        if ( (e % print_every_x) == 0 ) std::cout << "Processing Event " << e << "." << std::endl;
        
        input_tree->GetEntry(e);
        std::vector<fastjet::PseudoJet> input_particles;
        
        for ( int p = 0 ; p < input_n ; p++) {
            input_vec.SetPtEtaPhiM(input_pt[p], input_eta[p], input_phi[p], m_pion);
            input_particles.push_back(fastjet::PseudoJet( input_vec.Px(), input_vec.Py(), input_vec.Pz(), input_vec.E() ));
        }
        
        fastjet::ClusterSequenceArea jet_clusters(input_particles, jet_definition, area_definition);
        std::vector<fastjet::PseudoJet> jet_list = sorted_by_pt(jet_clusters.inclusive_jets(fj_jetptmin));
        // jet_clusters.inclusive_jets(fj_jetptmin)  // Use if NOT recording area
        
        jet_n = jet_list.size();
        if ( jet_n >= max_jets ) continue;
        if ( print_out && ( e % print_every_x ) == 0 ) std::cout << "Event #" << e << " has " << jet_n << " jets." << std::endl;
        
        if ( jet_n == 0 ) continue;
        
        for ( int j = 0 ; j < jet_n ; j++ ) {
            jet_pt[j]   = jet_list[j].pt();
            jet_y[j]    = jet_list[j].rap();
            jet_phi[j]  = jet_list[j].phi();
            jet_E[j]    = jet_list[j].E();
            jet_mass[j] = jet_list[j].m();
            jet_area[j] = jet_list[j].area();
            jet_area_err[j] = jet_list[j].area_error();
            
            std::vector<fastjet::PseudoJet> jet_constituents = jet_list[j].constituents();
            jet_const_n[j]  = jet_constituents.size();
                        
            for ( int p = 0 ; p < jet_const_n[j] ; p++) {
                jet_const_pt[j][p]   = jet_constituents[p].pt();
                jet_const_eta[j][p]  = jet_constituents[p].eta();
                jet_const_phi[j][p]  = jet_constituents[p].phi();
                jet_const_E[j][p]    = jet_constituents[p].E();
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
//    output_file->Close();
    
    delete input_file;
//    delete output_file;
        
    std::cout << "Files closed." << std::endl;
}



// ----- EVENT GENERATOR -----



void Event_Generator(char* file_name, int func_eventCount, float func_beamPower,
                     float func_ptBiasPow, float func_slimMin, float func_slimMax, float func_slimRap) {
    
    float func_ptHatMin = ptHatMin_ptMin_ratio * func_slimMin;
    
    std::cout << ">>> Generate PYTHIA Events <<<" << std::endl;
    Pythia_Generator(file_name, func_eventCount, func_beamPower, func_ptHatMin,
                     func_ptBiasPow, func_slimMin, func_slimMax, func_slimRap);

    std::cout << ">>> Generate Thermal Events <<<" << std::endl;
    Thermal_Generator(file_name, func_eventCount);

    std::cout << ">>> Combine PYTHIA and Thermal Events <<<" << std::endl;
    Combine_Events(file_name);
    
    std::cout << ">>> Cluster Combined Jets <<<" << std::endl;
    Jet_Clusterer(file_name);
    
    std::cout << ">>> Event Generation and Jet Clustering Complete! <<<" << std::endl;
}



// ----- MAIN -----



int main() {
    bool print_out = ::print_out;
    
    std::__fs::filesystem::create_directories(dir_master);
    std::__fs::filesystem::create_directories(dir_data);
    std::__fs::filesystem::create_directories(dir_plots);
    
    Event_Generator(
        "10_90_B8_Test_Trees.root", // file name
        500000,     // number of events to generate
        beamPower,  // beam power
        8.,        // pt bias power (pt^x), set to -1. to disable bias
        10.,        // pt min for slimming
        90.,        // pt max for slimming
        slim_rap);  // max rapidity for slimming
    
//    Event_Generator(
//        "10_90_B0_Train_Trees.root", // file name
//        500000,     // number of events to generate
//        beamPower,  // beam power
//        -1.,        // pt bias power (pt^x), set to -1. to disable bias
//        10.,        // pt min for slimming
//        90.,        // pt max for slimming
//        slim_rap);  // max rapidity for slimming
    
//    Event_Generator(
//        "10_90_B4_Train_Trees.root", // file name
//        500000,     // number of events to generate
//        beamPower,  // beam power
//        4.,        // pt bias power (pt^x), set to -1. to disable bias
//        10.,        // pt min for slimming
//        90.,        // pt max for slimming
//        slim_rap);  // max rapidity for slimming
//
//    Event_Generator(
//        "10_90_B8_Train_Trees.root", // file name
//        500000,     // number of events to generate
//        beamPower,  // beam power
//        8.,        // pt bias power (pt^x), set to -1. to disable bias
//        10.,        // pt min for slimming
//        90.,        // pt max for slimming
//        slim_rap);  // max rapidity for slimming
//
//    Event_Generator(
//        "40_60_B0_Test_Trees.root", // file name
//        100000,     // number of events to generate
//        beamPower,  // beam power
//        -1,        // pt bias power (pt^x), set to -1. to disable bias
//        40.,        // pt min for slimming
//        60.,        // pt max for slimming
//        slim_rap);  // max rapidity for slimming
//
//    Event_Generator(
//        "40_60_B4_Test_Trees.root", // file name
//        100000,     // number of events to generate
//        beamPower,  // beam power
//        4,        // pt bias power (pt^x), set to -1. to disable bias
//        40.,        // pt min for slimming
//        60.,        // pt max for slimming
//        slim_rap);  // max rapidity for slimming
//
//    Event_Generator(
//        "40_60_B8_Test_Trees.root", // file name
//        100000,     // number of events to generate
//        beamPower,  // beam power
//        8,        // pt bias power (pt^x), set to -1. to disable bias
//        40.,        // pt min for slimming
//        60.,        // pt max for slimming
//        slim_rap);  // max rapidity for slimming
    
    return 0;
}
