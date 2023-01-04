#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TGraph.h"

#include "_header.h"
using namespace Project_Constants;

const bool printout = true;
const bool debug = true;

void MinMaxBinning(float* array, float lowerBound, float upperBound) {
    int arrayLength = sizeof(array)/sizeof(array[0]);
    float binRange = upperBound - lowerBound;
    float binSize = binRange / arrayLength;
    for (int i = 0; i < arrayLength; i++) {
        array[i] = i * binSize + lowerBound;
    }
}





void Particle_Analyzer( char* input_file_name, char* input_tree_name, int label_arr_size, char* label_arr[10], bool bool_make_plots ) {
    // Sets file prefix based on input tree name
    char particle_dir[200];
    sprintf(particle_dir, "%s/Plots_Particles", dir_plots);
    std::__fs::filesystem::create_directories(particle_dir);
    
    char file_prefix[4];
    char file_title[100];
    if ( input_tree_name == "Pythia" ) {
        sprintf(file_prefix, "Pyth");
        sprintf(file_title, "PYTHIA Particles");
    }
    else if ( input_tree_name == "Thermal" ) {
        sprintf(file_prefix, "Ther");
        sprintf(file_title, "Thermal Particles");
    }
    else if ( input_tree_name == "Combined" ) {
        sprintf(file_prefix, "Comb");
        sprintf(file_title, "PYTHIA + Thermal Particles");
    }
    
    // Opens and reads the Root file
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s_%s_Trees.root", dir_data, file_prefix, input_file_name);
    TFile* input_file = new TFile(input_file_path, "READ");
    
    char input_tree_path[200];
    sprintf(input_tree_path, "%s_Tree", input_tree_name);
    TTree* input_tree = (TTree*) input_file->Get(input_tree_path);
    
    // Opens and makes the output file
    char output_file_path[200];
    sprintf(output_file_path, "%s/%s_%s_Plots.root", dir_data, file_prefix, input_file_name);
    TFile* output_file = new TFile(output_file_path, "UPDATE");
    
    std::cout << "Input file read. TTree accessed." << std::endl;
    
    // Variables
    int     particle_ntotal = 0;
    int     particle_n;
    float   particle_pt[4200];
    float   particle_eta[4200];
    float   particle_phi[4200];
    
    input_tree->SetBranchAddress("particle_n",      &particle_n);
    input_tree->SetBranchAddress("particle_pt",     particle_pt);
    input_tree->SetBranchAddress("particle_eta",    particle_eta);
    input_tree->SetBranchAddress("particle_phi",    particle_phi);
    
    int event_n = input_tree->GetEntries();
    
    // Declare histograms
    // ("Name in Root file", "Plot Title; X-Axis; Y-Axis", #Bins, Left edge of lowest bin (included), Right edge of highest bin (not included) )
    
    TH1D* th1d_particle_n;
    TH1D* th1d_particle_pt;
    TH1D* th1d_particle_eta;
    TH1D* th1d_particle_phi;
    
    char plot_id[200];
    char plot_title[200];
    
    if ( input_tree_name == "Pythia" ) {
        sprintf(plot_id, "th1d_%s_particle_n", input_tree_name);
        sprintf(plot_title, "%s, Particle Count; N_{ch}; N_{event}", file_title);
        th1d_particle_n = new TH1D(plot_id, plot_title, 31, -0.5, 92.5);
        
        sprintf(plot_id, "th1d_%s_particle_pt", input_tree_name);
        sprintf(plot_title, "%s, Particle p_{T}; p_{T, ch} [GeV/c]; dN_{ch}/dp_{T}", file_title);
        th1d_particle_pt = new TH1D(plot_id, plot_title, 20, 0, 80);
    }
    
    else {
        sprintf(plot_id, "th1d_%s_particle_n", input_tree_name);
        sprintf(plot_title, "%s, Particle Count; N_{ch}; N_{event}", file_title);
        th1d_particle_n = new TH1D(plot_id, plot_title, 50, 800, 2800);
        
        sprintf(plot_id, "th1d_%s_particle_pt", input_tree_name);
        sprintf(plot_title, "%s, Particle p_{T}; p_{T, ch} [GeV/c]; dN_{ch}/dp_{T}", file_title);
        th1d_particle_pt = new TH1D(plot_id, plot_title, 40, 0, 80);
    }
    
    sprintf(plot_id, "th1d_%s_particle_eta", input_tree_name);
    sprintf(plot_title, "%s, Particle Pseudorapidity; #eta_{ch}; dN_{ch}/d#eta", file_title);
    th1d_particle_eta = new TH1D(plot_id, plot_title, 20, -0.9, 0.9);
    
    sprintf(plot_id, "th1d_%s_particle_phi", input_tree_name);
    sprintf(plot_title, "%s, Particle Angle; #phi_{ch} [rad]; dN_{ch}/d#phi", file_title);
    th1d_particle_phi = new TH1D(plot_id, plot_title, 20, -math_pi, math_pi);
    
    // Automatically calculate error
    th1d_particle_n     ->Sumw2();
    th1d_particle_pt    ->Sumw2();
    th1d_particle_eta   ->Sumw2();
    th1d_particle_phi   ->Sumw2();
    
    // Iterates through events
    for (int e = 0; e < event_n; e++) {
        input_tree->GetEntry(e);
        
        // Fills tree with numbers of particles and jets per event
        th1d_particle_n->Fill(particle_n);
        particle_ntotal += particle_n;
        
        // Fills tree with particle info
        for (int p = 0; p < particle_n; p++) {
            // Fills particle pt, eta, phi
            th1d_particle_pt    ->Fill(particle_pt[p]);
            th1d_particle_eta   ->Fill(particle_eta[p]);
            th1d_particle_phi   ->Fill(particle_phi[p]);
        }
    }
    
    // Normalizes histograms
    th1d_particle_pt    ->Scale(1, "width");
    th1d_particle_eta   ->Scale(1, "width");
    th1d_particle_phi   ->Scale(1, "width");
    
    // Writes histograms to output file
    th1d_particle_n     ->Write("", TObject::kOverwrite);
    th1d_particle_pt    ->Write("", TObject::kOverwrite);
    th1d_particle_eta   ->Write("", TObject::kOverwrite);
    th1d_particle_phi   ->Write("", TObject::kOverwrite);
    
    output_file ->Write();
    
    // --- PLOT PARTICLE HISTOGRAMS ---
    
    if ( bool_make_plots ) {
        
        char plot_name[200];
        
        // Make hist_particle_n
        sprintf(plot_name, "Plots_Particles/th1d_%s_particle_n.pdf", input_tree_name);
        th1d_plotter(
            th1d_particle_n,
            plot_name,
            particle_line_color, particle_marker_color, particle_marker_style, 1.,
            0, false,
            label_arr, label_arr_size, 0.65, 0.75);
        
        // Make hist_particle_pt
        sprintf(plot_name, "Plots_Particles/th1d_%s_particle_pt.pdf", input_tree_name);
        th1d_plotter(
            th1d_particle_pt,
            plot_name,
            particle_line_color, particle_marker_color, particle_marker_style, 1.,
            .1, true,
            label_arr, label_arr_size, 0.65, 0.75);
        
        // Make hist_particle_eta
        sprintf(plot_name, "Plots_Particles/th1d_%s_particle_eta.pdf", input_tree_name);
        th1d_plotter(
            th1d_particle_eta,
            plot_name,
            particle_line_color, particle_marker_color, particle_marker_style, 1.,
            0, false,
            label_arr, label_arr_size, 0.65, 0.25);
        
        // Make hist_particle_phi
        sprintf(plot_name, "Plots_Particles/th1d_%s_particle_phi.pdf", input_tree_name);
        th1d_plotter(
            th1d_particle_phi,
            plot_name,
            particle_line_color, particle_marker_color, particle_marker_style, 1.,
            0, false,
            label_arr, label_arr_size, 0.65, 0.25);
    }

    // Closes files
    output_file->Close();
    input_file->Close();
    
    delete output_file;
    delete input_file;
}





void Particle_PID_Analyzer( char* input_file_name, char* input_tree_name, int label_arr_size, char* label_arr[10], bool bool_make_plots ) {
    
    // Opens and reads the Root file
    char input_file_path[200];
    sprintf(input_file_path, "%s/Pyth_%s_Trees.root", dir_data, input_file_name);
    TFile* input_file = new TFile(input_file_path, "READ");
    
    TTree* input_tree = (TTree*) input_file->Get("Pythia_Tree");
    
    std::cout << "Input file read. TTree accessed." << std::endl;
    
    // Variables
    int     particle_n;
    float   particle_pid[200];
    
    input_tree->SetBranchAddress("particle_n",      &particle_n);
    input_tree->SetBranchAddress("particle_pid",    particle_pid);
    
    int particle_pid_codes[] = {11, -11, 211, -211, 321, -321, 2212, -2212};
    int particle_pid_count[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    TH1F* th1d_particle_pid = new TH1F("th1d_Pythia_particle_pid", "Particle ID; PID; Number of Particles", 9, 0, 9);
    
    int event_n = input_tree->GetEntries();
    
    for ( int e = 0; e < event_n ; e++ ) {
        input_tree->GetEntry(e);
        
        // Fills tree with particle info
        for (int p = 0; p < particle_n; p++) {
            // Fills particle ID
            bool p_is_other = true;
            for (int i = 0; i < 8; i++) {
                if (particle_pid[p] == particle_pid_codes[i]) {
                    particle_pid_count[i] += 1;
                    th1d_particle_pid->AddBinContent(i+1, 1);
                    p_is_other = false;
                }
            }
            if (p_is_other) {
                particle_pid_count[8] += 1;
                th1d_particle_pid->AddBinContent(9, 1); // particle_pid_count[8] += 1;
            }
        }
    }
    
    std::cout << "Particle ID Counts: " << particle_pid_count[0] << ", " << particle_pid_count[1] << ", " << particle_pid_count[2] << ", " << particle_pid_count[3] << ", " << particle_pid_count[4] << ", " << particle_pid_count[5] << ", " << particle_pid_count[6] << ", " << particle_pid_count[7] << ", " << particle_pid_count[8] << ", " << std::endl;
    
    th1d_particle_pid -> Scale(1, "width");
    th1d_particle_pid -> Write();
    
    // Make hist_particle_pid
    char* particle_pid_names[9] = {"e^{-}","e^{+}","#pi^{+}","#pi^{-}","K^{+}","K^{-}","p^{+}","p^{-}","other"};
    TCanvas* canvas_pid = new TCanvas("canvas", "", 800, 600);

    gPad->SetTicks();
    gStyle->SetOptStat(0);

    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);

    th1d_particle_pid->SetLineColor(particle_line_color);
    th1d_particle_pid->SetMarkerColor(particle_marker_color);
    th1d_particle_pid->SetMarkerStyle(particle_marker_style);

    gPad->SetLogy(1);
    th1d_particle_pid->SetMinimum(1.);
    for (int i = 0; i < 9; i++) {
        th1d_particle_pid->GetXaxis()->SetBinLabel(i + 1, particle_pid_names[i]);
    }
    th1d_particle_pid->Draw();

    latex->DrawLatex(0.65, (0.75 + 0.06), Form("#scale[0.6]{#bf{PYTHIA8, n = %i}}", event_n));
    latex->DrawLatex(0.65, (0.75 + 0.02), "#scale[0.6]{#bf{p+p at 2.76 TeV}}");
    latex->DrawLatex(0.65, (0.75 - 0.02), "#scale[0.6]{#bf{#Sigma p_{T} > 50 GeV}}");
    latex->DrawLatex(0.65, (0.75 - 0.06), "#scale[0.6]{#bf{p_{T, jet} > 20 GeV}}");

    canvas_pid->Print(Form("%s/th1d_Pythia_particle_pid.pdf", dir_plots)); // Resets the canvas after printing
    gPad->SetLogy(0);
}





void Jet_Analyzer( char* input_file_name, char* input_tree_name, int label_arr_size, char* label_arr[10], bool bool_make_plots ) {
    char jet_dir[200];
    sprintf(jet_dir, "%s/Plots_Jets", dir_plots);
    std::__fs::filesystem::create_directories(jet_dir);
    
    // Sets file prefix based on input tree name
    char file_prefix[4];
    char file_title[100];
    if ( input_tree_name == "Pythia" ) {
        sprintf(file_prefix, "Pyth");
        sprintf(file_title, "PYTHIA Jets");
    }
    else if ( input_tree_name == "Combined" ) {
        sprintf(file_prefix, "Comb");
        sprintf(file_title, "PYTHIA + Thermal Jets");
    }
    
    // Opens and reads the Root file
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s_%s_Trees.root", dir_data, file_prefix, input_file_name);
//    sprintf(input_file_path, "%s/Jet_Comb_40_60_Test_Trees.root", dir_data);
    TFile* input_file = new TFile(input_file_path, "READ");
    
    char input_tree_path[200];
    sprintf(input_tree_path, "FastJet_Tree");
    TTree* input_tree = (TTree*) input_file->Get(input_tree_path);
    
    // Opens and makes the output file
    char output_file_path[200];
    sprintf(output_file_path, "%s/%s_%s_Plots.root", dir_data, file_prefix, input_file_name);
    TFile* output_file = new TFile(output_file_path, "UPDATE");
    
    std::cout << "Input file read. TTree accessed." << std::endl;
    
    int     jet_n;
    float   jet_pt[100];
    float   jet_y[100];
    float   jet_phi[100];
    float   jet_area[100];
    
    input_tree->SetBranchAddress("jet_n",   &jet_n);
    input_tree->SetBranchAddress("jet_pt",  jet_pt);
    input_tree->SetBranchAddress("jet_y",   jet_y);
    input_tree->SetBranchAddress("jet_phi", jet_phi);
    input_tree->SetBranchAddress("jet_area", jet_area);
    
    char plot_id[200];
    char plot_title[200];
    
    // Standard distribution plots
    sprintf(plot_id, "th1d_%s_jet_n", input_tree_name);
    sprintf(plot_title, "%s, Jet Count; N_{jet}; N_{event}", file_title);
    TH1D* th1d_jet_n = new TH1D(plot_id, plot_title, 6, 0.5, 6.5);
    
    sprintf(plot_id, "th1d_%s_jet_pt", input_tree_name);
    sprintf(plot_title, "%s, Jet p_{T}; p_{T, jet} [GeV/c]; dN_{jet}/dp_{T}", file_title);
    TH1D* th1d_jet_pt = new TH1D(plot_id, plot_title, 40, 0, 200);
    
    sprintf(plot_id, "th1d_%s_jet_y", input_tree_name);
    sprintf(plot_title, "%s, Jet Rapidity; y_{jet}; dN_{jet}/dy", file_title);
    TH1D* th1d_jet_y = new TH1D(plot_id, plot_title, 20, -0.9, 0.9);
    
    sprintf(plot_id, "th1d_%s_jet_phi", input_tree_name);
    sprintf(plot_title, "%s, Jet Angle; #phi_{jet} [rad]; dN_{jet}/d#phi", file_title);
    TH1D* th1d_jet_phi = new TH1D(plot_id, plot_title, 20, 0, 2*math_pi);
    
    sprintf(plot_id, "th1d_%s_jet_area", input_tree_name);
    sprintf(plot_title, "%s, Jet Area; Area [units^{2}]; dN_{jet}/d#phi", file_title);
    TH1D* th1d_jet_area = new TH1D(plot_id, plot_title, 20, 0, 1);
    
    th1d_jet_n      ->Sumw2();
    th1d_jet_pt     ->Sumw2();
    th1d_jet_y      ->Sumw2();
    th1d_jet_phi    ->Sumw2();
    th1d_jet_area   ->Sumw2();
    
    // Iterates through events
    for (int e = 0; e < input_tree->GetEntries(); e++) {
        input_tree->GetEntry(e);
        
        th1d_jet_n->Fill(jet_n);
        
        // Fills tree with jet info
        for (int j = 0; j < jet_n; j++) {
            // Basic Jet Info
            th1d_jet_pt     ->Fill(jet_pt[j]);
            th1d_jet_y      ->Fill(jet_y[j]);
            th1d_jet_phi    ->Fill(jet_phi[j]);
            th1d_jet_area   ->Fill(jet_area[j]);
        }
    }
    
    // Normalizes histograms
    th1d_jet_pt     ->Scale(1, "width");
    th1d_jet_y      ->Scale(1, "width");
    th1d_jet_phi    ->Scale(1, "width");
    th1d_jet_area   ->Scale(1, "width");
    
    // Writes histograms to output file
    th1d_jet_n      ->Write("", TObject::kOverwrite);
    th1d_jet_pt     ->Write("", TObject::kOverwrite);
    th1d_jet_y      ->Write("", TObject::kOverwrite);
    th1d_jet_phi    ->Write("", TObject::kOverwrite);
    th1d_jet_area   ->Write("", TObject::kOverwrite);
    
    // --- PLOT JET HISTOGRAMS ---
    
    if ( bool_make_plots ) {
        
        char plot_name[200];
    
        // Make th1d_jet_n
        sprintf(plot_name, "Plots_Jets/th1d_%s_jet_n.pdf", input_tree_name);
        th1d_plotter(
            th1d_jet_n,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            0, false,
            label_arr, label_arr_size, 0.65, 0.75);
        
        // Make th1d_jet_pt
        sprintf(plot_name, "Plots_Jets/th1d_%s_jet_pt.pdf", input_tree_name);
        th1d_plotter(
            th1d_jet_pt,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            10, true,
            label_arr, label_arr_size, 0.65, 0.75);
        
        // Make th1d_jet_y
        sprintf(plot_name, "Plots_Jets/th1d_%s_jet_y.pdf", input_tree_name);
        th1d_plotter(
            th1d_jet_y,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            0, false,
            label_arr, label_arr_size, 0.65, 0.25);
        
        // Make th1d_jet_phi
        sprintf(plot_name, "Plots_Jets/th1d_%s_jet_phi.pdf", input_tree_name);
        th1d_plotter(
            th1d_jet_phi,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            0, false,
            label_arr, label_arr_size, 0.65, 0.25);
        
        // Make hist_jet_area
        sprintf(plot_name, "Plots_Jets/th1d_%s_jet_area.pdf", input_tree_name);
        th1d_plotter(
            th1d_jet_area,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            0, false,
            label_arr, label_arr_size, 0.65, 0.25);
    }
    
    // Closes files
    output_file->Close();
    input_file->Close();
    
    delete output_file;
    delete input_file;
}



// MODIFY THIS!!!!
// Want to look at PYTHIA and Combined jets - cross reference to study top truth jets
// Maybe compare truth top 2 to top 2 combine jet distributions?
void Jet_Analyzer_pt_y_Cut( char* input_file_name, int label_arr_size, char* label_arr[10], bool bool_make_plots ) {
    
    char jet_dir[200];
    sprintf(jet_dir, "%s/Plots_Jets/Hard_Jets", dir_plots);
    std::__fs::filesystem::create_directories(jet_dir);
    
    // Opens and reads the Root file
    char pythia_file_path[200];
    sprintf(pythia_file_path, "%s/Pyth_%s_Trees.root", dir_data, input_file_name);
    TFile* pythia_file = new TFile(pythia_file_path, "READ");
    TTree* pythia_tree = (TTree*) pythia_file->Get("FastJet_Tree");
    
    // Opens and reads the Root file
    char combined_file_path[200];
    sprintf(combined_file_path, "%s/Comb_%s_Trees.root", dir_data, input_file_name);
    TFile* combined_file = new TFile(combined_file_path, "READ");
    TTree* combined_tree = (TTree*) combined_file->Get("FastJet_Tree");
    
    // Opens and makes the output file
    char output_file_path[200];
    sprintf(output_file_path, "%s/Comb+Pyth_%s_Plots.root", dir_data, input_file_name);
    TFile* output_file = new TFile(output_file_path, "UPDATE");
    
    std::cout << "Input file read. TTree accessed." << std::endl;
    
    int     pyth_jet_n;
    float   pyth_jet_pt[100];
    float   pyth_jet_y[100];
    float   pyth_jet_phi[100];
    float   pyth_jet_area[100];
    
    pythia_tree->SetBranchAddress("jet_n",    &pyth_jet_n);
    pythia_tree->SetBranchAddress("jet_pt",   pyth_jet_pt);
    pythia_tree->SetBranchAddress("jet_y",    pyth_jet_y);
    pythia_tree->SetBranchAddress("jet_phi",  pyth_jet_phi);
    pythia_tree->SetBranchAddress("jet_area", pyth_jet_area);
    
    int     comb_jet_n;
    float   comb_jet_pt[100];
    float   comb_jet_y[100];
    float   comb_jet_phi[100];
    float   comb_jet_area[100];
    
    combined_tree->SetBranchAddress("jet_n",    &comb_jet_n);
    combined_tree->SetBranchAddress("jet_pt",   comb_jet_pt);
    combined_tree->SetBranchAddress("jet_y",    comb_jet_y);
    combined_tree->SetBranchAddress("jet_phi",  comb_jet_phi);
    combined_tree->SetBranchAddress("jet_area", comb_jet_area);
    
    char plot_id[200];
    char plot_title[200];
    
    // Plots for rapidity and jet pT cuts
    sprintf(plot_id, "th1d_%s_jet_pt_09", "Comb+Pyth");
    sprintf(plot_title, "Jet p_{T}^{Raw} for |y| < 0.9; p_{T, jet}^{Raw} [GeV/c]; dN_{jet}/dp_{T}");
    TH1D* th1d_jet_pt_09 = new TH1D(plot_id, plot_title, 40, 0, 200);
    
    sprintf(plot_id, "th1d_%s_jet_pt_09_t2", "Comb+Pyth");
    sprintf(plot_title, "Jet p_{T}^{Raw} for |y| < 0.9, Two Highest-p_{T} Jets; p_{T, jet}^{Raw} [GeV/c]; dN_{jet}/dp_{T}");
    TH1D* th1d_jet_pt_09_t2 = new TH1D(plot_id, plot_title, 40, 0, 200);
    
    sprintf(plot_id, "th1d_%s_jet_pt_05", "Comb+Pyth");
    sprintf(plot_title, "Jet p_{T}^{Raw} for |y| < 0.5; p_{T, jet}^{Raw} [GeV/c]; dN_{jet}/dp_{T}");
    TH1D* th1d_jet_pt_05 = new TH1D(plot_id, plot_title, 40, 0, 200);
    
    sprintf(plot_id, "th1d_%s_jet_pt_05_t2", "Comb+Pyth");
    sprintf(plot_title, "Jet p_{T}^{Raw} for |y| < 0.5, Two Highest-p_{T} Jets; p_{T, jet}^{Raw} [GeV/c]; dN_{jet}/dp_{T}");
    TH1D* th1d_jet_pt_05_t2 = new TH1D(plot_id, plot_title, 40, 0, 200);
    
    sprintf(plot_id, "th1d_%s_jet_pt_05", "Pyth");
    sprintf(plot_title, "Jet p_{T}^{True} for |y| < 0.9, Two Highest-p_{T} Jets; p_{T, jet}^{True} [GeV/c]; dN_{jet}/dp_{T}");
    TH1D* th1d_pyth_jet_pt_09_t2 = new TH1D(plot_id, plot_title, 40, 0, 200);
    
    sprintf(plot_id, "th1d_%s_jet_pt_05_t2", "Pyth");
    sprintf(plot_title, "Jet p_{T}^{True} for |y| < 0.5, Two Highest-p_{T} Jets; p_{T, jet}^{True} [GeV/c]; dN_{jet}/dp_{T}");
    TH1D* th1d_pyth_jet_pt_05_t2 = new TH1D(plot_id, plot_title, 40, 0, 200);
    
    th1d_jet_pt_09      ->Sumw2();
    th1d_jet_pt_09_t2   ->Sumw2();
    th1d_jet_pt_05      ->Sumw2();
    th1d_jet_pt_05_t2   ->Sumw2();
    th1d_pyth_jet_pt_09_t2 ->Sumw2();
    th1d_pyth_jet_pt_05_t2 ->Sumw2();
    
    // Iterates through events
    for (int e = 0; e < combined_tree->GetEntries(); e++) {
        combined_tree->GetEntry(e);
        pythia_tree->GetEntry(e);
        
        // Fills tree with jet info
        for (int j = 0; j < comb_jet_n; j++) {
            // checks for a matching truth jet
            int pythia_match = -1;
            for ( int k = 0 ; k < pyth_jet_n ; k++ ) {
                if ( ( pow((pyth_jet_y[k] - comb_jet_y[j]), 2) + pow((pyth_jet_phi[k] - comb_jet_phi[j]), 2) ) < fj_rSquared) {
                    pythia_match = k;
                }
            }
            if ( pythia_match >= 0 ) {
                // Jet y and pT cuts
                if ( comb_jet_y[j] < 0.9 && comb_jet_y[j] > -0.9 ) {
                    if (j < 2) {
                        th1d_jet_pt_09_t2->Fill(comb_jet_pt[j]);
                        th1d_pyth_jet_pt_09_t2->Fill(pyth_jet_pt[pythia_match]);
                    }
                    th1d_jet_pt_09->Fill(comb_jet_pt[j]);
                }
                if ( comb_jet_y[j] < 0.5 && comb_jet_y[j] > -0.5  ) {
                    if (j < 2) {
                        th1d_jet_pt_05_t2->Fill(comb_jet_pt[j]);
                        th1d_pyth_jet_pt_05_t2->Fill(pyth_jet_pt[pythia_match]);
                    }
                    th1d_jet_pt_05->Fill(comb_jet_pt[j]);
                }
            }
        }
    }
    
    // Normalizes histograms
    th1d_jet_pt_09      ->Scale(1, "width");
    th1d_jet_pt_09_t2   ->Scale(1, "width");
    th1d_jet_pt_05      ->Scale(1, "width");
    th1d_jet_pt_05_t2   ->Scale(1, "width");
    th1d_pyth_jet_pt_09_t2 ->Scale(1, "width");
    th1d_pyth_jet_pt_05_t2 ->Scale(1, "width");
    
    // Writes histograms to output file
    th1d_jet_pt_09      ->Write("", TObject::kOverwrite);
    th1d_jet_pt_09_t2   ->Write("", TObject::kOverwrite);
    th1d_jet_pt_05      ->Write("", TObject::kOverwrite);
    th1d_jet_pt_05_t2   ->Write("", TObject::kOverwrite);
    th1d_pyth_jet_pt_09_t2 ->Write("", TObject::kOverwrite);
    th1d_pyth_jet_pt_05_t2 ->Write("", TObject::kOverwrite);
    
    // --- PLOT JET HISTOGRAMS ---
    
    if ( bool_make_plots ) {
        
        char plot_name[200];
        
        // Make th1d_jet_pt_09
        sprintf(plot_name, "Plots_Jets/Hard_Jets/th1d_%s_jet_pt_09.pdf", input_file_name);
        th1d_plotter(
            th1d_jet_pt_09,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            10, true,
            label_arr, label_arr_size, 0.60, 0.35);
        
        // Make th1d_jet_pt_09_t2
        sprintf(plot_name, "Plots_Jets/Hard_Jets/th1d_%s_jet_pt_09_t2.pdf", input_file_name);
        th1d_plotter(
            th1d_jet_pt_09_t2,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            10, true,
            label_arr, label_arr_size, 0.60, 0.35);
        
        // Make th1d_jet_pt_05
        sprintf(plot_name, "Plots_Jets/Hard_Jets/th1d_%s_jet_pt_05.pdf", input_file_name);
        th1d_plotter(
            th1d_jet_pt_05,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            10, true,
            label_arr, label_arr_size, 0.60, 0.35);
        
        // Make th1d_jet_pt_05_t2
        sprintf(plot_name, "Plots_Jets/Hard_Jets/th1d_%s_jet_pt_05_t2.pdf", input_file_name);
        th1d_plotter(
            th1d_jet_pt_05_t2,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            10, true,
            label_arr, label_arr_size, 0.60, 0.35);
        
        // Make th1d_pyth_jet_pt_09_t2
        sprintf(plot_name, "Plots_Jets/Hard_Jets/th1d_%s_pyth_jet_pt_09_t2.pdf", input_file_name);
        th1d_plotter(
            th1d_pyth_jet_pt_09_t2,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            10, true,
            label_arr, label_arr_size, 0.60, 0.35);
        
        // Make th1d_pyth_jet_pt_05_t2
        sprintf(plot_name, "Plots_Jets/Hard_Jets/th1d_%s_pyth_jet_pt_05_t2.pdf", input_file_name);
        th1d_plotter(
            th1d_pyth_jet_pt_05_t2,
            plot_name,
            jet_line_color, jet_marker_color, jet_marker_style, 1.,
            10, true,
            label_arr, label_arr_size, 0.60, 0.35);
    }
    
    // Closes files
    output_file->Close();
    pythia_file->Close();
    combined_file->Close();
    
    delete output_file;
    delete pythia_file;
    delete combined_file;
}





void Event_Analyzer( char* input_file_name, int label_arr_size, char* label_arr[10] ) {
    
//    Particle_Analyzer(input_file_name, "Pythia",   label_arr_size, label_arr, true);
//    Particle_Analyzer(input_file_name, "Thermal",  label_arr_size, label_arr, true);
//    Particle_Analyzer(input_file_name, "Combined", label_arr_size, label_arr, true);
//    Jet_Analyzer(input_file_name, "Pythia",   label_arr_size, label_arr, true);
    Jet_Analyzer_pt_y_Cut(input_file_name, label_arr_size, label_arr, true);
}





void Event_Analyzer_ROOT() {
    
    const int label_arr_size = 4;
    char* label_arr[label_arr_size] = {
        "PYTHIA + Thermal Background", "Modeling p+p at 2.76 TeV", "10 < p_{T, jet}^{True} < 90 GeV", "p_{T}HatMin > 7.5 GeV",};
    Event_Analyzer("10_90_Train", label_arr_size, label_arr);
//    Event_Analyzer("40_60_Train", label_arr_size, label_arr);
}
