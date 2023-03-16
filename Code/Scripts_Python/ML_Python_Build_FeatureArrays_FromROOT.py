import ROOT
from os import path
import numpy as np
from statistics import mean, median
import math
from array import array
import csv



def Build_FeatureArrays_FromROOT(
    input_file_path, 
    input_tree_name,
    csv_file_path,
    pt_true_min, 
    pt_true_max,
    top_n_jets = 100
    ) :
    """
    Creates arrays of all features exported from ROOT.
    Applies a cut using pT_True values between pt_true_min and pt_true_max.
    
    Returns:
    X_values       Array of arrays of input features (one array of features per jet)
    y_values       Array of target truth values
    sc_values_arr  Array of simple correction values
    """
    
    # Checks for existence of input file
    input_file = None;
    if (ROOT.gSystem.AccessPathName(input_file_path)) :
        print("Input file path does not exist:", input_file)
        sys.exit()
    else :
        input_file = ROOT.TFile.Open(input_file_path, "READ")
        print("Input file accessed successfully. Output file generated.")
    
    print("Accessing input tree...")
    input_tree = input_file.Get(input_tree_name)
    print("Input tree accessed successfully.")
    
    print("Creating .csv backup file...")
    output_csv  = open(csv_file_path, 'w', newline='')
    csv_writer  = csv.writer(output_csv)
    csv_header = [  'jet_pt_raw',      'jet_pt_corr',     'jet_mass',        'jet_area', 
                    'jet_area_err',    'jet_const_n',     'const_pt_mean',   'const_pt_median', 
                    'const_1_pt',      'const_2_pt',      'const_3_pt',      'const_4_pt',
                    'const_5_pt',      'const_6_pt',      'const_7_pt',      'const_8_pt',
                    'const_9_pt',      'const_10_pt',     'jet_y',           'jet_phi',
                    'jet_rho',         'jet_pt_true',     'jet_index']
    csv_writer.writerow(csv_header)
    print("Backup file started.")
    
    # Setup Arrays
    X_values  = []  # Array of arrays of inputs corresponding to pT_true as PYTHIA jet pT
    y_values  = []  # Array of targets for regression, pT_true is PYTHIA jet pT

    # Predictors
    jet_pt_raw       = None  # Raw/uncorrected jet pt
    jet_pt_corr      = None  # Corrected jet pt
    jet_mass         = None
    jet_area         = None
    jet_area_err     = None
    jet_const_n      = None
    const_pt_mean    = None  # Mean pt of jet constituents
    const_pt_median  = None  # Mean pt of jet constituents
    const_1_pt       = None  # pt of jet constituent particle 1
    const_2_pt       = None  # pt of jet constituent particle 2
    const_3_pt       = None  # pt of jet constituent particle 3
    const_4_pt       = None  # pt of jet constituent particle 4
    const_5_pt       = None  # pt of jet constituent particle 5
    const_6_pt       = None  # pt of jet constituent particle 6
    const_7_pt       = None  # pt of jet constituent particle 7
    const_8_pt       = None  # pt of jet constituent particle 8
    const_9_pt       = None  # pt of jet constituent particle 9
    const_10_pt      = None  # pt of jet constituent particle 10
    jet_y            = None
    jet_phi          = None
    jet_rho          = None

    # Targets
    jet_pt_true      = None  # True jet pt (determined from PYTHIA jets)

    # Helper Variables
    jet_counter      = 0     # Counts total number of jets in file
    jet_index        = 0     # Index of each jet in its event
    jet_const_pt_arr = []    # Array of jet constituents and their values
    sc_values        = []    # Array of simple correction values
    
    print("Preparing to collect data from TTree...")

    # Collecting from TTree
    for jet in input_tree :

        jet_index       = input_tree.jet_index
        if jet_index > top_n_jets: continue

        jet_pt_raw      = input_tree.jet_pt_raw
        jet_pt_corr     = input_tree.jet_pt_corr
        jet_mass        = input_tree.jet_mass
        jet_area        = input_tree.jet_area
        jet_area_err    = input_tree.jet_area_err
        jet_const_n     = input_tree.jet_const_n
        const_pt_mean   = input_tree.const_pt_mean
        const_pt_median = input_tree.const_pt_median
        const_1_pt      = input_tree.const_1_pt
        const_2_pt      = input_tree.const_2_pt
        const_3_pt      = input_tree.const_3_pt
        const_4_pt      = input_tree.const_4_pt
        const_5_pt      = input_tree.const_5_pt
        const_6_pt      = input_tree.const_6_pt
        const_7_pt      = input_tree.const_7_pt
        const_8_pt      = input_tree.const_8_pt
        const_9_pt      = input_tree.const_9_pt
        const_10_pt     = input_tree.const_10_pt
        jet_y           = input_tree.jet_y
        jet_phi         = input_tree.jet_phi
        jet_rho         = input_tree.jet_rho

        jet_pt_true     = input_tree.jet_pt_true

        temp_jet_arr = [
                jet_pt_raw,      jet_pt_corr,     jet_mass,        jet_area, 
                jet_area_err,    jet_const_n,     const_pt_mean,   const_pt_median, 
                const_1_pt,      const_2_pt,      const_3_pt,      const_4_pt,
                const_5_pt,      const_6_pt,      const_7_pt,      const_8_pt,
                const_9_pt,      const_10_pt,     jet_y,           jet_phi,
                jet_rho]

        if (jet_pt_true != 0.0) and (jet_pt_true > pt_true_min) and (jet_pt_true < pt_true_max) :
            X_values.append(temp_jet_arr)

            y_values.append(jet_pt_true)

            sc_values.append(jet_pt_corr)

            # Writes data to a .csv file for faster recall
            csv_row = temp_jet_arr
            csv_row.append(jet_pt_true)
            csv_row.append(jet_index)
            csv_writer.writerow(csv_row)

            jet_counter   = jet_counter + 1

            if jet_counter % 10000 == 0 : print(f"Jet: {jet_counter:2.0f} | pTraw: {jet_pt_raw:3.3f} | pTcorr: {jet_pt_corr: 3.3f} | pTtrue: {jet_pt_true: 5.3f}")
        
    output_csv.close()
    print("Backup .csv file closed.")

    print(f"All data transferred to array. Testing with {jet_counter} jets.\n")
    print(f"Data set lengths: {len(X_values)} / {len(y_values)} / {len(sc_values)}")

    input_file.Close()
    print("Input file closed.")
    
    return X_values, y_values, sc_values





def Build_FeatureArrays_FromROOT_ByEvent(
    input_file_path, 
    input_tree_name,
    csv_file_path,
    pt_true_min, 
    pt_true_max,
    top_n_jets = 100
    ) :
    """
    Creates arrays of all features exported from ROOT.
    Applies a cut using pT_True values between pt_true_min and pt_true_max.
    
    Returns:
    X_values       Array of arrays of input features (one array of features per jet)
    y_values       Array of target truth values
    sc_values_arr  Array of simple correction values
    """
    
    # Checks for existence of input file
    input_file = None;
    if (ROOT.gSystem.AccessPathName(input_file_path)) :
        print("Input file path does not exist:", input_file)
        return
    else :
        input_file = ROOT.TFile.Open(input_file_path, "READ")
        print("Input file accessed successfully. Output file generated.")
    
    print("Accessing input tree...")
    input_tree = input_file.Get(input_tree_name)
    print("Input tree accessed successfully.")
    
    print("Creating .csv backup file...")
    output_csv  = open(csv_file_path, 'w', newline='')
    csv_writer  = csv.writer(output_csv)
    csv_header = [  'jet_pt_raw',      'jet_pt_corr',     'jet_mass',        'jet_area', 
                    'jet_area_err',    'jet_const_n',     'const_pt_mean',   'const_pt_median', 
                    'const_1_pt',      'const_2_pt',      'const_3_pt',      'const_4_pt',
                    'const_5_pt',      'const_6_pt',      'const_7_pt',      'const_8_pt',
                    'const_9_pt',      'const_10_pt',     'jet_y',           'jet_phi',
                    'jet_rho',         'jet_pt_true',     'jet_index']
    csv_writer.writerow(csv_header)
    print("Backup file started.")
    
    # Setup Arrays
    X_values  = []  # Array of arrays of inputs corresponding to pT_true as PYTHIA jet pT
    y_values  = []  # Array of targets for regression, pT_true is PYTHIA jet pT

    # Predictors
    jet_pt_raw       = None  # Raw/uncorrected jet pt
    jet_pt_corr      = None  # Corrected jet pt
    jet_mass         = None
    jet_area         = None
    jet_area_err     = None
    jet_const_n      = None
    const_pt_mean    = None  # Mean pt of jet constituents
    const_pt_median  = None  # Mean pt of jet constituents
    const_1_pt       = None  # pt of jet constituent particle 1
    const_2_pt       = None  # pt of jet constituent particle 2
    const_3_pt       = None  # pt of jet constituent particle 3
    const_4_pt       = None  # pt of jet constituent particle 4
    const_5_pt       = None  # pt of jet constituent particle 5
    const_6_pt       = None  # pt of jet constituent particle 6
    const_7_pt       = None  # pt of jet constituent particle 7
    const_8_pt       = None  # pt of jet constituent particle 8
    const_9_pt       = None  # pt of jet constituent particle 9
    const_10_pt      = None  # pt of jet constituent particle 10
    jet_y            = None
    jet_phi          = None
    jet_rho          = None

    # Targets
    jet_pt_true      = None  # True jet pt (determined from PYTHIA jets)

    # Helper Variables
    jet_counter      = 0     # Counts total number of jets in file
    jet_index        = 0     # Index of each jet in its event
    jet_const_pt_arr = []    # Array of jet constituents and their values
    sc_values        = []    # Array of simple correction values
    
    print("Preparing to collect data from TTree...")

    # Collecting from TTree
    for event in input_tree :
        jet_n = event.jet_n
        
        for jet in range(0, jet_n):
            jet_pt_raw      = input_tree.jet_pt_raw[jet]
            jet_pt_corr     = input_tree.jet_pt_corr[jet]
            jet_mass        = input_tree.jet_mass[jet]
            jet_area        = input_tree.jet_area[jet]
            jet_area_err    = 0
            jet_const_n     = input_tree.jet_const_n[jet]
            const_pt_mean   = input_tree.const_pt_mean[jet]
            const_pt_median = input_tree.const_pt_median[jet]
            const_1_pt      = input_tree.const_1_pt[jet]
            const_2_pt      = input_tree.const_2_pt[jet]
            const_3_pt      = input_tree.const_3_pt[jet]
            const_4_pt      = input_tree.const_4_pt[jet]
            const_5_pt      = input_tree.const_5_pt[jet]
            const_6_pt      = input_tree.const_6_pt[jet]
            const_7_pt      = input_tree.const_7_pt[jet]
            const_8_pt      = input_tree.const_8_pt[jet]
            const_9_pt      = input_tree.const_9_pt[jet]
            const_10_pt     = input_tree.const_10_pt[jet]
            jet_y           = input_tree.jet_y[jet]
            jet_phi         = input_tree.jet_phi[jet]
            jet_rho         = input_tree.jet_rho[jet]

            jet_pt_true     = input_tree.jet_pt_true_pythia[jet]

            temp_jet_arr = [
                    jet_pt_raw,      jet_pt_corr,     jet_mass,        jet_area, 
                    jet_area_err,    jet_const_n,     const_pt_mean,   const_pt_median, 
                    const_1_pt,      const_2_pt,      const_3_pt,      const_4_pt,
                    const_5_pt,      const_6_pt,      const_7_pt,      const_8_pt,
                    const_9_pt,      const_10_pt,     jet_y,           jet_phi,
                    jet_rho]

            X_values.append(temp_jet_arr)

            y_values.append(jet_pt_true)

            sc_values.append(jet_pt_corr)

            # Writes data to a .csv file for faster recall
            csv_row = temp_jet_arr
            csv_row.append(jet_pt_true)
            csv_row.append(jet_index)
            csv_writer.writerow(csv_row)

            jet_counter   = jet_counter + 1

            if jet_counter % 10000 == 0 : print(f"Jet: {jet_counter:2.0f} | pTraw: {jet_pt_raw:3.3f} | pTcorr: {jet_pt_corr: 3.3f} | pTtrue: {jet_pt_true: 5.3f}")
        
    output_csv.close()
    print("Backup .csv file closed.")

    print(f"All data transferred to array. Testing with {jet_counter} jets.\n")
    print(f"Data set lengths: {len(X_values)} / {len(y_values)} / {len(sc_values)}")

    input_file.Close()
    print("Input file closed.")
    
    return X_values, y_values, sc_values
