from os import (path, mkdir)
import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from statistics import mean, median
from array import array
from datetime import datetime
from joblib import dump, load
import numpy as np
import math
import csv




def Build_FeatureArrays_FromCSV(
    csv_file_path
    ) :
    """
    Simple function to build an array from a csv that has the following columns (number indicates index):
    
    0  jet_pt_raw,      1  jet_pt_corr,     2  jet_mass,        3  jet_area, 
    4  jet_area_err,    5  jet_const_n,     6  const_pt_mean,   7  const_pt_median, 
    8  const_1_pt,      9  const_2_pt,      10 const_3_pt,      11 const_4_pt,
    12 const_5_pt,      13 const_6_pt,      14 const_7_pt,      15 const_8_pt,
    16 const_9_pt,      17 const_10_pt,     18 jet_y,           19 jet_phi,
    20 jet_rho,         21 jet_pt_true,     22 jet_index
    """
    
    input_csv  = open(csv_file_path)
    csv_reader = csv.reader(input_csv)
    
    X_values  = []
    y_values  = []
    sc_values = []
    
    skip_header_toggle = True  # Skips the header row of the csv
    jet_counter = 0            # Counts total number of jets in file
    
    print("Preparing to collect data from csv backup file...")
    
    for row in csv_reader:
        if skip_header_toggle:
            skip_header_toggle = False
            continue
        X_values.append([
            float(row[0]),  float(row[1]),  float(row[2]),  float(row[3]), 
            float(row[4]),  float(row[5]),  float(row[6]),  float(row[7]),
            float(row[8]),  float(row[9]),  float(row[10]), float(row[11]),
            float(row[12]), float(row[13]), float(row[14]), float(row[15]),
            float(row[16]), float(row[17]), float(row[18]), float(row[19]),
            float(row[20])])
        y_values.append(float(row[21]))
        sc_values.append(float(row[1]))
        
        jet_counter   = jet_counter + 1
        
#        if jet_counter % 10000 == 0 : print(f"Jet: {jet_counter:2.0f} | pTraw: {float(row[0]):3.3f} | pTcorr: {float(row[1]): 3.3f} | pTtrue: {float(row[21]): 5.3f}")
    
    print("Data collected!")
    input_csv.close()
    print("Backup .csv file closed.")

    print(f"All data transferred to array. Testing with {jet_counter} jets.\n")
    print(f"Data set lengths: {len(X_values)} / {len(y_values)} / {len(sc_values)}")
    
    return X_values, y_values, sc_values



def Build_SelectFeatureArray(
    X_features,
    feature_index
    ) :
    """
    Builds training and testing data sets
    """
    
    print("Selecting data from master array...")
    
    X_features_select = []
    for i in range(len(X_features)):
        X_temp = []
        for j in range(len(feature_index)):
            X_temp.append(X_features[i][feature_index[j]])
        X_features_select.append(X_temp)
        
    print("Data ready. Feature array length:", len(X_features_select), "\n")
    
    return X_features_select



def Write_MLResults_ToCSV(
    output_filename,  # File path with name (MUST include '.csv' at end)
    pt_true_array,    # Array of pt_true values (probably 'y_test' array)
    pt_sc_array,      # Array of pt values from rho*area method (probably 'sc_test' array)
    pt_reco_array,    # Array of pt_reco values from ML (probably 'results' array)
    feature_array,    # Array of pt_true values (probably 'X_test_select' array)
    feature_label     # Array of feature labels
    ) :
    """
    Writes ML results to a csv file
    """
    # Create csv file
    output_csv  = open(output_filename, 'w', newline='')
    csv_writer  = csv.writer(output_csv)
    
    # Build and write header
    csv_header = ['jet_pt_true', 'jet_pt_reco']
    no_pt_corr = False
    
    if ( 'jet_pt_corr' not in feature_label ):
        csv_header.append('jet_pt_corr')
        no_pt_corr = True
    
    for label in feature_label:
        csv_header.append(label)
    csv_writer.writerow(csv_header)
    
    # Add data
    for i in range(len(pt_true_array)):
        csv_row = [pt_true_array[i], pt_reco_array[i]]
        if no_pt_corr:
            csv_row.append(pt_sc_array[i])
        for feature in feature_array[i]:
            csv_row.append(feature)
        csv_writer.writerow(csv_row)
        
    output_csv.close()
    print("ML results .csv file closed.")
    
    return


def Write_MLCoefficients_ToCSV(
    output_filename,  # File path with name (MUST include '.csv' at end)
    coeff_array,      # Array of feature weight values (probably 'lr_coeffs' or 'rf_importance')
    feature_label     # Array of feature labels
    ) :
    """
    Writes ML results to a csv file
    """
    print("Writing coefficients to CSV...")
    # Create csv file
    output_csv  = open(output_filename, 'w', newline='')
    csv_writer  = csv.writer(output_csv)
    
    # Build and write header
    csv_header = []
    for label in feature_label:
        csv_header.append(label)
    csv_writer.writerow(csv_header)
    
    # Add data
    csv_row = []
    for coeff in coeff_array:
        csv_row.append(coeff)
    csv_writer.writerow(csv_row)
        
    output_csv.close()
    print("Coefficient CSV file complete.")
    
    return



def Train_LinearRegression(
    X_train,            # Array of array of input features
    y_train,            # Array of target values
    features_labels,    # Array of string labels for each feature
    use_scaler = True   # If True, uses StandardScalar
    ):
    print("\n----- Fitting Linear Regression Estimator -----\n")
    
    # Sets Scaler
    if use_scaler:
        scaler = StandardScaler(with_mean=True, with_std=True)
        print("\nUsing StandardScaler. Data will be recentered and normalized.\n")
    else:
        scaler = StandardScaler(with_mean=False, with_std=False)
        print("\nData will not be rescaled.\n")
        
    # Creates a linear regression model in a pipeline
    lr_estimator = LinearRegression()
    lr_pipeline = make_pipeline(
        scaler,
        lr_estimator )
    lr_coeffs = 0
    
    # Fits the regression model
    output = lr_pipeline.fit(X_train, y_train)
    print("\nLinear Regression Fit:\n", output)

    # Outputs regression coefficients
    lr_coeffs = lr_estimator.coef_
    intercept = lr_estimator.intercept_
    print(intercept)
    print(lr_coeffs)
    lr_coeffs = np.append(lr_coeffs, intercept)
    print(lr_coeffs)

    print("Regression Coefficients:")
    for i in range(len(lr_coeffs)) :
        print(features_labels[i], lr_coeffs[i])
    
    print(type(lr_pipeline))
    
    return lr_pipeline, lr_coeffs



def Train_RandomForestRegression(
    X_train,            # Array of array of input features
    y_train,            # Array of target values
    features_labels,    # Array of string labels for each feature
    use_scaler = True,  # If True, uses StandardScalar
    n_estimators = 100,
    max_depth = 5,
    n_jobs = 3
    ):
    print("\n----- Fitting Random Forest Regression Estimator -----\n")    
    
    # Sets Scaler
    if use_scaler:
        scaler = StandardScaler(with_mean=True, with_std=True)
        print("\nUsing StandardScaler. Data will be recentered and normalized.\n")
    else:
        scaler = StandardScaler(with_mean=False, with_std=False)
        print("\nData will not be rescaled.\n")
        
    # Creates a random forest estimator in a pipeline
    rf_estimator = RandomForestRegressor(
        n_jobs = n_jobs,
        max_depth = max_depth)
    rf_pipeline = make_pipeline(
        scaler,
        rf_estimator )

    # Fits the regression model
    output = rf_pipeline.fit(X_train, y_train)
    print("\nRandom Tree Regression Fit:\n", output)

    # Outputs feature importances
    rf_features = rf_estimator.feature_importances_

    print("Feature Importance:")
    for i in range(len(rf_features)) :
        print(features_labels[i], rf_features[i])
    
    return rf_pipeline, rf_features



def Train_MLPRegression(
    X_train,                     # Array of array of input features
    y_train,                     # Array of target values
    features_labels,             # Array of string labels for each feature
    use_scaler = True,   # If True, uses StandardScalar
    max_iter = 200,         # Maximum number of iterations to use for training
    hidden_layer_sizes = 100
    ):
    print("\n----- Fitting Neural Network Regression Estimator -----\n")   
    
    # Sets Scaler
    if use_scaler:
        scaler = StandardScaler(with_mean=True, with_std=True)
        print("\nUsing StandardScaler. Data will be recentered and normalized.\n")
    else:
        scaler = StandardScaler(with_mean=False, with_std=False)
        print("\nData will not be rescaled.\n")
        
    # Creates a MLP estimator in a pipeline
    mlp_estimator = MLPRegressor(
        max_iter = max_iter,
        hidden_layer_sizes = hidden_layer_sizes)
    mlp_pipeline = make_pipeline(
        scaler,
        mlp_estimator)

    # Fits the regression model
    output = mlp_pipeline.fit(X_train, y_train)
    print("\nMultilayer Perceptron Regression Fit:\n", output)
    
    return mlp_pipeline



def Train_All_Estimators(
    X_train,             # Array of array of input features
    y_train,             # Array of target values
    features_labels,     # Array of string labels for each feature
    use_scaler = True,   # If True, uses StandardScalar
    use_lr = True,       # If True, trains a linear regression estimator
    use_rf = True,       # If True, trains a random forest regression estimator
    use_mlp = True       # If True, trains a multilayer perceptron regression estimator
    ):
    """
    Function for training Machine Learning Estimators.
    Takes in feature array (X_train) - this is actually an array of arrays, 
    a target array (y_train), and an array of feature labels (feature_arr_labels).
    
    Note that X_train and y_train MUST be the same length, 
    and feature_arr_lables should be 1 longer than the number of features used.
    
    Returns:
    lr_pipeline      Trained sklearn LinearRegression estimator
    rf_pipeline      Trained sklearn RandomForestRegressor estimator
    mlp_pipeline     Trained sklearn MLPRegressor estimator
    lr_coeffs_arr    Array of coefficients used by lr_pipeline, corresponding to each input feature
    rf_features_arr  Array of feature importances used by rf_pipeline, corresponding to each input feature
    """
    
    print("\n===== Running Train_Estimators Function =====\n")
    
    # --- LINEAR REGRESSION ---
    
    if use_lr:
        Train_LinearRegression(X_train, y_train, features_labels, use_scaler = use_scaler)
    else:
        lr_pipeline = False
        lr_coeffs_arr = False

    
    # --- RANDOM FOREST REGRESSION ---
    
    if use_rf:
        Train_RandomForestRegression(X_train, y_train, features_labels, use_scaler = use_scaler)
    else:
        rf_pipeline = False
        rf_features_arr = False
        
        
    # --- MULTILAYER PERCEPTRON REGRESSION ---
    
    if use_mlp:
        Train_MLPRegression(X_train, y_train, features_labels, use_scaler = use_scaler)
    else:
        mlp_pipeline = False
    
    return lr_pipeline, rf_pipeline, mlp_pipeline, lr_coeffs_arr, rf_features_arr

    
    
def Test_Estimator(
    ml_pipeline, # Trained sklearn machine learning estimator pipeline
    X_test,      # Array of arrays of input features
    y_test       # Array of pT truth values for comparison
    ) :
    """
    Simple function to test a trained sklearn estimator on new data.
    
    Returns:
    results        Array of predicted values from input features
    results_delta  Array of differences between prediction and truth (prediction - truth)
    """
    
    if not ml_pipeline:
        print("No pipeline provided!")
        return False, False
    
    results = ml_pipeline.predict(X_test)
    results_delta = []
    
    for j in range(len(results)):
        results_delta.append(results[j] - y_test[j])
    
    return results, results_delta

    

def TestAndSave_Estimator(
    feature_label,      # Array of labels corresponding to each feature
    feature_index,      # Array of indices for each feature used in X_train
    ml_pipeline,        # Trained Machine Learning Regression Pipeline
    X_test_select,      # Array of testing data features
    y_test,             # Array of testing data targets
    sc_test,            # Array of testing data simple corrections
    pt_test_min,        # Float of min pT to test with
    pt_test_max,        # Float of max pT to test with
    output_filename,    # Directory path + name for output csv file
    use_scaler = True   # If true, rescales data
    ) :
    
    X_test_temp  = []
    y_test_temp  = []
    sc_test_temp = []
    
    for i in range(len(y_test)):
        if y_test[i] > pt_test_min and y_test[i] < pt_test_max:
            X_test_temp.append(X_test_select[i])
            y_test_temp.append(y_test[i])
            sc_test_temp.append(sc_test[i])
        else: continue
    
    # Tests estimator
    ml_results, ml_results_delta = Test_Estimator(
        ml_pipeline,
        X_test_temp,
        y_test_temp
        )
    
    # Writes outputs to a csv file
    Write_MLResults_ToCSV(
        output_filename,
        y_test_temp,
        sc_test_temp,
        ml_results,
        X_test_temp,
        feature_label
        )
    
    return



def Test_All_Estimators(
    X_test,          # Array of arrays of input features
    y_test,          # Array of pT truth values for comparison
    lr_pipeline,     # Trained Linear Regression estimator pipeline
    rf_pipeline,     # Trained Random Forest estimator pipeline
    mlp_pipeline     # Trained Multilayer Perceptron estimator pipeline
    ) :
    """
    Function for testing the Machine Learning estimators on testing data sets.
    
    Note: The estimators used for this MUST be trained with the same features used to test!
    You cannot change (increase, decrease, or swap) the estimators between training and testing.
    
    Note: The number of datapoints for testing CAN be different between training and testing.
    For example, you can train with 100k points and test 500k points.
    
    Returns:
    lr_results         Array of predictions from LinearRegression
    lr_results_delta   Array of differences between predictions from LinearRegression and truth
    rf_results         Array of predictions from RandomForestRegression
    rf_results_delta   Array of differences between predictions from RandomForestRegression and truth
    mlp_results        Array of predictions from MLPRegression
    mlp_results_delta  Array of differences between predictions from MLPRegression and truth
    """
    
    # Initializes ML results variables
    lr_results, lr_results_delta   = False, False
    rf_results, rf_results_delta   = False, False
    mlp_results, mlp_results_delta = False, False
    
    # Tests each estimator and assigns outcomes to results
    if lr_pipeline:
        lr_results, lr_results_delta   = Test_Estimator(lr_pipeline, X_test, y_test)
    
    if rf_pipeline:
        rf_results, rf_results_delta   = Test_Estimator(rf_pipeline, X_test, y_test)
    
    if mlp_pipeline:
        mlp_results, mlp_results_delta = Test_Estimator(mlp_pipeline, X_test, y_test)
    
    return lr_results, lr_results_delta, rf_results, rf_results_delta, mlp_results, mlp_results_delta



def Full_TrainTest(
    train_file_bundle,  # Array of tuples of train files:       [(0:"File_Name.csv", 1:"Base_Name", 2:"Bias")]
    test_file_bundle,   # Array of tuples of test files:        [(0:"File_Name.csv", 1:"Base_Name", 2:"Bias")]
    feature_bundle,     # Array of tuples of feature labels:    [(0:feature_label_arr, 1:feature_index_arr)]
    test_bin_array,     # Array of tuples of test bin min/max:  [("Test_Bin_Set_Label", ((test_pt_min, test_pt_min),...))]
    traintest_bin_array, # Array of tuples to train and test over the same range:  [("Bin_Set_Label",  ((traintest_pt_min, traintest_pt_max),...))]. If traintest_bin_array = None, skips this step of the function
    output_directory,   # Output directory. Subfolders will be created accordingly
    train_pt_min,
    train_pt_max,
    use_lr  = True,     # Use linear regression. Defaults to True since LR is fast
    use_rf  = True,     # Use random forest regression. Defaults to False since RF is slow
    use_mlp = True,     # Use multi-layer perceptron regression. Defaults to False since MLP is slow
    use_lr_tt  = True,  # Use linear regression. Defaults to True since LR is fast
    use_rf_tt  = False, # Use random forest regression. Defaults to False since RF is slow
    use_mlp_tt = False, # Use multi-layer perceptron regression. Defaults to False since MLP is slow
    rf_n_estimators = 100,
    rf_max_depth = 5,
    rf_n_jobs = 3,
    mlp_max_iter = 200,
    mlp_hidden_layer_sizes = 100
    ) :
    """
    """
    
    time_start = datetime.now()
    
    # Tries to makes the output directory if it doesn't exist yet
    if ( not path.exists(output_directory) ): mkdir(output_directory)
    output_joblib_directory = output_directory + "JOBLIB_Backup/"
    if ( not path.exists(output_joblib_directory) ): mkdir(output_joblib_directory)

    
    # Loops over each training file
    for train_file_info in train_file_bundle:
        train_csv_path = train_file_info[0]
        X_train, y_train, sc_train = Build_FeatureArrays_FromCSV(train_csv_path)
        train_bias  = train_file_info[2]
        # Makes train bias directory
        output_train_directory = output_directory + "Train_" + train_bias + "/"
        if ( not path.exists(output_train_directory) ): mkdir(output_train_directory)
        
        # Loops over each feature set to train and test with
        for feature_set in feature_bundle:
            feature_label = feature_set[0]
            feature_label_lr = feature_set[0].copy()
            feature_label_lr.append("lr_intercept") # Adds field for linear regression y-intercept
            feature_index = feature_set[1]
            
            # Makes feature set directory
            output_feature_directory = output_train_directory
            if ( not path.exists(output_feature_directory) ): mkdir(output_feature_directory)
            
            #############################################################
            #                                                           #
            #   TRAIN ML OVER FULL GeV RANGE, TEST OVER SPECIFIC BINS   #
            #                                                           #
            #############################################################
            
            # Builds training and testing arrays
            print("\nBuilding training and testing selected feature arrays...")
            X_train_select  = Build_SelectFeatureArray(X_train, feature_index)
            lr_pipeline     = None
            lr_coeffs       = None
            rf_pipeline     = None
            rf_coeffs       = None
            mlp_pipeline    = None
            output_csv_base = "Train_" + train_bias + "_F" + str(len(feature_label)) + "_" + str(int(train_pt_min)) + "_" + str(int(train_pt_max))
            
            # Trains estimators
            if use_lr:
                print("\nTraining linear regression estimator...")
                pipeline_backup = output_joblib_directory + output_csv_base + "_Pipeline_LR.joblib"
                coeff_backup = output_joblib_directory + output_csv_base + "_Pipeline_LR_Coeffs.joblib"
                coeff_csv = output_feature_directory + output_csv_base + "_LR_Coeffs.csv"
                if ( path.exists(pipeline_backup) and path.exists(coeff_csv) ):
                    print("Loading", pipeline_backup)
                    print("Loading", coeff_backup)
                    lr_pipeline = load(pipeline_backup)
                    lr_coeffs = load(coeff_backup)
                else:
                    lr_pipeline, lr_coeffs = Train_LinearRegression(
                        X_train_select,
                        y_train,
                        feature_label_lr,
                        use_scaler = True)
                    Write_MLCoefficients_ToCSV(
                        coeff_csv,
                        lr_coeffs,
                        feature_label_lr
                        )
                    dump(lr_pipeline, pipeline_backup)
                    dump(lr_coeffs, coeff_backup)
                    
            if use_rf:
                print("\nTraining random forest regression estimator...")
                pipeline_backup = output_joblib_directory + output_csv_base + "_Pipeline_RF.joblib"
                coeff_backup = output_joblib_directory + output_csv_base + "_Pipeline_RF_Coeffs.joblib"
                coeff_csv = output_feature_directory + output_csv_base + "_RF_Coeffs.csv"
                if ( path.exists(pipeline_backup) and path.exists(coeff_backup) ):
                    print("Loading", pipeline_backup)
                    print("Loading", coeff_backup)
                    rf_pipeline = load(pipeline_backup)
                    rf_coeffs = load(coeff_backup)
                    Write_MLCoefficients_ToCSV(
                        coeff_csv,
                        rf_coeffs,
                        feature_label
                        )
                else:
                    rf_pipeline, rf_coeffs = Train_RandomForestRegression(
                        X_train_select,
                        y_train,
                        feature_label,
                        use_scaler = True,
                        n_estimators = rf_n_estimators,
                        max_depth = rf_max_depth,
                        n_jobs = rf_n_jobs
                        )
                    Write_MLCoefficients_ToCSV(
                        coeff_csv,
                        rf_coeffs,
                        feature_label
                        )
                    dump(rf_pipeline, pipeline_backup)
                    dump(rf_coeffs, coeff_backup)
            
            if use_mlp:
                print("\nTraining multilayer perceptron (neural net) regression estimator...")
                pipeline_backup = output_joblib_directory + output_csv_base + "_Pipeline_MLP.joblib"
                if ( path.exists(pipeline_backup) ):
                    print("Loading", pipeline_backup)
                    print("Loading", coeff_backup)
                    mlp_pipeline = load(pipeline_backup)
                else:
                    mlp_pipeline = Train_MLPRegression(
                        X_train_select,
                        y_train,
                        feature_label,
                        use_scaler = True,
                        max_iter = mlp_max_iter,
                        hidden_layer_sizes = mlp_hidden_layer_sizes
                        )
                    dump(mlp_pipeline, pipeline_backup)
            
            # Loops over each test file
            for test_file_info in test_file_bundle:
                test_csv_path = test_file_info[0]
                test_bias     = test_file_info[2]
                X_test,  y_test,  sc_test = Build_FeatureArrays_FromCSV(test_csv_path)
                X_test_select = Build_SelectFeatureArray(X_test, feature_index)
                
                # Makes test bias directory
                output_test_directory = output_feature_directory + "Test_" + test_bias + "/"
                if ( not path.exists(output_test_directory) ): mkdir(output_test_directory)
                
                for test_bin_set in test_bin_array:
                    test_bin_label = test_bin_set[0]
                    
                    # Builds Train Bin outputs directories
                    output_directory_temp = output_test_directory + "F" + str(len(feature_label)) + "_" + test_bin_label + "/"
                    if ( not path.exists(output_directory_temp) ): mkdir(output_directory_temp)
                            
                    for test_bin in test_bin_set[1]:
                        test_pt_min = test_bin[0]
                        test_pt_max = test_bin[1]

                        output = "\nTesting " + str(len(feature_index)) + " features on " + str(test_pt_min) + "-" + str(test_pt_max) + " GeV..."
                        print(output)
                        
                        output_csv_path = output_directory_temp + output_csv_base + "_Test_" + str(int(test_pt_min)) + "_" + str(int(test_pt_max))
                        
                        if use_lr:
                            TestAndSave_Estimator(
                                feature_label,
                                feature_index,
                                lr_pipeline,
                                X_test_select,
                                y_test,
                                sc_test,
                                test_pt_min,
                                test_pt_max,
                                output_csv_path + "_LR.csv",
                                use_scaler = True
                                )
                        if use_rf:
                            TestAndSave_Estimator(
                                feature_label,
                                feature_index,
                                rf_pipeline,
                                X_test_select,
                                y_test,
                                sc_test,
                                test_pt_min,
                                test_pt_max,
                                output_csv_path + "_RF.csv",
                                use_scaler = True
                                )
                        if use_mlp:
                            TestAndSave_Estimator(
                                feature_label,
                                feature_index,
                                mlp_pipeline,
                                X_test_select,
                                y_test,
                                sc_test,
                                test_pt_min,
                                test_pt_max,
                                output_csv_path + "_MLP.csv",
                                use_scaler = True
                                )

                        print("Test and save complete!\n")
            
            #################################################
            #                                               #
            #   TRAIN AND TEST OVER SMALLER BINS USING LR   #
            #                                               #
            #################################################
            
            if traintest_bin_array :
                for traintest_bin_set in traintest_bin_array:
                    traintest_bin_label = traintest_bin_set[0]
                    
                    # Builds output directory
                    output_directory_temp = output_test_directory + "F" + str(len(feature_label)) + "_" + traintest_bin_label + "/"
                    if ( not path.exists(output_directory_temp) ): mkdir(output_directory_temp)
                    
                    for traintest_bin in traintest_bin_set[1]:
                        train_pt_min_temp = traintest_bin[0]
                        train_pt_max_temp = traintest_bin[1]
                        test_pt_min = traintest_bin[0]
                        test_pt_max = traintest_bin[1]

                        X_train_cut  = []
                        y_train_cut  = []
                        sc_train_cut = []
                        for i in range(len(X_train)):
                            if (y_train[i] > train_pt_min_temp) and (y_train[i] < train_pt_max_temp):
                                X_train_cut.append(X_train[i])
                                y_train_cut.append(y_train[i])
                                sc_train_cut.append(sc_train[i])

                        X_test_cut  = []
                        y_test_cut  = []
                        sc_test_cut = []
                        for i in range(len(X_test)):
                            if (y_test[i] > test_pt_min) and (y_test[i] < test_pt_max):
                                X_test_cut.append(X_test[i])
                                y_test_cut.append(y_test[i])
                                sc_test_cut.append(sc_test[i])

                        # Builds training and testing arrays
                        print("\nBuilding training and testing selected feature arrays...")
                        X_train_select = Build_SelectFeatureArray(X_train_cut, feature_index)
                        X_test_select  = Build_SelectFeatureArray(X_test_cut, feature_index)

                        # Trains estimator
                        output_joblib_train = output_joblib_directory + "Train_" + train_bias + "_F" + str(len(feature_label)) + "_" + str(int(train_pt_min_temp)) + "_" + str(int(train_pt_max_temp))
                        output_csv_train = output_directory_temp + "Train_" + train_bias + "_F" + str(len(feature_label)) + "_" + str(int(train_pt_min_temp)) + "_" + str(int(train_pt_max_temp))
                        output_csv_test = output_csv_train + "_Test_" + str(int(test_pt_min)) + "_" + str(int(test_pt_max))
                        
                        if ( use_lr_tt ):
                            print("\nTraining and testing linear regression estimator...")
                            pipeline_backup = output_joblib_train + "_Pipeline_LR.joblib"
                            coeff_backup = output_joblib_train + "_Pipeline_LR_Coeffs.joblib"
                            lr_pipeline = None
                            lr_coeffs = None
                            if ( path.exists(pipeline_backup) and path.exists(coeff_backup) ):
                                print("Loading", pipeline_backup)
                                print("Loading", coeff_backup)
                                lr_pipeline = load(pipeline_backup)
                                lr_coeffs = load(coeff_backup)
                            else:
                                lr_pipeline, lr_coeffs = Train_LinearRegression(
                                    X_train_select,
                                    y_train_cut,
                                    feature_label_lr,
                                    use_scaler = True)
                                Write_MLCoefficients_ToCSV(
                                    output_csv_train + "_LR_Coeffs.csv",
                                    lr_coeffs,
                                    feature_label_lr
                                    )
                                dump(lr_pipeline, pipeline_backup)
                                dump(lr_coeffs, coeff_backup)
                            TestAndSave_Estimator(
                                feature_label,
                                feature_index,
                                lr_pipeline,
                                X_test_select,
                                y_test_cut,
                                sc_test_cut,
                                test_pt_min,
                                test_pt_max,
                                output_csv_test + "_LR.csv",
                                use_scaler = True
                                )
                        if ( use_rf_tt ):
                            print("\nTraining and testing random forest estimator...")
                            pipeline_backup = output_joblib_train + "_Pipeline_RF.joblib"
                            coeff_backup = output_joblib_train + "_Pipeline_RF_Coeffs.joblib"
                            rf_pipeline = None
                            rf_coeffs = None
                            if ( path.exists(pipeline_backup) and path.exists(coeff_backup) ):
                                print("Loading", pipeline_backup)
                                print("Loading", coeff_backup)
                                rf_pipeline = load(pipeline_backup)
                                rf_coeffs = load(coeff_backup)
                            else:
                                rf_pipeline, rf_coeffs = Train_RandomForestRegression(
                                    X_train_select,
                                    y_train_cut,
                                    feature_label,
                                    use_scaler = True,
                                    n_estimators = rf_n_estimators,
                                    max_depth = rf_max_depth,
                                    n_jobs = rf_n_jobs
                                    )
                                Write_MLCoefficients_ToCSV(
                                    output_csv_train + "_RF_Coeffs.csv",
                                    rf_coeffs,
                                    feature_label
                                    )
                                dump(lr_pipeline, pipeline_backup)
                                dump(lr_coeffs, coeff_backup)
                            TestAndSave_Estimator(
                                feature_label,
                                feature_index,
                                rf_pipeline,
                                X_test_select,
                                y_test_cut,
                                sc_test_cut,
                                test_pt_min,
                                test_pt_max,
                                output_csv_test + "_RF.csv",
                                use_scaler = True
                                )
                        if ( use_mlp_tt ): # Can execute with MLP
                            print("\nTraining and testing MLP estimator...")
                            pipeline_backup = output_joblib_train + "_Pipeline_MLP.joblib"
                            mlp_pipeline = None
                            if ( path.exists(pipeline_backup) ):
                                print("Loading", pipeline_backup)
                                mlp_pipeline = load(pipeline_backup)
                            else:
                                mlp_pipeline = Train_MLPRegression(
                                    X_train_select,
                                    y_train_cut,
                                    feature_label,
                                    use_scaler = True,
                                    max_iter = mlp_max_iter,
                                    hidden_layer_sizes = mlp_hidden_layer_sizes
                                    )
                                dump(mlp_pipeline, pipeline_backup)
                            TestAndSave_Estimator(
                                feature_label,
                                feature_index,
                                mlp_pipeline,
                                X_test_select,
                                y_test_cut,
                                sc_test_cut,
                                test_pt_min,
                                test_pt_max,
                                output_csv_test + "_MLP.csv",
                                use_scaler = True
                                )

                        print("Test and save complete!\n")
                        
    # Prints timestamp of completion
    time_end = datetime.now()
    time_delta = time_end - time_start
    time_end_string = time_end.strftime("%Y/%m/%d %H:%M:%S")
    print("\nComplete!", time_end_string)
    print("Process duration:", time_delta)
