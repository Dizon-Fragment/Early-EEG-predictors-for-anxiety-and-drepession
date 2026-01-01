%% FINAL Interaction & Cross-Network Analysis Pipeline (v8 - Explicit Targets & Betas)
% Fixed pathing issue where 'output_dir' was not updating inside the loop.
clear; close all; clc;

%% 1. Parameters
main_path = '/data/home/EEG001/Longitudinal/datafromMac/'; 
result_plot = '/data/home/EEG001/Longitudinal/ForPlot/';
symptom_path = fullfile(main_path, 'longitudinal_prefinish', 'subjset_symptom.mat');

folders = {'csd', 'norm'};
p_thresh = [0.05, 0.01];
year_data = {'2015', '2017', '2019'};
idx = {'coh', 'plv'};

% Ensure opname is a char array (string), not a cell array, to match Summary Script
opname_base = 'preLooPERM2idx_cb_v5_mse'; 
opname = [opname_base, '_edge']; 

no_iterations = 1000; 

%% 2. Main Analysis Loop
% We combine directory setup and analysis into one loop to ensure paths are correct.
for dd = 1:length(folders)
    
    % --- Define Output Directory for current folder ---
    % FIX: output_dir must be defined inside the loop for 'dd'
    output_dir = fullfile(result_plot, 'result', folders{dd}, opname);
    
    % Create directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
        fprintf('Created directory: %s\n', output_dir);
    end

    for yy = 1:length(year_data)
        
        fprintf('\nProcessing: %s - %s...\n', folders{dd}, year_data{yy});
        
        % --- Load Data ---
        input_dir = fullfile(main_path, 'result', folders{dd}, 'CPM');
        
        % Check if input directory exists to prevent crash
        if ~exist(input_dir, 'dir')
            error('Input directory not found: %s', input_dir);
        end
        
        load(symptom_path);
        
        % Define Targets (Y)
        % Ensure variable names match your loaded .mat content
        Y_Anxiety = table2array(symptom_score(:, 'SAS'));
        Y_Depression = table2array(symptom_score(:, 'SDS'));
        Y_Comorbidity = Y_Anxiety + Y_Depression; % Sum score
        
        Targets = struct('Anxiety', Y_Anxiety, 'Depression', Y_Depression, 'Comorbidity', Y_Comorbidity);
        target_names = fieldnames(Targets);
        
        % Load Matrices
        alpha_mats = {}; beta1_mats = {};
        for ii = 1:length(idx)
            % Construct filenames
            f_alpha = fullfile(input_dir, [idx{ii},'_alpha_',year_data{yy},'_ForCPM.mat']);
            f_beta  = fullfile(input_dir, [idx{ii},'_beta1_',year_data{yy},'_ForCPM.mat']);
            
            % Load and assign
            load(f_alpha);
            alpha_mats{ii} = eval(['sub_',idx{ii},'_alpha_',year_data{yy}]);
            
            load(f_beta);
            beta1_mats{ii} = eval(['sub_',idx{ii},'_beta1_',year_data{yy}]);
        end
        clear sub_*; % Clear temporary variables to save memory
        
        % =================================================================
        %  STEP A: PRE-COMPUTE STAGE 1 (Feature Extraction)
        %  Extract Brain Scores for Alpha_Anx and Beta_Dep
        % =================================================================
        fprintf('  > Running Stage 1 (Feature Extraction)...\n');
        
        % Hypothesis 1: Alpha predicts Anxiety
        s1_Alpha_Anx_p05 = run_single_cpm(alpha_mats, Y_Anxiety, 0.05, idx);
        s1_Alpha_Anx_p01 = run_single_cpm(alpha_mats, Y_Anxiety, 0.01, idx);
        
        % Hypothesis 2: Beta predicts Depression
        s1_Beta_Dep_p05 = run_single_cpm(beta1_mats, Y_Depression, 0.05, idx);
        s1_Beta_Dep_p01 = run_single_cpm(beta1_mats, Y_Depression, 0.01, idx);
        
        % =================================================================
        %  STEP B: RUN STAGE 2, 3, 4 FOR EACH TARGET VARIABLE
        %  Predict Anx, Dep, and Comorbidity using Brain Scores
        % =================================================================
        Extract_Results = struct();
        
        % Save Stage 1 Base info for reference
        Extract_Results.Stage1.Alpha_Anx_p05 = s1_Alpha_Anx_p05;
        Extract_Results.Stage1.Beta_Dep_p05  = s1_Beta_Dep_p05;
        Extract_Results.Stage1.Alpha_Anx_p01 = s1_Alpha_Anx_p01;
        Extract_Results.Stage1.Beta_Dep_p01  = s1_Beta_Dep_p01;

        for t = 1:length(target_names)
            t_name = target_names{t};
            Y_current = Targets.(t_name);
            
            fprintf('    > Analyzing Target: %s\n', t_name);
            
            % --- Stage 2: Interaction Test (Force Interaction) ---
            % Input 1: Alpha_Anx_Score, Input 2: Beta_Dep_Score -> Predict Y_current
            Extract_Results.Stage2.(t_name).p05_vs_p05 = ...
                test_interaction_effects(s1_Alpha_Anx_p05, s1_Beta_Dep_p05, Y_current, no_iterations);
            
            % --- Stage 3: Model Comparison (Best Model) ---
            Extract_Results.Stage3.(t_name).p05_vs_p05 = ...
                find_best_predictive_model(s1_Alpha_Anx_p05, s1_Beta_Dep_p05, Y_current, no_iterations);
            
            % --- Stage 4: Cross-Threshold (Robustness) ---
            Extract_Results.Stage4.(t_name).p05_vs_p01 = ...
                find_best_predictive_model(s1_Alpha_Anx_p05, s1_Beta_Dep_p01, Y_current, no_iterations);
                
            Extract_Results.Stage4.(t_name).p01_vs_p05 = ...
                find_best_predictive_model(s1_Alpha_Anx_p01, s1_Beta_Dep_p05, Y_current, no_iterations);
        end
        
        % Save Results
        % FIX: Now uses the correctly updated 'output_dir'
        output_file_name = [folders{dd}, '_', year_data{yy}, '_', opname, '.mat'];
        output_file_path = fullfile(output_dir, output_file_name);
        
        save(output_file_path, 'Extract_Results');
        fprintf('  > Saved: %s\n', output_file_path);
    end
end
disp('Analysis Complete.');