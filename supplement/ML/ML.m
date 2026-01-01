%% FINAL AUTOMATED SEARCH & PERMUTATION
clear; close all; clc;

% --- 1. SETTINGS ---
main_path = '/data/home/EEG001/Longitudinal/datafromMac/'; 
symptom_path = fullfile( main_path, 'longitudinal_prefinish', 'subjset_symptom.mat');
input_base_dir = fullfile(main_path, 'result', 'csd', 'CPM'); 

idx = {'coh', 'plv'};
target_years = {'2015', '2017', '2019'};
n_perms = 1000; % Permutation iterations

% Define Pairs
target_pairs = {
    'alpha', 'SAS';
    'beta1', 'SDS'
};

% Define Parameter Grid
model_grid = {
    'ols',          'None',         [0];                
    'ridge',        'Lambda',       [0.01, 0.1, 1, 10]; 
    'lasso',        'Lambda',       [0.001, 0.01, 0.1]; 
    'svm',          'BoxConstraint',[0.01, 0.1, 1, 10]; 
    'randomforest', 'MinLeafSize',  [3, 5, 10, 15]      
};

p_thresh_list = [0.05, 0.01, 0.005, 0.001]; 

% --- 2. LOAD DATA ---
fprintf('Loading behavior...\n');
load(symptom_path);
T_name = symptom_score.Properties.VariableNames;
T_result = table2array(symptom_score);

% Storage for final valid results
% Columns: {Year, Freq, Behav, P_Thr, Model, Param, Val, LOOCV_r, LOOCV_p, PERM_p}
final_results = {}; 

fprintf('\n========================================================================================================================\n');
fprintf('%-5s | %-5s | %-5s | %-5s | %-12s | %-8s | %-6s | %-7s | %-7s | %-7s\n', ...
    'Year', 'Freq', 'Behav', 'P_Thr', 'Model', 'Param', 'Val', 'LOOCV_r', 'LOOCV_p', 'PERM_p');
fprintf('========================================================================================================================\n');

% --- 3. MAIN LOOP ---
for pair_idx = 1:size(target_pairs, 1)
    ff = target_pairs{pair_idx, 1};
    beh_name_partial = target_pairs{pair_idx, 2};
    
    % Find behavior column
    beh_col = find(contains(T_name, beh_name_partial));
    if isempty(beh_col), continue; end
    beh_col = beh_col(1); 
    real_beh_name = T_name{beh_col};
    
    % Pre-load Data
    data_cache = struct();
    for y = 1:length(target_years)
        yy = target_years{y};
        try
            for ii = 1:length(idx)
                 d = load(fullfile(input_base_dir, [idx{ii},'_', ff, '_', yy,'_','ForCPM.mat']));
                 data_cache.(['y',yy]){ii} = d.(['sub_', idx{ii}, '_', ff, '_', yy]);
            end
        catch
            data_cache.(['y',yy]) = [];
        end
    end
    
    for y = 1:length(target_years)
        yy = target_years{y};
        mats = data_cache.(['y',yy]);
        if isempty(mats), continue; end
        
        behav = T_result(:, beh_col);
        
        % Remove NaNs
        nan_idx = isnan(behav);
        if any(nan_idx)
            behav(nan_idx) = [];
            for ii = 1:length(mats), mats{ii}(:,:,nan_idx) = []; end
        end
        
        % Loop Thresholds
        for pp = p_thresh_list
            % Loop Models
            for m_idx = 1:size(model_grid, 1)
                mdl_name = model_grid{m_idx, 1};
                param_name = model_grid{m_idx, 2};
                param_vals = model_grid{m_idx, 3};
                
                % Loop Hyperparameters
                for val = param_vals
                    full_mdl_name = [mdl_name, '_all'];
                    if strcmpi(mdl_name, 'ols'), pass_val = []; else, pass_val = val; end
                    
                    % 1. RUN LOOCV
                    [pred, ~] = predict_with_loocv(mats, behav, full_mdl_name, pp, pass_val);
                    [r_true, p_loocv] = corr(pred, behav);
                    
                    % 2. DECISION: If Significant, run PERMUTATION
                    perm_p = NaN; % Default if not run
                    sig_mark = '';
                    
                    if p_loocv < 0.05
                        % >>> START PERMUTATION <<<
                        % fprintf('    [Sig detected] Permuting %s %s...\n', mdl_name, num2str(val));
                        perm_rs = zeros(n_perms, 1);
                        n_sub = length(behav);
                        
                        % Use parfor if possible
                        parfor k = 1:n_perms
                            shuf_b = behav(randperm(n_sub));
                            [pred_p, ~] = predict_with_loocv(mats, shuf_b, full_mdl_name, pp, pass_val);
                            [r_p, ~] = corr(pred_p, shuf_b);
                            if isnan(r_p), r_p=0; end
                            perm_rs(k) = r_p;
                        end
                        
                        perm_p = (sum(perm_rs >= r_true) + 1) / (n_perms + 1);
                        
                        if perm_p < 0.05
                            sig_mark = '*'; 
                            % Store valid result
                            final_results(end+1,:) = {yy, ff, beh_name_partial, pp, mdl_name, param_name, val, r_true, p_loocv, perm_p};
                        end
                    end
                    
                    % Only print if LOOCV was significant (to save screen space)
                    if p_loocv < 0.05
                        fprintf('%-5s | %-5s | %-5s | %.3f | %-12s | %-8s | %-6.3f | %-7.3f | %-7.3f | %.4f %s\n', ...
                            yy, ff, beh_name_partial, pp, mdl_name, param_name, val, r_true, p_loocv, perm_p, sig_mark);
                    end
                    
                end
            end
        end
    end
end

fprintf('\nDone. "final_results" variable contains all robust models.\n');
save('All_Permutation_Significant_Results.mat', 'final_results');


% =========================================================================
% FUNCTION: predict_with_loocv (Updated)
% =========================================================================
function [behav_pred, all_models] = predict_with_loocv(all_mats, all_behav, model_name, p_thresh, hyperparam)
    if nargin < 5, hyperparam = []; end

    n_subjects = size(all_mats{1}, 3);
    n_nodes = size(all_mats{1}, 1);
    n_modalities = numel(all_mats);
    
    behav_pred = zeros(n_subjects, 1);
    all_models = cell(n_subjects, 1);
    
    for i = 1:n_subjects 
        train_indices = true(1, n_subjects);
        train_indices(i) = false;
        
        test_mat = cell(n_modalities, 1);
        train_mats = cell(n_modalities, 1);
        for m = 1:n_modalities
            train_mats{m} = all_mats{m}(:, :, train_indices);
            test_mat{m} = all_mats{m}(:, :, i);
        end
        train_behav = all_behav(train_indices);

        pos_masks = cell(n_modalities, 1);
        neg_masks = cell(n_modalities, 1);
        for m = 1:n_modalities
            train_vcts = reshape(train_mats{m}, n_nodes*n_nodes, []);
            [r, p] = corr(train_vcts', train_behav, 'rows', 'complete');
            r_mat = reshape(r, n_nodes, n_nodes);
            p_mat = reshape(p, n_nodes, n_nodes);
            pos_masks{m} = (r_mat > 0) & (p_mat < p_thresh);
            neg_masks{m} = (r_mat < 0) & (p_mat < p_thresh);
        end
        
        n_train = n_subjects - 1;
        X_train_full = zeros(n_train, n_modalities * 2);
        X_test_full = zeros(1, n_modalities * 2);

        for m = 1:n_modalities
            for s = 1:n_train
                mat = train_mats{m}(:,:,s);
                X_train_full(s, (m-1)*2+1) = sum(mat(pos_masks{m}), 'all');
                X_train_full(s, (m-1)*2+2) = sum(mat(neg_masks{m}), 'all');
            end
            X_test_full(1, (m-1)*2+1) = sum(test_mat{m}(pos_masks{m}), 'all');
            X_test_full(1, (m-1)*2+2) = sum(test_mat{m}(neg_masks{m}), 'all');
        end
        
        if endsWith(model_name, '_pos')
            X_train_sel = X_train_full(:, 1:2:end);
            X_test_sel = X_test_full(:, 1:2:end);
            base_model = extractBefore(model_name, '_pos');
        elseif endsWith(model_name, '_neg')
            X_train_sel = X_train_full(:, 2:2:end);
            X_test_sel = X_test_full(:, 2:2:end);
            base_model = extractBefore(model_name, '_neg');
        elseif endsWith(model_name, '_all')
            X_train_sel = X_train_full;
            X_test_sel = X_test_full;
            base_model = extractBefore(model_name, '_all');
        else
            X_train_sel = X_train_full;
            X_test_sel = X_test_full;
            base_model = model_name;
        end

        [y_pred, model] = train_and_predict(base_model, X_train_sel, train_behav, X_test_sel, hyperparam);
        behav_pred(i) = y_pred;
        all_models{i} = model;
    end
end

% =========================================================================
% FUNCTION: train_and_predict (Updated)
% =========================================================================
function [y_pred, model] = train_and_predict(model_name, X_train, y_train, X_test, hyperparam)
    if any(isnan(X_train(:))), X_train(isnan(X_train)) = 0; end
    if any(isnan(X_test(:))), X_test(isnan(X_test)) = 0; end
    
    if ~strcmpi(model_name, 'ols')
        mu = mean(X_train);
        sigma = std(X_train);
        sigma(sigma==0) = 1; 
        X_train = (X_train - mu) ./ sigma;
        X_test  = (X_test - mu) ./ sigma;
    end

    model = []; y_pred = NaN;

    switch lower(model_name)
        case 'ols'
            X_train_bias = [ones(size(X_train,1),1), X_train];
            X_test_bias  = [ones(size(X_test,1),1), X_test];
            b = regress(y_train, X_train_bias);
            y_pred = X_test_bias * b;

        case 'ridge'
            lam = 0.1; 
            if ~isempty(hyperparam), lam = hyperparam; end
            model = fitrlinear(X_train, y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lam);
            y_pred = predict(model, X_test);

        case 'lasso'
            if ~isempty(hyperparam)
                [B, FitInfo] = lasso(X_train, y_train, 'Alpha', 1, 'Lambda', hyperparam);
                b = B(:,1); intercept = FitInfo.Intercept(1);
            else
                [B, FitInfo] = lasso(X_train, y_train, 'Alpha', 1, 'NumLambda', 20);
                [~, idx] = min(FitInfo.MSE);
                if all(B(:, idx) == 0), non_zeros = find(any(B~=0,1)); if ~isempty(non_zeros), idx=non_zeros(end); else, idx=length(FitInfo.Lambda); end; end
                b = B(:, idx); intercept = FitInfo.Intercept(idx);
            end
            y_pred = X_test * b + intercept;

        case 'svm'
            C = 1;
            if ~isempty(hyperparam), C = hyperparam; end
            model = fitrsvm(X_train, y_train, 'KernelFunction', 'linear', 'Standardize', false, 'BoxConstraint', C);
            y_pred = predict(model, X_test);
            
        case 'randomforest'
            leaf = 5;
            if ~isempty(hyperparam), leaf = hyperparam; end
            model = TreeBagger(50, X_train, y_train, 'Method', 'regression', 'OOBPrediction', 'on', 'MinLeafSize', leaf, 'NumPredictorsToSample', 'all');
            y_pred = predict(model, X_test);

        otherwise
            error('Unknown model: %s', model_name);
    end
end