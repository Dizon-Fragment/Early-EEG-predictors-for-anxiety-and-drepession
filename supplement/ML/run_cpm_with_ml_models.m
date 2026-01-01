% =========================================================================
% FILE: run_cpm_with_ml_models.m (FINAL VERSION - Implements Grouped Testing Logic)
% =========================================================================
function all_results = run_cpm_with_ml_models(all_mats, all_behav, p_thresh, n_permutations)
% This final version implements a two-phase, grouped testing logic.
% For each base model (e.g., ridge), it first runs LOOCV for all variants (_all, _pos, _neg).
% If ANY variant is initially significant, it then runs permutation tests for ALL variants.
% All initial LOOCV results (metrics, predictions, masks) are always saved.

    % --- Parameter Definition ---
    base_models = {'ols','ridge', 'lasso', 'elasticnet', 'svm', 'randomforest', 'gp'};
    feature_variants = {'all', 'pos', 'neg'}; % Suffixes are added later
    permutation_alpha = 0.05;
    all_results = struct();

    fprintf('\n======= Starting Evaluation with Grouped Testing Logic =======\n');

    % --- Outer loop through each BASE MODEL (e.g., ridge, lasso) ---
    for i = 1:length(base_models)
        base_name = base_models{i};
        fprintf('\n\n--- Processing Model Family: %s ---\n', upper(base_name));
        
        % --- PHASE 1: Initial LOOCV Run for all variants in the family ---
        fprintf('  Phase 1: Running initial LOOCV for all variants (_all, _pos, _neg)...\n');
        temp_results = struct();
        run_permutations_for_this_family = false;

        for j = 1:length(feature_variants)
            variant_name = feature_variants{j};
            model_name_full = [base_name, '_', variant_name];
            
            % Run LOOCV and get predictions, models, and masks
            [behav_pred, ~, pos_masks, neg_masks] = predict_with_loocv(all_mats, all_behav, model_name_full, p_thresh);
            
            % Calculate and store metrics
            metrics = calculate_performance_metrics(all_behav, behav_pred);
            temp_results.(variant_name).metrics = metrics;
            temp_results.(variant_name).predicted_behavior = behav_pred;

            % Save masks (only needs to be done once per feature set)
            if ~isfield(all_results, 'edge_masks') || ~isfield(all_results.edge_masks, variant_name)
                 all_results.edge_masks.(variant_name).pos_masks_per_fold = pos_masks;
                 all_results.edge_masks.(variant_name).neg_masks_per_fold = neg_masks;
            end
            
            fprintf('    - Variant [%s]: Initial r=%.3f, p=%.4f\n', upper(variant_name), metrics.r, metrics.p);
            
            % Check if this variant's performance warrants permutation for the whole family
            if metrics.p < permutation_alpha
                run_permutations_for_this_family = true;
            end
        end

        % --- DECISION POINT ---
        if run_permutations_for_this_family
            fprintf('  Decision: At least one variant was significant. Proceeding to permutation tests for the ENTIRE family.\n');
        else
            fprintf('  Decision: No variant was significant. Skipping permutation tests for this family.\n');
        end

        % --- PHASE 2: Permutation Run (if warranted) ---
        for j = 1:length(feature_variants)
            variant_name = feature_variants{j};
            
            if run_permutations_for_this_family
                fprintf('  Phase 2: Running permutations for variant [%s]...\n', upper(variant_name));
                model_name_full = [base_name, '_', variant_name];
                
                perm_dist = zeros(n_permutations, 1);
                n_subjects = length(all_behav);
                
                parfor p_iter = 1:n_permutations
                    shuffled_behav = all_behav(randperm(n_subjects));
                    [behav_pred_perm, ~] = predict_with_loocv(all_mats, shuffled_behav, model_name_full, p_thresh);
                    [perm_r, ~] = corr(behav_pred_perm, shuffled_behav, 'rows', 'complete');
                    if isnan(perm_r), perm_r = 0; end
                    perm_dist(p_iter) = perm_r;
                end
                
                true_r = temp_results.(variant_name).metrics.r;
                p_value = (sum(abs(perm_dist) >= abs(true_r)) + 1) / (n_permutations + 1);
                
                temp_results.(variant_name).permutation_p_value = p_value;
                temp_results.(variant_name).permutation_distribution = perm_dist;
                fprintf('    - Permutation test finished. P-value = %.4f\n', p_value);
            else
                % If skipping, fill in with NaN values
                temp_results.(variant_name).permutation_p_value = NaN;
                temp_results.(variant_name).permutation_distribution = [];
            end
        end
        
        % --- Final Step: Copy temporary results to the main results struct ---
        all_results.(base_name) = temp_results;
    end
    fprintf('\n======= All Model Evaluations Completed =======\n');
end