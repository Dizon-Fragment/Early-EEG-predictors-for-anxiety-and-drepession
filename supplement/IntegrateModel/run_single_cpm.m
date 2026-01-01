% --- Save this file as "run_single_cpm.m" ---
function [results] = run_single_cpm(all_mats, all_behav, p_thresh, idx)
    % Runs a standard LOOCV CPM and returns detailed results for pos, neg, and glm networks.
    
    [~,~,~,behav_pred_pos, behav_pred_neg, behav_pred_glm, pos_matri,...
     neg_matri, ~, ~, ~, ~, P_pos, P_neg, P_glm, fit_pos, fit_neg, fit_glm, ...
     pos_matri_sum, neg_matri_sum] = ...
     predict_behavior_longitudinal_forplot(all_mats, all_behav, p_thresh, idx); % Assuming you use this function

    % Package results for each network type
    results.pos = package_s1(behav_pred_pos, all_behav, P_pos, fit_pos, pos_matri_sum, pos_matri);
    results.neg = package_s1(behav_pred_neg, all_behav, P_neg, fit_neg, neg_matri_sum, neg_matri);
    results.glm = package_s1(behav_pred_glm, all_behav, P_glm, fit_glm, pos_matri_sum + neg_matri_sum, cat(4, pos_matri, neg_matri));
end

function s = package_s1(pred, actual, p_pearson, fit, matri_sum, matri_fit)
    s = struct();
    s.behav_pred = pred;
    [r, ~] = corr(pred, actual, 'rows', 'complete');
    s.r = r;
    s.p_pearson = p_pearson; % The p-value from the simple correlation
    s.fit_coef = fit; % Regression coefficients from each fold
    s.mask_sum = matri_sum; % Sum of masks across all folds
    s.mask_fit = matri_fit; % Masks used in each fold
end