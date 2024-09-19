function [true_prediction_r_pos,true_prediction_r_neg,true_prediction_r_glm,behav_pred_pos, behav_pred_neg, behav_pred_glm, pos_matri,...
    neg_matri, true_MSE_pos, true_MSE_neg, true_MSE_glm, error_run, P_pos, P_neg,P_glm,fit_pos,fit_neg,fit_glm,pos_matri_sum,neg_matri_sum] = predict_behavior_longitudinal_forplot(all_mats,all_behav, thresh, idx)
% % Copyright 2024 Guangzhi Deng

% This code is rewrited by Guangzhi Deng for Longitudinal CPM.

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as:

% Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang, Chun MM,
% Papademetris X & Constable RT. (2015). Functional connectome
% fingerprinting: Identifying individuals using patterns of brain
% connectivity. Nature Neuroscience 18, 1664-1671.

% This code provides a framework for implementing functional
% connectivity-based behavioral prediction in a leave-one-subject-out
% cross-validation scheme, as described in Finn, Shen et al 2015 (see above
% for full reference). The first input ('all_mats') is a pre-calculated
% MxMxN matrix containing all individual-subject connectivity matrices,
% where M = number of nodes in the chosen brain atlas and N = number of
% subjects. Each element (i,j,k) in these matrices represents the
% correlation between the BOLD timecourses of nodes i and j in subject k
% during a single fMRI session. The second input ('all_behav') is the
% Nx1 vector of scores for the behavior of interest for all subjects.

% As in the reference paper, the predictive power of the model is assessed
% via correlation between predicted and observed scores across all
% subjects. Note that this assumes normal or near-normal distributions for
% both vectors, and does not assess absolute accuracy of predictions (only
% relative accuracy within the sample). It is recommended to explore
% additional/alternative metrics for assessing predictive power, such as
% prediction error sum of squares or prediction r^2.


% ------------ INPUTS -------------------

warning('off');
% ---------------------------------------

no_sub = size(all_mats{1},3);
no_node = size(all_mats{1},1);

behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);
behav_pred_glm = zeros(no_sub,1);
error_run = {};

for leftout = 1:no_sub
    for ii = 1:numel(all_mats)
        
        
        train_mats = all_mats{ii};
        train_mats(:,:,leftout) = [];
        train_vcts = reshape(train_mats,[],size(train_mats,3));
        if ~isempty(find(isnan(train_vcts), 1))
            if length(find(isnan(train_vcts))) == size(train_mats,3)*size(train_mats,2)
                train_vcts(isnan(train_vcts)) = 1;
            else
                error_run = [error_run; [yy, ff, ii, i]];
            end
        end
        
        train_behav = all_behav;
        train_behav(leftout) = [];
        
        % correlate all edges with behavior
        
        [r_mat,p_mat] = corr(train_vcts',train_behav, 'rows', 'complete');
        
        r_mat = reshape(r_mat,no_node,no_node);
        p_mat = reshape(p_mat,no_node,no_node);
        
        % set threshold and define masks
        
        pos_mask = zeros(no_node,no_node);
        neg_mask = zeros(no_node,no_node);
        
        pos_edges = find(r_mat > 0 & p_mat < thresh);
        neg_edges = find(r_mat < 0 & p_mat < thresh);
        
        pos_mask(pos_edges) = 1;
        neg_mask(neg_edges) = 1;
        eval(['pos_matri_', idx{ii},'(:,:,leftout) =  pos_mask;'])
        eval(['neg_matri_', idx{ii},'(:,:,leftout) =  neg_mask;'])
        
        
        
        % get sum of all edges in TRAIN subs (divide by 2 to control for the
        % fact that matrices are symmetric)
        
        train_sumpos = zeros(no_sub-1,1);
        train_sumneg = zeros(no_sub-1,1);
        
        for ss = 1:size(train_sumpos)
            train_sumpos(ss) = nansum(nansum(train_mats(:,:,ss).*pos_mask))/2;
            train_sumneg(ss) = nansum(nansum(train_mats(:,:,ss).*neg_mask))/2;
        end
        
        eval([idx{ii},'_pos = train_sumpos;'])
        eval([idx{ii},'_neg = train_sumneg;'])
        
        % run model on TEST sub
        test_mat = all_mats{ii}(:,:,leftout);
        test_sumpos = nansum(nansum(test_mat.*pos_mask))/2;
        test_sumneg = nansum(nansum(test_mat.*neg_mask))/2;
        eval([idx{ii},'_testpos = test_sumpos;'])
        eval([idx{ii},'_testneg = test_sumneg;'])
    end
    
    % build model on TRAIN subs
    [coh_pos_curve, coh_neg_curve, train_behav_curve]=prepareCurveData(coh_pos, coh_neg, train_behav);
    [plv_pos_curve, plv_neg_curve, train_behav_curve]=prepareCurveData(plv_pos, plv_neg, train_behav);
    
    X_pos=[coh_pos_curve,plv_pos_curve,ones(size(train_behav_curve))];
    [fit_pos,bint_pos, r_pos, rint_pos, stats_pos] = regress(train_behav_curve,X_pos);
    
    X_neg=[coh_neg_curve, plv_neg_curve,ones(size(train_behav_curve))];
    [fit_neg,bint_neg, r_neg, rint_neg, stats_neg] = regress(train_behav_curve,X_neg);
    
    X_glm=[coh_pos_curve, coh_neg_curve,plv_pos_curve, plv_neg_curve,ones(size(train_behav_curve))];
    [fit_glm,bint_glm, r_glm, rint_glm, stats_glm] = regress(train_behav_curve,X_glm);
    
    
    behav_pred_pos(leftout) = fit_pos(1)*coh_testpos + fit_pos(2)*plv_testpos + fit_pos(3);
    behav_pred_neg(leftout) = fit_neg(1)*coh_testneg + fit_neg(2)*plv_testneg + fit_neg(3);
    behav_pred_glm(leftout) = fit_glm(1)*coh_testpos + fit_glm(2)*coh_testneg + fit_glm(3)*plv_testpos + fit_glm(4)*plv_testneg +...
        fit_glm(5);
    
    % Calculate the pos_matrics and neg_matrics
    pos_matri(:,:,leftout) = fit_pos(1).*pos_matri_coh(:,:,leftout) + fit_pos(2).*pos_matri_plv(:,:,leftout);
    neg_matri(:,:,leftout) = fit_neg(1).*neg_matri_coh(:,:,leftout) + fit_neg(2).*neg_matri_plv(:,:,leftout);
    
end

pos_matri_sum = pos_matri_coh + pos_matri_plv;
neg_matri_sum = neg_matri_coh + neg_matri_plv;

% compare predicted and observed scores
% positive
[behav_pred_pos_curve, all_behav_curve]=prepareCurveData(behav_pred_pos, all_behav);
if ~isempty(behav_pred_pos_curve) && length(behav_pred_pos_curve) > 3
    [true_prediction_r_pos, P_pos] = corr(behav_pred_pos_curve,all_behav_curve);
else
    [true_prediction_r_pos, P_pos] = corr(behav_pred_pos,all_behav);
end

% negative
[behav_pred_neg_curve, all_behav_curve]=prepareCurveData(behav_pred_neg, all_behav);
if ~isempty(behav_pred_neg_curve) && length(behav_pred_neg_curve) > 3
    [true_prediction_r_neg, P_neg] = corr(behav_pred_neg_curve,all_behav_curve);
else
    [true_prediction_r_neg, P_neg] = corr(behav_pred_neg,all_behav);
end

% glm
[behav_pred_glm_curve, all_behav_curve]=prepareCurveData(behav_pred_glm, all_behav);
if ~isempty(behav_pred_glm_curve) && length(behav_pred_glm_curve) > 3
    [true_prediction_r_glm, P_glm] = corr(behav_pred_glm_curve,all_behav_curve);
else
    [true_prediction_r_glm, P_glm] = corr(behav_pred_glm,all_behav);
end

true_MSE_pos = sum((behav_pred_pos-all_behav).^2)/(no_sub-length(fit_pos)-1);
true_MSE_neg = sum((behav_pred_neg-all_behav).^2)/(no_sub-length(fit_neg)-1);
true_MSE_glm = sum((behav_pred_glm-all_behav).^2)/(no_sub-length(fit_glm)-1);
end





% figure(1); plot(behav_pred_pos,all_behav,'r.'); lsline
% figure(2); plot(behav_pred_neg,all_behav,'b.'); lsline

