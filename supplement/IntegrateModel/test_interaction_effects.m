function [results] = test_interaction_effects(s1_model1, s1_model2, Y_target, no_iterations)
    % v5.5 Modified + MSE: Adds LOOCV for interaction weights AND calculates MSE.
    
    results = struct();
    nets = {'pos', 'neg', 'glm'};
    n_sub = length(Y_target);
    
    for i = 1:length(nets)
        for j = 1:length(nets)
            combo_name = [nets{i}, '_vs_', nets{j}];
            
            % Get Brain Scores (Already from Stage 1 LOOCV)
            X1 = s1_model1.(nets{i}).behav_pred; 
            X2 = s1_model2.(nets{j}).behav_pred; 
            
            % --- 1. Calculate TRUE Predictive R & MSE using LOOCV ---
            y_pred_cv = zeros(n_sub, 1);
            
            for k = 1:n_sub
                train_idx = true(n_sub, 1); train_idx(k) = false;
                
                % Train data
                tbl_train = table(X1(train_idx), X2(train_idx), Y_target(train_idx), 'VariableNames', {'X1', 'X2', 'Y'});
                mdl_k = fitlm(tbl_train, 'Y ~ X1*X2');
                
                % Predict
                y_pred_cv(k) = predict(mdl_k, table(X1(k), X2(k), 'VariableNames', {'X1', 'X2'}));
            end
            
            % Calculate Predictive Metrics
            [res.r, res.p_ana] = corr(y_pred_cv, Y_target, 'rows', 'complete');
            
            % === ADDED: Calculate MSE ===
            res.mse = mean((y_pred_cv - Y_target).^2); 
            % ============================
            
            % --- 2. Extract Representative Betas (from Full Model) ---
            tbl_full = table(X1, X2, Y_target, 'VariableNames', {'X1', 'X2', 'Y'});
            mdl_full = fitlm(tbl_full, 'Y ~ X1*X2');
            
            res.model_formula = 'Y ~ X1 * X2';
            res.r_squared = mdl_full.Rsquared.Ordinary; 
            
            coeffs = mdl_full.Coefficients;
            res.beta_Intercept = coeffs.Estimate(1);
            res.beta_Net1 = coeffs.Estimate(2);      
            res.beta_Net2 = coeffs.Estimate(3);      
            if height(coeffs) >= 4
                res.beta_Interaction = coeffs.Estimate(4);
                res.p_Interaction = coeffs.pValue(4);
            else
                res.beta_Interaction = NaN; res.p_Interaction = NaN;
            end

            % --- 3. Permutation Test ---
            r_null = zeros(no_iterations, 1);
            X_mat = [ones(n_sub,1), X1, X2, X1.*X2]; % Pre-build for speed
            
            parfor it = 1:no_iterations
                Y_shuf = Y_target(randperm(n_sub));
                y_pred_null = zeros(n_sub, 1);
                for k = 1:n_sub
                    train_idx = true(n_sub, 1); train_idx(k) = false;
                    b_null = X_mat(train_idx, :) \ Y_shuf(train_idx);
                    y_pred_null(k) = X_mat(k, :) * b_null;
                end
                r_null(it) = corr(y_pred_null, Y_shuf);
            end
            
            res.perm_p = (sum(r_null >= res.r) + 1) / (no_iterations + 1);
            results.(combo_name) = res;
        end
    end
end
