function [results] = find_best_predictive_model(s1_model1, s1_model2, Y_target, no_iterations)
    % v5.5 Modified + MSE: Calculates MSE for the best model.
    
    results = struct();
    nets = {'pos', 'neg', 'glm'};
    n_sub = length(Y_target);
    
    for i = 1:length(nets)
        for j = 1:length(nets)
            combo_name = [nets{i}, '_vs_', nets{j}];
            
            X1 = s1_model1.(nets{i}).behav_pred;
            X2 = s1_model2.(nets{j}).behav_pred;
            
            model_formulas = {'Y ~ X1', 'Y ~ X2', 'Y ~ X1 + X2', 'Y ~ X1*X2'};
            r_cv_values = [-Inf, -Inf, -Inf, -Inf];
            mse_cv_values = [NaN, NaN, NaN, NaN]; % Store MSEs
            
            % --- 1. Evaluate all 4 models ---
            for m = 1:4
                y_pred_cv = zeros(n_sub, 1);
                
                % Prepare Design Matrix
                if m == 1, X_design = [ones(n_sub,1), X1];
                elseif m == 2, X_design = [ones(n_sub,1), X2];
                elseif m == 3, X_design = [ones(n_sub,1), X1, X2];
                elseif m == 4, X_design = [ones(n_sub,1), X1, X2, X1.*X2];
                end
                
                for k = 1:n_sub
                    train_idx = true(n_sub, 1); train_idx(k) = false;
                    b = X_design(train_idx, :) \ Y_target(train_idx);
                    y_pred_cv(k) = X_design(k, :) * b;
                end
                
                r_cv_values(m) = corr(y_pred_cv, Y_target, 'rows', 'complete');
                mse_cv_values(m) = mean((y_pred_cv - Y_target).^2); % Calculate MSE
            end
            
            % --- 2. Select Winner ---
            [best_r, best_idx] = max(r_cv_values);
            best_formula = model_formulas{best_idx};
            
            % --- 3. Save Results ---
            res.best_formula = best_formula;
            res.best_r = best_r;
            
            % === ADDED: Save Best MSE ===
            res.mse = mse_cv_values(best_idx);
            % ============================

            % Fit Full Model for Betas
            tbl_full = table(X1, X2, Y_target, 'VariableNames', {'X1', 'X2', 'Y'});
            best_mdl_full = fitlm(tbl_full, best_formula);
            [~, res.p_ana] = corr(predict(best_mdl_full, tbl_full), Y_target);
            
            % Extract Coefficients (Inline Logic)
            coeffs = best_mdl_full.Coefficients;
            res.beta_Intercept = coeffs.Estimate(1);
            res.beta_Net1 = NaN; res.beta_Net2 = NaN; res.beta_Interaction = NaN;
            row_names = coeffs.Properties.RowNames;
            
            idx_X1 = find(strcmp(row_names, 'X1')); if ~isempty(idx_X1), res.beta_Net1 = coeffs.Estimate(idx_X1); end
            idx_X2 = find(strcmp(row_names, 'X2')); if ~isempty(idx_X2), res.beta_Net2 = coeffs.Estimate(idx_X2); end
            idx_Int = find(strcmp(row_names, 'X1:X2')); 
            if ~isempty(idx_Int), res.beta_Interaction = coeffs.Estimate(idx_Int); res.p_Interaction = coeffs.pValue(idx_Int);
            else, res.p_Interaction = NaN; end
            
            % --- 4. Permutation ---
            r_null = zeros(no_iterations, 1);
            parfor it = 1:no_iterations
                Y_shuf = Y_target(randperm(n_sub));
                y_pred_null = zeros(n_sub, 1);
                
                if best_idx == 1, X_d = [ones(n_sub,1), X1];
                elseif best_idx == 2, X_d = [ones(n_sub,1), X2];
                elseif best_idx == 3, X_d = [ones(n_sub,1), X1, X2];
                elseif best_idx == 4, X_d = [ones(n_sub,1), X1, X2, X1.*X2];
                end
                
                for k = 1:n_sub
                    train_idx = true(n_sub, 1); train_idx(k) = false;
                    b_null = X_d(train_idx, :) \ Y_shuf(train_idx);
                    y_pred_null(k) = X_d(k, :) * b_null;
                end
                r_null(it) = corr(y_pred_null, Y_shuf);
            end
            res.perm_p = (sum(r_null >= best_r) + 1) / (no_iterations + 1);
            
            results.(combo_name) = res;
        end
    end
end