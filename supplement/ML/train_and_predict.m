% =========================================================================
% FILE: predict_with_loocv.m (Updated for Tuning)
% =========================================================================


% =========================================================================
% FILE: train_and_predict.m (Updated for Tuning)
% =========================================================================
function [y_pred, model] = train_and_predict(model_name, X_train, y_train, X_test, hyperparam)
% Accepts 'hyperparam' to override defaults.
% OLS still skips standardization. Others standardize.

    % 1. NaNs
    if any(isnan(X_train(:))), X_train(isnan(X_train)) = 0; end
    if any(isnan(X_test(:))), X_test(isnan(X_test)) = 0; end
    
    % 2. Standardization (Skip for OLS)
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
            % Default Lambda is 0.1 if not provided
            lam = 0.1; 
            if ~isempty(hyperparam), lam = hyperparam; end
            
            model = fitrlinear(X_train, y_train, 'Learner', 'leastsquares', ...
                'Regularization', 'ridge', 'Lambda', lam);
            y_pred = predict(model, X_test);

        case 'lasso'
            % For Lasso, we force it to use the provided Lambda if given
            if ~isempty(hyperparam)
                % Train with specific lambda
                [B, FitInfo] = lasso(X_train, y_train, 'Alpha', 1, 'Lambda', hyperparam);
                b = B(:,1); intercept = FitInfo.Intercept(1);
            else
                % Default: Auto-select min MSE
                [B, FitInfo] = lasso(X_train, y_train, 'Alpha', 1, 'NumLambda', 20);
                [~, idx] = min(FitInfo.MSE);
                if all(B(:, idx) == 0), non_zeros = find(any(B~=0,1)); if ~isempty(non_zeros), idx=non_zeros(end); else, idx=length(FitInfo.Lambda); end; end
                b = B(:, idx); intercept = FitInfo.Intercept(idx);
            end
            y_pred = X_test * b + intercept;

        case 'svm'
            % Default BoxConstraint (C) is 1
            C = 1;
            if ~isempty(hyperparam), C = hyperparam; end
            
            model = fitrsvm(X_train, y_train, 'KernelFunction', 'linear', ...
                'Standardize', false, 'BoxConstraint', C);
            y_pred = predict(model, X_test);
            
        case 'randomforest'
            % Default MinLeafSize is 5
            leaf = 5;
            if ~isempty(hyperparam), leaf = hyperparam; end
            
            model = TreeBagger(50, X_train, y_train, 'Method', 'regression', ...
                'OOBPrediction', 'on', 'MinLeafSize', leaf, 'NumPredictorsToSample', 'all');
            y_pred = predict(model, X_test);

        otherwise
            error('Unknown model: %s', model_name);
    end
end