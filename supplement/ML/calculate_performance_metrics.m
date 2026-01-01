function metrics = calculate_performance_metrics(y_true, y_pred)
% Calculates a series of performance metrics for regression models.
%
% Inputs:
%   y_true - (N x 1) vector of true values.
%   y_pred - (N x 1) vector of predicted values.
%
% Outputs:
%   metrics - A struct containing all performance metrics.

    % Remove any potential NaN values to ensure robust calculation.
    valid_idx = ~isnan(y_true) & ~isnan(y_pred);
    y_true = y_true(valid_idx);
    y_pred = y_pred(valid_idx);

    metrics = struct();

    % 1. Pearson Correlation and its p-value.
    [r, p] = corr(y_true, y_pred);
    metrics.r = r;
    metrics.p = p;

    % 2. Mean Squared Error (MSE).
    metrics.mse = mean((y_true - y_pred).^2);

    % 3. Root Mean Squared Error (RMSE).
    metrics.rmse = sqrt(metrics.mse);

    % 4. Mean Absolute Error (MAE).
    metrics.mae = mean(abs(y_true - y_pred));

    % 5. R-squared (Coefficient of Determination).
    ss_total = sum((y_true - mean(y_true)).^2);
    ss_residual = sum((y_true - y_pred).^2);
    metrics.r_squared = 1 - (ss_residual / ss_total);
end