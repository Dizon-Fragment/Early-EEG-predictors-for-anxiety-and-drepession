function logP = calculate_log_posterior(y, X, params, n, k)
    % This helper function calculates the log posterior, ensuring correct dimensions.
    
    % --- FIX: Reshape beta part of params to be a column vector ---
    beta = params(1:k);
    if isrow(beta)
        beta = beta'; % Transpose if it's a row vector
    end
    
    sigma = params(k+1);
    
    % Check for invalid parameter values which can occur during sampling
    if sigma <= 0
        logP = -inf;
        return;
    end
    
    logLikelihood = -n/2 * log(2*pi) - n * log(sigma) - 1/(2*sigma^2) * sum((y - X * beta).^2);
    logPrior = -log(sigma); % Part of the non-informative prior 1/sigma^2
    
    logP = logLikelihood + logPrior;
end