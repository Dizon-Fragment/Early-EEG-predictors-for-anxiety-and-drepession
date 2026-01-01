function [behav_pred, all_models] = predict_with_loocv(all_mats, all_behav, model_name, p_thresh, hyperparam)
% Added 'hyperparam' argument to pass specific parameters (e.g., Lambda, C)
    
    if nargin < 5, hyperparam = []; end % Handle case with no params

    n_subjects = size(all_mats{1}, 3);
    n_nodes = size(all_mats{1}, 1);
    n_modalities = numel(all_mats);
    
    behav_pred = zeros(n_subjects, 1);
    all_models = cell(n_subjects, 1);
    
    % Pre-allocate feature matrices
    % (Logic kept same as before, simplified for brevity here)
    
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

        % Feature Selection (CPM)
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
        
        % Summation
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
        
        % Feature Subset Selection (_pos, _neg, _all)
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

        % Train and Predict (Passing hyperparam)
        [y_pred, model] = train_and_predict(base_model, X_train_sel, train_behav, X_test_sel, hyperparam);
        
        behav_pred(i) = y_pred;
        all_models{i} = model;
    end
end