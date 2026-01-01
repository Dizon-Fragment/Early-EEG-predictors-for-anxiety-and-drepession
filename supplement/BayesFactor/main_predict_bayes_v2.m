%% Initialization
clear; close all; clc;

%% Rerun the longitudial permutation
%%% ---------------------- parameters ------------------------
main_path = '/data/home/EEG001/Longitudinal/datafromMac/'; 
result_plot = '/data/home/EEG001/Longitudinal/ForPlot/';
Fre_file = {'delta','theta','alpha','beta1', 'beta2','gamma'};
idx = {'coh','plv'};
idx_spctrm = {'coh','plv'};
fmri_path = fullfile( main_path, 'result_Zheyi', 'sig_func_conn_beh2');
symptom_path = fullfile( main_path, 'longitudinal_prefinish', 'subjset_symptom.mat');
folders = {'csd', 'norm'};
p_thresh = [0.05, 0.01];
p_thresh_name = {'p5', 'p1'};
year_data = {'2015', '2017', '2019'};
Fdata = {'area'};
ipname = {'preLooPERM2idx','FzCPMPERM2idx'};
opname = strcat(ipname, {'_BFedge2'});
figname = strcat(ipname, {'_BFpic2'});

%% create the result dir
for dd = 1:length(folders)
    for i = 1:numel(opname)
        output_dir = fullfile(result_plot, 'result', folders{dd}, opname{i});
        output_pic = fullfile(result_plot, 'result', folders{dd}, figname{i});
        bayesian_output_dir = fullfile(output_dir, 'BF');
        if ~exist(output_dir, 'dir'), mkdir(output_dir); end
        if ~exist(output_pic, 'dir'), mkdir(output_pic); end
        if ~exist(bayesian_output_dir, 'dir'), mkdir(bayesian_output_dir); end
    end
end
warning('off');

%% Rerun the significant permutation results
loop_opname = 1:2; % symptom, Fz
loop_year = 1:3; % 1, 2, 3 (2015, 2017, 2019)
loop_freq = [1:length(Fre_file)];
loop_pthresh = [1:length(p_thresh)];

for o = loop_opname
    opname_new = opname{o};
    figname_new = figname{o};
    isfmri = contains(lower(opname{o}), 'fmri');
    isFz = contains(lower(opname{o}), 'fz');
    if isfmri
        Fdata = {'area','net'};
    else
        Fdata = {'area'};
    end
    
    for fd = 1:numel(Fdata)
        clear -except [^isFz, isfmri,figname_new,opname_new, year_data, main_path, symptom_path, folders, Fre_file, idx, idx_spctrm, year_data, Fdata, ipname, opname, figname,p_thresh, loop_year, loop_freq, loop_pthresh, loop_opname, output_dir, bayesian_output_dir];
        if isfmri
            sig_fmri = load(fullfile(fmri_path, ['sig_func_', Fdata{fd}, '_beh_data_n34.mat']));
            tmpa = cellfun(@(x) strfind(sig_fmri.sig_func_corr_name, x, 'ForceCellOutput',false), {'SAS',  'SDS'}, 'UniformOutput', false);
            tmpb = cellfun(@(x) isempty(char(x)) == 0 ,tmpa{1}, 'UniformOutput', true);
            tmpc = cellfun(@(x) isempty(char(x)) == 0 ,tmpa{2}, 'UniformOutput', true);
            idxfmri = tmpb + tmpc;
            T_result = sig_fmri.sig_func_corr_data(:, logical(idxfmri));
            T_name = sig_fmri.sig_func_corr_name(:, logical(idxfmri));
            fprintf('----------------------------------------------------------\n');
            fprintf('%s\n','Matchng finished...');
        else
            load(fullfile(symptom_path)) 
            T_name = symptom_score.Properties.VariableNames;
            T_result = table2cell(symptom_score);
            fprintf('\n----------------------------------------------------------\n');
            fprintf('%s\n','Matchng finished...');
        end
        
        for dd = 1:length(folders) 
            if isFz
                input_dir = fullfile(main_path, 'result',  folders{dd}, 'CPM_diff');
                year_data_loop = {'fz75', 'fz95', 'fz97'};
            else
                input_dir = fullfile(main_path, 'result',  folders{dd}, 'CPM');
                year_data_loop = {'2015', '2017', '2019'};
            end
            error_run =[];
            
            for yy = loop_year
                for ff = loop_freq
                    for tt = loop_pthresh
                        
                        txt_dir = fullfile(main_path, 'result',  folders{dd}, ipname{o});
                        filecond = fullfile(txt_dir, [Fdata{fd}, '_',  Fre_file{ff}, '_', year_data_loop{yy}, '_', num2str(p_thresh(tt)), '_', ipname{o}, '.txt']);
                        output_dir = fullfile(result_plot, 'result',  folders{dd}, opname_new);
                        bayesian_output_dir = fullfile(output_dir, 'BF');
                        file_path = fullfile(output_dir, [Fdata{fd}, '_', Fre_file{ff}, '_', year_data_loop{yy}, '_', num2str(p_thresh(tt)), '_', opname_new, '.mat']);
                        
                        % --- IMPORTANT: Commented out to force overwrite ---
                        % if exist(file_path,'file')
                        %     continue
                        % end
                        
                        if exist( filecond, 'file')
                            abc = importdata(filecond);
                        else
                            continue
                        end
                        
                        if isempty(abc)
                            fclose('all')
                            continue
                        end
                        
                        tmp_abc = {};
                        tmp_type = {};
                        for tabc = 1:length(abc.textdata)
                            tmp_set = strsplit(abc.textdata{tabc,1}, {'_'});
                            if isfmri 
                                if numel(tmp_set(1:end-1)) == 3 
                                    tmp_abc{tabc,1} = [tmp_set{1}, '_', tmp_set{2}, '_', tmp_set{3}];
                                    tmp_type{tabc,1} = tmp_set{end}(1:end-2);
                                elseif numel(tmp_set(1:end-1)) == 4
                                    tmp_abc{tabc,1} = [tmp_set{1}, '_', tmp_set{2}, '_', tmp_set{3}, '_', tmp_set{4}];
                                    tmp_type{tabc,1} = tmp_set{end}(1:end-2);
                                else
                                    error(['something wrong in ', abc.textdata{tabc,1}]);
                                end
                            else 
                                if numel(tmp_set(1:end-1)) == 1 
                                    tmp_abc{tabc,1} = tmp_set{1};
                                    tmp_type{tabc,1} = tmp_set{end}(1:end-2);
                                elseif numel(tmp_set(1:end-1)) == 2 
                                    tmp_abc{tabc,1} = [tmp_set{1}, '_', tmp_set{2}];
                                    tmp_type{tabc,1} = tmp_set{end}(1:end-2);
                                else
                                    error(['something wrong in ', abc.textdata{tabc,1}]);
                                end
                            end
                        end
                        
                        sig_beh = tmp_abc(abc.data<=0.05,1);
                        sig_Beh = unique(sig_beh);
                        sig_scales = zeros(1,length(sig_Beh));
                        sig_p = abc.data(abc.data<=0.05,1);
                        sig_type = tmp_type(abc.data<=0.05,1);
                        
                        if isempty(sig_Beh)
                            fclose('all')
                            continue
                        end
                        
                        tmpscale = cellfun(@(k) strsplit(k), sig_Beh,'UniformOutput', false);
                        a = arrayfun(@(k) tmpscale{k}(1),1:length(tmpscale),'UniformOutput', true);
                        sig_scales = find(ismember(T_name, a) == 1);
                        
                        % find out the sig information
                        tmp_behidx = cellfun(@(x) find(ismember(sig_beh, x)), sig_Beh,  'UniformOutput', false);
                        sig_info = cellfun(@(x) strcat(sig_type(x(1)), {32}, sig_type(x(end)), {32}, num2str(length(x))), tmp_behidx,  'UniformOutput', true);
                        
                        if isempty(sig_scales)
                            fclose('all')
                            continue
                        end
                      
                        savei = 1;
                        Extractedge_perm = struct();
                        for i = sig_scales 
                            clear -regexp [^loop_pthresh loop_freq loop_year loop_opname Extractedge_perm savei, sig_p sig_type sig_info, isFz, isfmri,figname_new,opname_new, folders, year_data, filedraw file_path output_dir bayesian_output_dir sig_scales Main_dir, Data_Dir, Dir_path, csd_path, fmri_path, Dpath, Fre_file, idx, idx_spctrm, year_data, Fdata, ipname,opname, T_result, T_name, input_dir, error_run, subdata, p_thresh,figname, year_data_loop, yy, ff, tt, p_thresh_name];
                            if isFz
                                for ii = 1:length(idx)
                                    fprintf('----------------------------------------------------------\n');
                                    fprintf('%s\n',['Frequency: ''', Fre_file{ff}, ''' and connectome index: ''', idx{ii}, ''' is processing ......'])
                                    load(fullfile(input_dir, [idx{ii},'_', Fre_file{ff}, '_FisherZ.mat']));
                                    fz_list = who('-regexp', '^fz+\d{2}$');
                                    for ll = 1:numel(fz_list)
                                        eval([idx{ii},'_', fz_list{ll},' = ',fz_list{ll},';'])
                                    end
                                    all_fz(ii,:) = strcat(idx(ii), {'_'}, fz_list);
                                    if ~isempty(fz_list)
                                        clear(fz_list{:})
                                        clear fz_list
                                    end
                                end
                                fz_list = all_fz(:,yy);
                                fz_matrics = {eval(fz_list{1}),eval(fz_list{2})};
                                all_mats  = fz_matrics;
                            else
                                clearname = who('-regexp', 'coh_|plv_');
                                if ~isempty(clearname)
                                    clear(clearname{:})
                                end
                                for ii = 1:length(idx)
                                    fprintf('----------------------------------------------------------\n');
                                    fprintf('%s\n',['Frequency: ''', Fre_file{ff}, ''' and connectome index: ''', idx{ii}, ''' is processing ......']);
                                    
                                    load(fullfile(input_dir, [idx{ii},'_', Fre_file{ff}, '_', year_data_loop{yy},'_','ForCPM.mat'])); 
                                    eval([idx{ii},'_mtx = sub_',[idx{ii},'_', Fre_file{ff}, '_', year_data_loop{yy}],';'])
                                    clearname = who('-regexp', 'sub_');
                                    if ~isempty(clearname)
                                        clear(clearname{:})
                                    end
                                end
                                all_mats  = {coh_mtx, plv_mtx};
                            end
                            
                            if isfmri
                                all_behav = T_result(:, i);
                            else
                                all_behav = cell2mat(T_result(:, i));
                            end
                            
                            % 1. Calculate TRUE Prediction
                            [true_prediction_r_pos,true_prediction_r_neg,true_prediction_r_glm,behav_pred_pos, behav_pred_neg, behav_pred_glm, pos_matri,...
                                neg_matri, true_MSE_pos, true_MSE_neg, true_MSE_glm, error_run, P_pos, P_neg,P_glm,fit_pos,fit_neg,fit_glm,pos_matri_sum,neg_matri_sum] = predict_behavior_longitudinal_forplot(all_mats,all_behav, p_thresh(tt),idx);
                            
                            % Extract edges 
                            pos_matri_1 = pos_matri_sum ./2;
                            pos_mask = pos_matri_1(:,:,1);
                            for j =2:size(pos_matri_1,3)
                                pos_mask = pos_matri_1(:,:,j)+ pos_mask;
                            end
                            
                            neg_matri_1 = neg_matri_sum ./2;
                            neg_mask= neg_matri_1(:,:,1);
                            for j =2:size(neg_matri_1,3)
                                neg_mask = neg_matri_1(:,:,j)+ neg_mask;
                            end
                            
                            pos_e = tril(pos_mask,0);
                            neg_e = tril(neg_mask,0);
                            pos_xy = []; neg_xy = [];
                            [pos_xy(:,1), pos_xy(:,2)] = find(pos_e>=(size(pos_matri_1,3)*0.8)); 
                            [neg_xy(:,1), neg_xy(:,2)] = find(neg_e>=(size(neg_matri_1,3)*0.8));
                            
                            %%% 2. PERMUTATION (Corrected Logic)
                            no_iterations = 1000;
                            prediction_r = zeros(no_iterations,3);
                            prediction_MSE = zeros(no_iterations,3);
                            
                            % Place TRUE values in the first row
                            prediction_r(1,:) = [true_prediction_r_pos, true_prediction_r_neg, true_prediction_r_glm];
                            prediction_MSE(1,:) = [true_MSE_pos, true_MSE_neg, true_MSE_glm];
                            
                            parfor it=2:no_iterations
                                % Shuffle Behavior
                                n_idx = find(isnan(all_behav));
                                all_behav_tmp = all_behav;
                                all_behav_tmp(n_idx) = [];
                                all_mats_tmp = all_mats;
                                for aa = 1:numel(all_mats)
                                    all_mats_tmp{aa} = all_mats{aa};
                                    all_mats_tmp{aa}(:,:,n_idx) = [];
                                end
                                no_sub = size(all_mats_tmp{1},3);
                                
                                % fprintf(' Performing iteration %d out of %d\n', it, no_iterations);
                                new_behav = all_behav_tmp(randperm(no_sub));
                                
                                % Run CPM
                                [a,b,c,~,~,~,~,~,e,f,g] = predict_behavior_longitudinal_forplot(all_mats_tmp,new_behav, p_thresh(tt),idx);
                                prediction_r(it,:) = [a,b,c];
                                prediction_MSE(it,:) = [e,f,g];
                            end
                            
                            % --- Calculate P-values using Rank Method (Sorting) ---
                            % R values (Higher is better -> Descend)
                            [~, idx_pos] = sort(prediction_r(:,1), 'descend');
                            rank_pos = find(idx_pos == 1, 1);
                            pval_pos = rank_pos / no_iterations;
                            
                            [~, idx_neg] = sort(prediction_r(:,2), 'descend');
                            rank_neg = find(idx_neg == 1, 1);
                            pval_neg = rank_neg / no_iterations;
                            
                            [~, idx_glm] = sort(prediction_r(:,3), 'descend');
                            rank_glm = find(idx_glm == 1, 1);
                            pval_glm = rank_glm / no_iterations;
                            
                            % MSE values (Lower is better -> Ascend)
                            [~, idx_mse_pos] = sort(prediction_MSE(:,1), 'ascend');
                            rank_mse_pos = find(idx_mse_pos == 1, 1);
                            pval_pos_MSE = rank_mse_pos / no_iterations;
                            
                            [~, idx_mse_neg] = sort(prediction_MSE(:,2), 'ascend');
                            rank_mse_neg = find(idx_mse_neg == 1, 1);
                            pval_neg_MSE = rank_mse_neg / no_iterations;
                            
                            [~, idx_mse_glm] = sort(prediction_MSE(:,3), 'ascend');
                            rank_mse_glm = find(idx_mse_glm == 1, 1);
                            pval_glm_MSE = rank_mse_glm / no_iterations;
                            
                            % Keep position variables for struct compatibility (Rank = Position)
                            position_pos = rank_pos; 
                            position_neg = rank_neg; 
                            position_glm = rank_glm;
                            position_pos_MSE = rank_mse_pos; 
                            position_neg_MSE = rank_mse_neg; 
                            position_glm_MSE = rank_mse_glm;

                            % Save Results
                            tmp_name = T_name{i};
                            eval(['Extractedge_perm.', tmp_name,'{1}=pval_pos;']); 
                            eval(['Extractedge_perm.', tmp_name,'{2}=pval_neg;']); 
                            eval(['Extractedge_perm.', tmp_name,'{3}=pval_glm;']); 
                            eval(['Extractedge_perm.', tmp_name,'{4}=position_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{5}=position_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{6}=position_glm;']);
                            eval(['Extractedge_perm.', tmp_name,'{7}=prediction_r;']); 
                            eval(['Extractedge_perm.', tmp_name,'{8}={loop_opname, loop_year, loop_pthresh};']);
                            eval(['Extractedge_perm.', tmp_name,'{9}=sig_info{savei};']);                         
                            eval(['Extractedge_perm.', tmp_name,'{10}=pos_xy;']);
                            eval(['Extractedge_perm.', tmp_name,'{11}=neg_xy;']);
                            eval(['Extractedge_perm.', tmp_name,'{12}=pos_e;']);
                            eval(['Extractedge_perm.', tmp_name,'{13}=neg_e;']);
                            
                            % Save TRUE MSE
                            eval(['Extractedge_perm.', tmp_name,'{14}=true_MSE_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{15}=true_MSE_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{16}=true_MSE_glm;']);
                            
                            eval(['Extractedge_perm.', tmp_name,'{17}=P_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{18}=P_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{19}=P_glm;']);
                            
                            % Save TRUE R (Indices 20-22)
                            eval(['Extractedge_perm.', tmp_name,'{20}=true_prediction_r_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{21}=true_prediction_r_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{22}=true_prediction_r_glm;']);
                            
                            eval(['Extractedge_perm.', tmp_name,'{23}=behav_pred_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{24}=behav_pred_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{25}=behav_pred_glm;']);
                            eval(['Extractedge_perm.', tmp_name,'{26}=fit_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{27}=fit_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{28}=fit_glm;']);
                            eval(['Extractedge_perm.', tmp_name,'{29}=pos_matri_1;']); 
                            eval(['Extractedge_perm.', tmp_name,'{30}=neg_matri_1;']); 
                            eval(['Extractedge_perm.', tmp_name,'{31}=pos_matri;']); 
                            eval(['Extractedge_perm.', tmp_name,'{32}=neg_matri;']); 
                            eval(['Extractedge_perm.', tmp_name,'{33}=T_result;']);
                            eval(['Extractedge_perm.', tmp_name,'{34}=T_name;']);
                            
                            % Save MSE P-values
                            eval(['Extractedge_perm.', tmp_name,'{35}=pval_pos_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{36}=pval_neg_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{37}=pval_glm_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{38}=position_pos_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{39}=position_neg_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{40}=position_glm_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{41}=prediction_MSE;']);
                            savei = savei + 1;
                            
                            % 3. Bayesian Analysis
                            model_types = {'pos', 'neg', 'glm'};
                            p_perms = [pval_pos, pval_neg, pval_glm];
                            
                            nan_mask_behav = isnan(all_behav);
                            
                            for mt_idx = 1:length(model_types)
                                m_type = model_types{mt_idx};
                                p_perm = p_perms(mt_idx);
                                
                                if ~isnan(p_perm) && p_perm < 0.05
                                    fprintf('  >> CPM Permutation significant for [%s]. Running Bayesian analysis using MCMC...\n', m_type);
                                    
                                    true_behav_final = all_behav(~nan_mask_behav);
                                    eval(['predicted_behav_final = behav_pred_', m_type, ';']);
                                    
                                    % Prepare Data for JASP/Plots
                                    bayesian_table = table(true_behav_final, predicted_behav_final, 'VariableNames', {'True_Scores', 'Predicted_Scores'});
                                    p_thresh_str = p_thresh_name{tt};
                                    jasp_filename = sprintf('ForJASP_CPM_%s_%s_Freq%s_p%s_%s_%s.csv', ...
                                        folders{dd}, year_data_loop{yy}, Fre_file{ff}, p_thresh_str, tmp_name, m_type);
                                    writetable(bayesian_table, fullfile(bayesian_output_dir, jasp_filename));
                                    
                                    % MCMC Sampling
                                    X = [ones(size(predicted_behav_final)), predicted_behav_final];
                                    y = true_behav_final;
                                    n = length(y);
                                    k = size(X,2); 
                                    
                                    logPosterior = @(params) calculate_log_posterior(y, X, params, n, k);
                                    
                                    ols_beta = X\y;
                                    ols_sigma = std(y - X*ols_beta);
                                    start_point = [ols_beta; ols_sigma];
                                    if isinf(logPosterior(start_point)) || isnan(logPosterior(start_point))
                                        fprintf('      Warning: Cannot evaluate log-posterior at start point. Skipping.\n');
                                        continue;
                                    end
                                    
                                    fprintf('      Running MCMC sampler (slicesample)...\n');
                                    num_samples = 5000;
                                    burn_in = 1000;
                                    logPostFunc = @(p) logPosterior(p);
                                    
                                    try
                                        trace = slicesample(start_point, num_samples, 'logpdf', logPostFunc, 'burnin', burn_in, 'thin', 5);
                                    catch
                                        warning('MCMC failed to converge.');
                                        continue;
                                    end
                                    
                                    % Save Bayesian Results
                                    bayes_results = struct();
                                    bayes_results.mcmc_trace = trace; 
                                    bayes_results.slope_posterior_mean = mean(trace(:,2));
                                    bayes_results.intercept_posterior_mean = mean(trace(:,1));
                                    bayes_results.slope_95_credible_interval = quantile(trace(:,2), [0.025, 0.975]);
                                    
                                    eval(['Extractedge_perm.', tmp_name, '{42} = bayes_results;']);
                                    eval(['Extractedge_perm.', tmp_name, '{43} = m_type;']);
                                end
                            end
                            
                        end
                        save(fullfile(file_path), 'Extractedge_perm')
                        fprintf('Saved: %s\n', file_path);
                    end
                end
            end
        end
    end
end

disp('Modeling finished.')


% --- Helper Functions ---

function logP = calculate_log_posterior(y, X, params, n, k)
    beta = params(1:k);
    if isrow(beta), beta = beta'; end
    sigma = params(k+1);
    
    if sigma <= 0
        logP = -inf;
        return;
    end
    
    logLikelihood = -n/2 * log(2*pi) - n * log(sigma) - 1/(2*sigma^2) * sum((y - X * beta).^2);
    logPrior = -log(sigma); 
    
    logP = logLikelihood + logPrior;
end

function [true_prediction_r_pos,true_prediction_r_neg,true_prediction_r_glm,behav_pred_pos, behav_pred_neg, behav_pred_glm, pos_matri,...
    neg_matri, true_MSE_pos, true_MSE_neg, true_MSE_glm, error_run, P_pos, P_neg,P_glm,fit_pos,fit_neg,fit_glm,pos_matri_sum,neg_matri_sum] = predict_behavior_longitudinal_forplot(all_mats,all_behav, thresh, idx)

    warning('off');
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
                    error_run = [error_run; [0, 0, ii, 0]];
                end
            end
            
            train_behav = all_behav;
            train_behav(leftout) = [];
            
            [r_mat,p_mat] = corr(train_vcts',train_behav, 'rows', 'complete');
            r_mat = reshape(r_mat,no_node,no_node);
            p_mat = reshape(p_mat,no_node,no_node);
            
            pos_mask = zeros(no_node,no_node);
            neg_mask = zeros(no_node,no_node);
            
            pos_edges = find(r_mat > 0 & p_mat < thresh);
            neg_edges = find(r_mat < 0 & p_mat < thresh);
            
            pos_mask(pos_edges) = 1;
            neg_mask(neg_edges) = 1;
            eval(['pos_matri_', idx{ii},'(:,:,leftout) =  pos_mask;'])
            eval(['neg_matri_', idx{ii},'(:,:,leftout) =  neg_mask;'])
            
            train_sumpos = zeros(no_sub-1,1);
            train_sumneg = zeros(no_sub-1,1);
            
            for ss = 1:size(train_sumpos)
                train_sumpos(ss) = nansum(nansum(train_mats(:,:,ss).*pos_mask))/2;
                train_sumneg(ss) = nansum(nansum(train_mats(:,:,ss).*neg_mask))/2;
            end
            
            eval([idx{ii},'_pos = train_sumpos;'])
            eval([idx{ii},'_neg = train_sumneg;'])
            
            test_mat = all_mats{ii}(:,:,leftout);
            test_sumpos = nansum(nansum(test_mat.*pos_mask))/2;
            test_sumneg = nansum(nansum(test_mat.*neg_mask))/2;
            eval([idx{ii},'_testpos = test_sumpos;'])
            eval([idx{ii},'_testneg = test_sumneg;'])
        end
        
        [coh_pos_curve, coh_neg_curve, train_behav_curve]=prepareCurveData(coh_pos, coh_neg, train_behav);
        [plv_pos_curve, plv_neg_curve, train_behav_curve]=prepareCurveData(plv_pos, plv_neg, train_behav);
        
        X_pos=[coh_pos_curve,plv_pos_curve,ones(size(train_behav_curve))];
        [fit_pos,~, ~, ~, ~] = regress(train_behav_curve,X_pos);
        
        X_neg=[coh_neg_curve, plv_neg_curve,ones(size(train_behav_curve))];
        [fit_neg,~, ~, ~, ~] = regress(train_behav_curve,X_neg);
        
        X_glm=[coh_pos_curve, coh_neg_curve,plv_pos_curve, plv_neg_curve,ones(size(train_behav_curve))];
        [fit_glm,~, ~, ~, ~] = regress(train_behav_curve,X_glm);
        
        behav_pred_pos(leftout) = fit_pos(1)*coh_testpos + fit_pos(2)*plv_testpos + fit_pos(3);
        behav_pred_neg(leftout) = fit_neg(1)*coh_testneg + fit_neg(2)*plv_testneg + fit_neg(3);
        behav_pred_glm(leftout) = fit_glm(1)*coh_testpos + fit_glm(2)*coh_testneg + fit_glm(3)*plv_testpos + fit_glm(4)*plv_testneg +...
            fit_glm(5);
        
        pos_matri(:,:,leftout) = fit_pos(1).*pos_matri_coh(:,:,leftout) + fit_pos(2).*pos_matri_plv(:,:,leftout);
        neg_matri(:,:,leftout) = fit_neg(1).*neg_matri_coh(:,:,leftout) + fit_neg(2).*neg_matri_plv(:,:,leftout);
    end
    
    pos_matri_sum = pos_matri_coh + pos_matri_plv;
    neg_matri_sum = neg_matri_coh + neg_matri_plv;
                                
    [behav_pred_pos_curve, all_behav_curve]=prepareCurveData(behav_pred_pos, all_behav);
    if ~isempty(behav_pred_pos_curve) && length(behav_pred_pos_curve) > 3
        [true_prediction_r_pos, P_pos] = corr(behav_pred_pos_curve,all_behav_curve);
    else
        [true_prediction_r_pos, P_pos] = corr(behav_pred_pos,all_behav);
    end
    
    [behav_pred_neg_curve, all_behav_curve]=prepareCurveData(behav_pred_neg, all_behav);
    if ~isempty(behav_pred_neg_curve) && length(behav_pred_neg_curve) > 3
        [true_prediction_r_neg, P_neg] = corr(behav_pred_neg_curve,all_behav_curve);
    else
        [true_prediction_r_neg, P_neg] = corr(behav_pred_neg,all_behav);
    end
    
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