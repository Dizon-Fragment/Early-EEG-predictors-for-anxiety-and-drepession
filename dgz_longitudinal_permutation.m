%% Initialization
clear; close all; clc;

%% Run the longitudial permutation
%%% ---------------------- parameters ------------------------
main_path = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data'; 
result_plot = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data/longitudinalData/csd';
% Fre_file = {'delta','theta','alpha','beta1', 'beta2','gamma'};
Fre_file = {'alpha','beta1'};
idx = {'coh','plv'};
idx_spctrm = {'coh','plv'};
symptom_path = fullfile( main_path,  'subjset_symptom.mat');
folders = {'csd'};
p_thresh = [0.05];
year_data = {'2015', '2017', '2019'};
Fdata = {'area'};
ipname = {'preLoo2idx','FzCPM2idx'};
opname = strcat(ipname, {'_PERMedge'});
figname = strcat(ipname, {'_PERMpic'});

%% create the result dir
for dd = 1:length(folders)
    for i = 1:numel(opname)
        output_dir = fullfile(result_plot, 'result',  opname{i});
        output_pic = fullfile(result_plot, 'result', figname{i});
        if ~exist(fullfile(output_dir), 'dir')
            mkdir(fullfile(output_dir))
        end
        if ~exist(fullfile(output_pic), 'dir')
            mkdir(fullfile(output_pic))
        end
    end
end
warning('off');

%% Rerun the significant permutation results
loop_opname = 1:2; % symptom, Fz
loop_year = 1:3; % 1 2 3
loop_freq = [1:length(Fre_file)];
loop_pthresh = [1:length(p_thresh)];
% for o = 1:numel(opname)
for o = loop_opname
    opname_new = opname{o};
    figname_new = figname{o};
    isFz = contains(lower(opname{o}), 'fz');
    Fdata = {'area'};

    
    for fd = 1:numel(Fdata)
        clear -except [^loop_pthresh loop_freq loop_year loop_opname isFz,figname_new,opname_new, year_data, main_path, symptom_path, folders, Fre_file, idx, idx_spctrm, year_data, Fdata, ipname, opname, figname,p_thresh];
        
        load(fullfile(symptom_path)) % get the 'symptom_score' 2023/1/19
        T_name = symptom_score.Properties.VariableNames;
        T_result = table2cell(symptom_score);
        fprintf('\n----------------------------------------------------------\n');
        fprintf('%s\n','Matchng finished...');
        
        for dd = 1:length(folders) %  csd
            if isFz
                input_dir = fullfile(result_plot, 'CPM_diff');
                year_data = {'fz75', 'fz95', 'fz97'};
            else
                input_dir = fullfile(result_plot,  'CPM'); % path for eeg data
                year_data = {'2015', '2017', '2019'};
            end
            error_run =[];
            %             for yy = 1:length(year_data) % 2015, 2017 or 2019
            for yy = loop_year
%                 for ff = 1:length(Fre_file) % delta, theta and so on
                for ff = loop_freq
                    % loop for two threshold
%                     for tt = 1:length(p_thresh)
                    for tt = loop_pthresh
                        % output result
                        clear -regexp [^loop_pthresh loop_freq loop_year loop_opname isFz,figname_new,opname_new, folders, year_data, main_path, Data_Dir, Dir_path, csd_path, symptom_path, Dpath, Fre_file, idx, idx_spctrm, year_data, Fdata, ipname,opname, T_result, T_name, input_dir, error_run, subdata, p_thresh, figname];
                        txt_dir = fullfile(result_plot, 'result',   ipname{o});
                        filecond = fullfile(txt_dir, [Fdata{fd}, '_',  Fre_file{ff}, '_', year_data{yy}, '_', num2str(p_thresh(tt)), '_', ipname{o}, '.txt']); % input txt
                        output_dir = fullfile(result_plot, 'result',  opname_new);
                        output_fig = fullfile(result_plot, 'result',   figname_new);
                        file_path = fullfile(output_dir, [Fdata{fd}, '_', Fre_file{ff}, '_', year_data{yy}, '_', num2str(p_thresh(tt)), '_', opname_new, '.mat']); % save the result as txt file
                        %                         fig_path = fullfile(output_fig, [Fdata{fd}, '_', Fre_file{ff}, '_', year_data{yy}, '_', num2str(p_thresh(tt)), '_', figname_new, '.txt']); % save the result as txt file
                        
                        if exist( filecond, 'file')
                            abc = importdata(filecond);% change the name of file
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
                            % for symptom
                            if numel(tmp_set(1:end-1)) == 1 % only name without 'std'
                                tmp_abc{tabc,1} = tmp_set{1};
                                tmp_type{tabc,1} = tmp_set{end}(1:end-2);
                            elseif numel(tmp_set(1:end-1)) == 2 % std
                                tmp_abc{tabc,1} = [tmp_set{1}, '_', tmp_set{2}];
                                tmp_type{tabc,1} = tmp_set{end}(1:end-2);
                            else
                                error(['something wrong in ', abc.textdata{tabc,1}]);
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
                        for i = sig_scales % change the scales you need to analyze
                            clear -regexp [^loop_pthresh loop_freq loop_year loop_opname Extractedge_perm savei, sig_p sig_type sig_info, isFz,figname_new,opname_new, folders, year_data, filedraw file_path output_dir sig_scales Main_dir, Data_Dir, Dir_path, csd_path, fmri_path, Dpath, Fre_file, idx, idx_spctrm, year_data, Fdata, ipname,opname, T_result, T_name, input_dir, error_run, subdata, p_thresh,figname];
                            if isFz
                                for ii = 1:length(idx)
                                    fprintf('----------------------------------------------------------\n');
                                    fprintf('%s\n',['Frequency: ''', Fre_file{ff}, ''' and connectome index: ''', idx{ii}, ''' is processing ......'])
                                    load(fullfile(input_dir, [idx{ii},'_', Fre_file{ff}, '_FisherZ.mat']));
                                    fz_list = who('-regexp', '^fz+\d{2}$');
                                    % get 3 connectome idx matrics
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
                                    
                                    load(fullfile(input_dir, [idx{ii},'_', Fre_file{ff}, '_', year_data{yy},'_','ForCPM.mat'])); % load the ForCPM mat
                                    % get 3 connectome idx matrics
                                    eval([idx{ii},'_mtx = sub_',[idx{ii},'_', Fre_file{ff}, '_', year_data{yy}],';'])
                                    clearname = who('-regexp', 'sub_');
                                    if ~isempty(clearname)
                                        clear(clearname{:})
                                    end
                                end
                                all_mats  = {coh_mtx, plv_mtx};
                            end
                            
                            
                            all_behav = cell2mat(T_result(:, i));
                            
                            
                            
                            %                             [true_prediction_r_pos,true_prediction_r_neg,true_prediction_r_glm] = predict_behavior_longitudinal_cnt2idx(all_mats,all_behav, p_thresh(tt),idx);
                            [true_prediction_r_pos,true_prediction_r_neg,true_prediction_r_glm,behav_pred_pos, behav_pred_neg, behav_pred_glm, pos_matri,...
                                neg_matri, true_MSE_pos, true_MSE_neg, true_MSE_glm, error_run, P_pos, P_neg,P_glm,fit_pos,fit_neg,fit_glm,pos_matri_sum,neg_matri_sum] = predict_behavior_longitudinal_forplot(all_mats,all_behav, p_thresh(tt),idx);
                            
                            % Extract edges %GzD, 2021/6/14
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
                            glm_e = (pos_e + neg_e);
                            pos_xy = []; neg_xy = [];
                            [pos_xy(:,1), pos_xy(:,2)] = find(pos_e>=(size(pos_matri_1,3)*0.8)); % 80% subjects
                            [neg_xy(:,1), neg_xy(:,2)] = find(neg_e>=(size(neg_matri_1,3)*0.8));
                            
                            %%% permutation
                            no_iterations = 1000;
                            prediction_r = zeros(no_iterations,3);
                            prediction_r(1,1) = true_prediction_r_pos;
                            prediction_r(1,2) = true_prediction_r_neg;
                            prediction_r(1,3) = true_prediction_r_glm;
                            % MSE
                            prediction_MSE = zeros(no_iterations,3);
                            prediction_MSE(1,1) = true_MSE_pos;
                            prediction_MSE(1,2) = true_MSE_neg;
                            prediction_MSE(1,3) = true_MSE_glm;
                            
                            
                            for it=2:no_iterations
                                %
                                n_idx = find(isnan(all_behav));
                                all_behav_tmp = all_behav;
                                all_behav_tmp(n_idx) = [];
                                all_mats_tmp = all_mats;
                                for aa = 1:numel(all_mats)
                                    all_mats_tmp{aa} = all_mats{aa};
                                    all_mats_tmp{aa}(:,:,n_idx) = [];
                                end
                                no_sub = size(all_mats_tmp{1},3);
                                %
                                fprintf(' Performing iteration %d out of %d\n', it, no_iterations);
                                new_behav                = all_behav_tmp(randperm(no_sub));
                                [prediction_r(it,1),prediction_r(it,2),prediction_r(it,3),~,~,~,~,~,prediction_MSE(it,1),prediction_MSE(it,2),prediction_MSE(it,3)] = ...
                                    predict_behavior_longitudinal_forplot(all_mats_tmp,new_behav, p_thresh(tt),idx);
                            end
                            sorted_prediction_r_pos = sort(prediction_r(:,1),'descend'); % positive
                            sorted_prediction_MSE_pos = sort(prediction_MSE(:,1),'descend');
                            position_pos            = find(sorted_prediction_r_pos==true_prediction_r_pos);
                            position_pos_MSE            = find(sorted_prediction_MSE_pos==true_MSE_pos);
                            if ~isempty(position_pos)
                                pval_pos                = position_pos(1)/no_iterations;
                            else
                                pval_pos = nan;
                            end
                            
                            if ~isempty(position_pos_MSE)
                                pval_pos_MSE                = position_pos_MSE(1)/no_iterations;
                            else
                                pval_pos_MSE = nan;
                            end
                            
                            sorted_prediction_r_neg = sort(prediction_r(:,2),'descend');
                            sorted_prediction_MSE_neg = sort(prediction_MSE(:,2),'descend');
                            position_neg            = find(sorted_prediction_r_neg==true_prediction_r_neg);
                            position_neg_MSE            = find(sorted_prediction_MSE_neg==true_MSE_neg);
                            if ~isempty(position_neg)
                                pval_neg                = position_neg(1)/no_iterations;
                            else
                                pval_neg = nan;
                            end
                            
                            if ~isempty(position_neg_MSE)
                                pval_neg_MSE                = position_neg_MSE(1)/no_iterations;
                            else
                                pval_neg_MSE = nan;
                            end
                            
                            sorted_prediction_r_glm = sort(prediction_r(:,3),'descend');
                            sorted_prediction_MSE_glm = sort(prediction_MSE(:,3),'descend');
                            position_glm            = find(sorted_prediction_r_glm==true_prediction_r_glm);
                            position_glm_MSE            = find(sorted_prediction_MSE_glm==true_MSE_glm);
                            if ~isempty(position_glm)
                                pval_glm                = position_glm(1)/no_iterations;
                            else
                                pval_glm = nan;
                            end
                            
                            if ~isempty(position_glm_MSE)
                                pval_glm_MSE                = position_glm_MSE(1)/no_iterations;
                            else
                                pval_glm_MSE = nan;
                            end
                            
                            tmp_name = T_name{i};
                            eval(['Extractedge_perm.', tmp_name,'{1}=pval_pos;']); % p_value positive
                            eval(['Extractedge_perm.', tmp_name,'{2}=pval_neg;']); % p_value negative
                            eval(['Extractedge_perm.', tmp_name,'{3}=pval_glm;']); % p_value glm
                            eval(['Extractedge_perm.', tmp_name,'{4}=position_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{5}=position_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{6}=position_glm;']);
                            eval(['Extractedge_perm.', tmp_name,'{7}=prediction_r;']); % all distributions (1 2 3)
                            eval(['Extractedge_perm.', tmp_name,'{8}={loop_opname, loop_year, loop_pthresh};']);
                            eval(['Extractedge_perm.', tmp_name,'{9}=sig_info{savei};']);                         
                            eval(['Extractedge_perm.', tmp_name,'{10}=pos_xy;']);
                            eval(['Extractedge_perm.', tmp_name,'{11}=neg_xy;']);
                            eval(['Extractedge_perm.', tmp_name,'{12}=pos_e;']);
                            eval(['Extractedge_perm.', tmp_name,'{13}=neg_e;']);
                            eval(['Extractedge_perm.', tmp_name,'{14}=true_MSE_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{15}=true_MSE_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{16}=true_MSE_glm;']);
                            eval(['Extractedge_perm.', tmp_name,'{17}=P_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{18}=P_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{19}=P_glm;']);
                            eval(['Extractedge_perm.', tmp_name,'{20}=true_prediction_r_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{21}=true_prediction_r_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{22}=true_prediction_r_glm;']);
                            eval(['Extractedge_perm.', tmp_name,'{23}=behav_pred_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{24}=behav_pred_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{25}=behav_pred_glm;']);
                            eval(['Extractedge_perm.', tmp_name,'{26}=fit_pos;']);
                            eval(['Extractedge_perm.', tmp_name,'{27}=fit_neg;']);
                            eval(['Extractedge_perm.', tmp_name,'{28}=fit_glm;']);
                            eval(['Extractedge_perm.', tmp_name,'{29}=pos_matri_1;']); % pos_matri_sum
                            eval(['Extractedge_perm.', tmp_name,'{30}=neg_matri_1;']); % neg_matri_sum
                            eval(['Extractedge_perm.', tmp_name,'{31}=pos_matri;']); % pos_matri_fit
                            eval(['Extractedge_perm.', tmp_name,'{32}=neg_matri;']); % neg_matri_fit
                            eval(['Extractedge_perm.', tmp_name,'{33}=T_result;']);
                            eval(['Extractedge_perm.', tmp_name,'{34}=T_name;']);
                            eval(['Extractedge_perm.', tmp_name,'{35}=pval_pos_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{36}=pval_neg_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{37}=pval_glm_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{38}=position_pos_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{39}=position_neg_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{40}=position_glm_MSE;']);
                            eval(['Extractedge_perm.', tmp_name,'{41}=prediction_MSE;']);
                            savei = savei + 1;
                        end
                        save(fullfile(file_path), 'Extractedge_perm')
                    end
                end
            end
        end
    end
end