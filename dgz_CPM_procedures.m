%% longitudinal CPM prediction
clear;close all;clc;
% -------------------------------------------- Parameter --------------------------------------------------- %
Main_dir = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data'; 
Data_Dir = fullfile(Main_dir, 'longitudinalData');
csd_path = fullfile( Data_Dir, 'csd');
symptom_path = fullfile( Main_dir, 'subjset_symptom.mat');
Dpath = { csd_path};
Fre_file = {'alpha','beta1'};
idx = {'coh','plv'};
idx_spctrm = {'coh','plv'};
year_data = {'2015', '2017', '2019'};
Fdata = {'area','net'};
ipname = ['ForCPM'];
opname = ['preLoo2idx']; 
thresh = [0.05]; 

% --------------------------------------------------------------------------------------------------------------- %

%% make dir for output path
for dd = 1:length(Dpath)
    output_dir = fullfile(Dpath{dd}, 'result',opname);
    if ~exist(fullfile(output_dir), 'dir')
        mkdir(fullfile(output_dir))
    end
end
warning('off');

%%
for fd = 1
    load(fullfile(symptom_path)) % get the 'symptom_score' 2023/1/19
    T_name = symptom_score.Properties.VariableNames;
    T_result = table2cell(symptom_score);
    fprintf('\n----------------------------------------------------------\n');
    fprintf('%s\n','Matchng finished...');
    
    for dd = 1:length(Dpath)
        input_dir = fullfile(Dpath{dd}, 'CPM'); % path for eeg data
        error_run =[];
        for yy = 1:length(year_data)
            for ff = 1:length(Fre_file)
                 for tt = 1:length(thresh)
                        % output result
                        output_dir = fullfile(Dpath{dd}, 'result',opname);
                        file_path = fullfile(output_dir, [Fdata{fd}, '_', Fre_file{ff}, '_', year_data{yy}, '_', num2str(thresh(tt)), '_', opname, '.txt']); % save the result as txt file
                        filedraw = fopen(file_path,'w');
                        
                        for i = 1:size(T_result, 2) % change the scales you need to analyze
                            clearname = who('-regexp', 'coh_|plv_');
                            if ~isempty(clearname)
                                clear(clearname{:})
                            end
                            for ii = 1:length(idx)
                                fprintf('----------------------------------------------------------\n');
                                fprintf('%s\n',['Frequency: ''', Fre_file{ff}, ''' and connectome index: ''', idx{ii}, ''' is processing ......']);
                                
                                load(fullfile(input_dir, [idx{ii},'_', Fre_file{ff}, '_', year_data{yy},'_','ForCPM.mat'])); % load the ForCPM mat
                                % get 2 connectome idx matrics
                                eval([idx{ii},'_mtx = sub_',[idx{ii},'_', Fre_file{ff}, '_', year_data{yy}],';'])
                                clearname = who('-regexp', 'sub_');
                                if ~isempty(clearname)
                                    clear(clearname{:})
                                end
                            end
                            
                            all_behav = cell2mat(T_result(:, i)); 
                            all_mats  = {coh_mtx, plv_mtx};
                            
                            no_sub = size(all_mats{1},3);
                            no_node = size(all_mats{1},1);
                            behav_pred_pos = zeros(no_sub,1);
                            behav_pred_neg = zeros(no_sub,1);
                            behav_pred_glm = zeros(no_sub,1);
                            for leftout = 1:no_sub
                                for ii = 1:numel(all_mats)
                                    fprintf('\n%s\n Leaving out subj # %6.3f for ques %s', ['Frequency: ', Fre_file{ff}, ' and connectome index: ', idx{ii}],leftout, T_name{i});
                                    % ---------------------------------------------------------------
                                    % threshold for feature selection
                                    
                                    % ---------------------------------------------------------------
                                    
                                    % leave out subject from matrices and behavior
                                    
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
                                    
                                    pos_edges = find(r_mat > 0 & p_mat < thresh(tt));
                                    neg_edges = find(r_mat < 0 & p_mat < thresh(tt));
                                    
                                    pos_mask(pos_edges) = 1;
                                    neg_mask(neg_edges) = 1;
                                    
                                    
                                    % get sum of all edges in TRAIN subs (divide by 2 to control for the
                                    % fact that matrices are symmetric)
                                    
                                    train_sumpos = zeros(no_sub-1,1);
                                    train_sumneg = zeros(no_sub-1,1);
                                    
                                    for ss = 1:size(train_sumpos, 1)
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
                            end
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
                            
                            fprintf(filedraw,'%s r:%.3f p:%.3f\n',[T_name{i}, '_positive'], true_prediction_r_pos,P_pos);
                            fprintf(filedraw,'%s r:%.3f p:%.3f\n',[T_name{i}, '_negative'], true_prediction_r_neg,P_neg);
                            fprintf(filedraw,'%s r:%.3f p:%.3f\n',[T_name{i}, '_glm'], true_prediction_r_glm,P_glm);
                            
                            % fclose(filedraw);
                        end
                        fclose('all');
                        clear all_behav all_behav_curve all_mats behav_pred_glm behav_pred_glm_curve behav_pred_neg behav_pred_neg_curve behav_pred_pos behav_pred_pos_curve
                        clear dataP fit_glm fit_neg fit_pos neg_edges neg_mask P_glm p_mat P_neg P_pos pos_edges pos_mask r_mat  test_mat tmp_name
                        clear train_behav train_behav_curve train_mats train_sumneg train_sumneg_curve train_sumpos train_sumpos_curve train_vcts X
                        
                 end
            end
        end
    end
    clear T_result T_name
end
error_run
% check the txt output and error_run


%% Developmental CPM prediction
clear;close all;clc;
% -------------------------------------------- Parameter --------------------------------------------------- %
Main_dir = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data'; 
Data_Dir = fullfile(Main_dir, 'longitudinalData');
csd_path = fullfile( Data_Dir, 'csd');
symptom_path = fullfile( Main_dir, 'subjset_symptom.mat');
Dpath = {csd_path};
Fre_file = {'alpha','beta1'};
idx = {'coh','plv'};
idx_spctrm = {'coh','plv'};
Fdata = {'area','net'};
ipname = ['FisherZ'];
opname = ['FzCPM2idx']; 
thresh = [0.05]; 


% --------------------------------------------------------------------------------------------------------------- %

%% make dir for output path
for dd = 1:length(Dpath)
    output_dir = fullfile(Dpath{dd}, 'result',opname);
    if ~exist(fullfile(output_dir), 'dir')
        mkdir(fullfile(output_dir))
    end
end
warning('off');

%%
% T = readtable(csv_path);
error_run =[];
for fd = 1 % area, only one
    load(fullfile(symptom_path)) % get the 'symptom_score' 2023/1/19
    T_name = symptom_score.Properties.VariableNames;
    T_result = table2cell(symptom_score);
    fprintf('\n----------------------------------------------------------\n');
    fprintf('%s\n','Matchng finished...');
    
    for dd = 1:length(Dpath)
        input_dir = fullfile(Dpath{dd}, 'CPM_diff'); % path for eeg data
        
        for ff = 1:length(Fre_file)
            for tt = 1:length(thresh)
                clear -regexp [^Main_dir, Data_Dir, Dir_path, csd_path, fmri_path, Dpath, Fre_file, idx, idx_spctrm, year_data, Fdata, ipname,opname, T_result, T_name, input_dir, error_run, sub_data, thresh];
                for ii = 1:length(idx)
                    fprintf('----------------------------------------------------------\n');
                    
                    fprintf('%s\n',['Frequency: ''', Fre_file{ff}, ''' and connectome index: ''', idx{ii}, ''' is processing ......']);
                    load(fullfile(input_dir, [idx{ii},'_', Fre_file{ff}, '_', 'FisherZ.mat']));
                    fz_list = who('-regexp', '^fz+\d{2}$');
                    % get 2 connectome idx matrics
                    for ll = 1:numel(fz_list)
                        eval([idx{ii},'_', fz_list{ll},' = ',fz_list{ll},';'])
                    end
                    all_fz(ii,:) = strcat(idx(ii), {'_'}, fz_list);
                    if ~isempty(fz_list)
                        clear(fz_list{:})
                        clear fz_list
                    end
                end
                
                for yy = 1:size(all_fz,2)
                    fz_list = all_fz(:,yy);
                    fz_matrics = {eval(fz_list{1}),eval(fz_list{2})};
                    % output result
                    output_dir = fullfile(Dpath{dd}, 'result',opname);
                    file_path = fullfile(output_dir, [Fdata{fd}, '_',  Fre_file{ff}, '_', fz_list{1}(end-3:end), '_', num2str(thresh(tt)), '_', opname, '.txt']); % save the result as txt file
                    filedraw = fopen(file_path,'w');
                    
                    for i = 1:size(T_result, 2) % change the scales you need to analyze
                        all_behav = cell2mat(T_result(:, i)); % new change 2023/1/19
                        all_mats  = fz_matrics;
             
                        
                        % ---------------------------------------------------------------
                        % threshold for feature selection
                        
                        no_sub = size(all_mats{1},3);
                        no_node = size(all_mats{1},1);
                        
                        behav_pred_pos = zeros(no_sub,1);
                        behav_pred_neg = zeros(no_sub,1);
                        behav_pred_glm = zeros(no_sub,1);
                        % ---------------------------------------------------------------
                        
                        
                        for leftout = 1:no_sub
                            for ii = 1:length(idx)
                                fprintf('%s\n Leaving out subj # %6.3f for ques %s\n', ['Frequency: ', Fre_file{ff}, ' and connectome index: ', idx{ii}],leftout, T_name{i});
                                
                                % leave out subject from matrices and behavior
                                
                                train_mats = all_mats{ii};
                                train_mats(:,:,leftout) = [];
                                train_vcts = reshape(train_mats,[],size(train_mats,3));
                                if ~isempty(find(isnan(train_vcts), 1))
                                    if length(find(isnan(train_vcts))) == size(train_mats,3)*size(train_mats,2)
                                        train_vcts(isnan(train_vcts)) = 1;
                                    else
                                        error_run = [error_run; [dd, yy, ff, ii, i]];
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
                                
                                pos_edges = find(r_mat > 0 & p_mat < thresh(tt));
                                neg_edges = find(r_mat < 0 & p_mat < thresh(tt));
                                
                                pos_mask(pos_edges) = 1;
                                neg_mask(neg_edges) = 1;
                                
                                
                                % get sum of all edges in TRAIN subs (divide by 2 to control for the
                                % fact that matrices are symmetric)
                                
                                train_sumpos = zeros(no_sub-1,1);
                                train_sumneg = zeros(no_sub-1,1);
                                
                                for ss = 1:size(train_sumpos, 1)
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
                        end
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
                        
                        fprintf(filedraw,'%s r:%.3f p:%.3f\n',[T_name{i}, '_positive'], true_prediction_r_pos,P_pos);
                        fprintf(filedraw,'%s r:%.3f p:%.3f\n',[T_name{i}, '_negative'], true_prediction_r_neg,P_neg);
                        fprintf(filedraw,'%s r:%.3f p:%.3f\n',[T_name{i}, '_glm'], true_prediction_r_glm,P_glm);
                        
                        % fclose(filedraw);
                    end
                    fclose('all');
                    clear all_behav all_behav_curve all_mats behav_pred_glm behav_pred_glm_curve behav_pred_neg behav_pred_neg_curve behav_pred_pos behav_pred_pos_curve
                    clear dataP fit_glm fit_neg fit_pos neg_edges neg_mask P_glm p_mat P_neg P_pos pos_edges pos_mask r_mat  test_mat tmp_name
                    clear train_behav train_behav_curve train_mats train_sumneg train_sumneg_curve train_sumpos train_sumpos_curve train_vcts X
                end
            end
        end
    end
    clear T_result T_name
end
error_run
% check the txt output and error_run
