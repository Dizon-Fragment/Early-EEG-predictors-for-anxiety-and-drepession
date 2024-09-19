%% Plot the Raincloud pic
clear; close all; clc;
%%% ---------------------- parameters ------------------------
main_path = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes'; 
result_plot = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data/picresult';
symptom_path = fullfile( main_path,  'data','subjset_symptom.mat');
pic_path = fullfile( result_plot,  'Raincloud');
mkdir( pic_path)
code_path = fullfile( result_plot, 'toolbox', 'raincloud');
color_path = fullfile( result_plot, 'toolbox', 'slanCM');
restoredefaultpath;
addpath( code_path)
addpath( genpath( color_path))
ipname = {'preLoo2idx_PERMedge', 'FzCPM2idx_PERMedge'};

%% Specific parameters
symptom = {'SAS','SDS'};
year = {'2017', '2019',  'fz97'}; % change dgz 20240613
freq = {'alpha', 'beta1'};
% subedge = 'neg';
closepic = 0;
rgbcolor = [31,184,207];

%% load the predictive data
datapath = fullfile( main_path, 'data','longitudinalData','csd', 'result', ipname{1});
datapath2 = fullfile( main_path, 'data','longitudinalData','csd', 'result', ipname{2});
load(fullfile(symptom_path)) % get the 'symptom_score' 2023/1/19
T_name = symptom_score.Properties.VariableNames;
T_result = table2cell(symptom_score);
fprintf('\n----------------------------------------------------------\n');
fprintf('%s\n%s\n','Loading symptom and predictive data ......', ' finished.');
data_struct = struct();
for yy = 1:numel(year)
    for s = 1:numel(symptom)
        clear -regexp [^datapath2 data_struct T_result T_name symptom_path rgbcolor freq year symptom main_path result_plot pic_path color_path ipname datapath]
        % load the 'Extractedge_perm'
        if yy <= 2
            load(fullfile(datapath, ['area_', freq{s}, '_', year{yy}, '_0.05_', ipname{1},'.mat']));
        else
            load(fullfile(datapath2, ['area_', freq{s}, '_', year{yy}, '_0.05_', ipname{2},'.mat']));
        end
        tmp_data = eval(['Extractedge_perm.', symptom{s},'(23:25)']); % pos neg glm
        sig_data = eval(['cat(2,Extractedge_perm.', symptom{s},'{1:3})']); % permutation results
        idx = find( sig_data <= 0.05);
        if numel(idx) == 1
            eval(['data_struct.', symptom{s}, '_', freq{s}, '_', year{yy},' = tmp_data{1,idx};']);
        else
            if find(ismember(idx, 3)) % glm
                eval(['data_struct.', symptom{s}, '_', freq{s}, '_', year{yy},' = tmp_data{1,3};']);
            else
                eval(['data_struct.', symptom{s}, '_', freq{s}, '_', year{yy},'_more = tmp_data(1,idx);']);
            end
        end
    end
end
% merge the data
dataToPlot = {};% SAS and SDS
dataPlotName = {};% SAS and SDS
dataFN = fieldnames( data_struct);
addname = {'2015', 'fz75'};
cpmfilename = {'CPM', 'CPM_diff'};
matname = {'ForCPM', 'FisherZ'};
if isempty(find( contains( dataFN, 'more')))
    for s = 1:numel(symptom)
        % add the 2015 fz75 for try
        for aa = 1:numel( addname)
        idx = {'coh','plv'};
            input_dir = fullfile(main_path, 'data','longitudinalData','csd', cpmfilename{aa});
        for ii = 1:length(idx)
            fprintf('----------------------------------------------------------\n');
            fprintf('%s\n',['Frequency: ''', freq{s}, ''' and connectome index: ''', idx{ii}, ''' is processing ......']);
             % get 3 connectome idx matrics
            if aa == 1
                load(fullfile(input_dir, [idx{ii},'_', freq{s}, '_',addname{aa},'_',matname{aa},'.mat'])); % load the ForCPM mat
                eval([idx{ii},'_mtx = sub_',[idx{ii},'_', freq{s}, '_2015'],';'])
            elseif aa == 2
                load(fullfile(input_dir, [idx{ii},'_', freq{s}, '_',matname{aa},'.mat'])); % load the ForCPM mat
                eval([idx{ii},'_mtx = fz75;'])
            end
            clearname = who('-regexp', 'sub_');
            if ~isempty(clearname)
                clear(clearname{:})
            end
            clearname2 = who('-regexp', 'fz');
            if ~isempty(clearname2)
                clear(clearname2{:})
            end
        end
        all_mats  = {coh_mtx, plv_mtx};
        all_behav = cell2mat(T_result(:, ismember(T_name,symptom{s})));
        [~,~,~,behav_pred_pos, behav_pred_neg, behav_pred_glm] = ...
            predict_behavior_longitudinal_forplot(all_mats,all_behav, 0.05,idx);
        if aa == 1
            switch symptom{s}
                case 'SAS'
                    data_2015 = behav_pred_neg;
                case 'SDS'
                    data_2015 = behav_pred_glm;
            end
        elseif aa == 2
            switch symptom{s}
                case 'SAS'
                    data_fz75 = behav_pred_neg;
                case 'SDS'
                    data_fz75 = behav_pred_glm;
            end
        end
        end
        % save the data
        dataToPlot{s} = {eval(['data_struct.', symptom{s}, '_', freq{s}, '_2017'])', ...
            eval(['data_struct.', symptom{s}, '_', freq{s}, '_fz97'])', ...
            eval(['data_struct.', symptom{s}, '_', freq{s}, '_2019'])', ...
            cell2mat(T_result(:, find(ismember(T_name, symptom{s})))')};
        dataPlotName{s} = {['Age 9'], ['Age 9 + age 11'],['Age 11'],  ['Age 13']};

    end
else
    error( 'Please check your <data_struct>, comfirming one column dataset for further studies.')
end

%% Plot the rainclouds
close all;
for s = 1:numel(symptom)
    dataInput = dataToPlot{s};
    dataInputName = dataPlotName{s};
    DB_fillrate = 3.5;
    DotSize = 60;
    markerEdgeAlpha = 0.2;
    markerFaceAlpha = 0.15;
    halfline = 0.15;
    savepath = fullfile(pic_path, [symptom{s}, '3ageBlank.svg']);
    rainCloudsFor_longitudinal(dataInput, dataInputName, DB_fillrate,...
        DotSize, markerEdgeAlpha, markerFaceAlpha,halfline,1, 0,0, 0, 1,savepath)
    % if95ci, ifClose, ifMaxScreen, ifSave, ifBlank, savename
end



