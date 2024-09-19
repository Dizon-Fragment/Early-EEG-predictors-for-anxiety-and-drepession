%% Initialization
clear; close all; clc;

%% Extract the data
%%% ---------------------- parameters ------------------------
main_path = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data'; 
result_plot = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data/longitudinalData/csd';
Fre_file = {'alpha','beta1'};
idx = {'coh','plv'};
idx_spctrm = {'coh','plv'};
symptom_path = fullfile( main_path,  'subjset_symptom.mat');
folders = {'csd'};
p_thresh = [0.05];
year_data = {'2015', '2017', '2019'};
Fdata = {'area'};
opname = {'preLoo2idx','FzCPM2idx'};
ipname = strcat(opname, {'_PERMedge'});




%% Extract the data and plot the correlation result
corr_data = {};
perm_data = {};
for o = 1:length(opname)
    opname_new = ipname{o};
    isFz = contains(lower(ipname{o}), 'fz');
    Fdata = {'area'};
    
    for fd = 1:numel(Fdata)
        clear -except [^isFz,figname_new,opname_new, year_data, main_path, symptom_path, folders, Fre_file, idx, idx_spctrm, year_data, Fdata, ipname, opname, figname,p_thresh];
        
        load(fullfile(symptom_path)) % get the 'symptom_score' 2023/1/19
        T_name = symptom_score.Properties.VariableNames;
        T_result = table2cell(symptom_score);
        fprintf('\n----------------------------------------------------------\n');
        fprintf('%s\n','Matchng finished...');

        for dd = 1:length(folders) %  csd
            if isFz
                input_dir = fullfile(main_path, 'result',  folders{dd}, 'CPM_diff');
                year_data = {'fz75', 'fz95', 'fz97'};
            else
                input_dir = fullfile(main_path, 'result',  folders{dd}, 'CPM'); % path for eeg data
                year_data = {'2015', '2017', '2019'};
            end
            for yy = 1:length(year_data)
                for ff = 1:length(Fre_file)
                    for tt = 1:length(p_thresh)
                        clear -regexp [^isFz,figname_new,opname_new, folders, year_data, main_path, Data_Dir, Dir_path, csd_path, symptom_path, Dpath, Fre_file, idx, idx_spctrm, year_data, Fdata, ipname,opname, T_result, T_name, input_dir, error_run, subdata, p_thresh, figname];
                        output_dir = fullfile(result_plot, 'result', opname_new);
                        filecond = fullfile(output_dir, [Fdata{fd}, '_', Fre_file{ff}, '_', year_data{yy}, '_', num2str(p_thresh(tt)), '_', opname_new, '.mat']); % input txt
                        if exist(filecond, 'file')
                            load(filecond) % get the Extractedge_perm
                        else
                            continue
                        end
                        % extract the data
                        symptomname = fieldnames(Extractedge_perm);
                        timename = [Fdata{fd}, '_', Fre_file{ff}, '_', year_data{yy}, '_', num2str(p_thresh(tt))];
                        for ss = 1:length(symptomname)
                            fprintf('-----------------------------------------------------------\n')
                            fprintf('Loading and extracting ......\n')
                            tmp_symptom = cell2mat(T_result(:, ~cellfun(@isempty, (regexp(T_name, [symptomname{ss},'$'])))));
                            tmp_perm = eval( ['Extractedge_perm.', symptomname{ss}]);
                            corr_data = [corr_data; [{opname_new}, {timename}, symptomname(ss), {cat(2,tmp_perm{17:19})}, ...
                                {cat(2,tmp_perm{20:22})},[tmp_symptom, cat(2,tmp_perm{23:25})], {cat(2,tmp_perm(10:13))}]];
                            corr_info = {'Pvalue_r', 'Rvalue', 'Behav_pred', 'SigEdge'};
                            perm_data = [perm_data; [{opname_new}, {timename}, symptomname(ss), {cat(2,tmp_perm{1:3})}, ...
                                {cat(2,tmp_perm{4:6})},tmp_perm(7), {cat(2,tmp_perm{14:16})}]];
                            perm_info = {'Pvalue_perm', 'Pposition', 'Rdistribution','MSE'};
                        
                        end
                    end
                end
            end
        end
    end
    outpdir = fullfile( result_plot, 'result', 'dataExtraction'); mkdir(outpdir);
    save(fullfile( outpdir, [opname{o},'_Extraction.mat']), 'corr_data', 'corr_info', 'perm_data', 'perm_info')
end

