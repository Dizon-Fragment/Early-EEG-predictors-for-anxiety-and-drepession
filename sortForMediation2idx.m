%% Data sort
clear; close all; clc;
% -----------------------------------------------------------------------
% Copyright 2024 Guangzhi Deng.
% -----------------------------------------------------------------------
% spm12_path = '/home/GuangzhiDeng/spm12/';
m3_path = 'D:\Installation_package\MediationToolbox-master';
m3_supt_path = 'D:\Installation_package\CanlabCore-master';

restoredefaultpath;
addpath(genpath(m3_path));
addpath(genpath(m3_supt_path));
% addpath(genpath(spm12_path));

% a  b  c'  c  ab	
%% find out some specific significant result
clear; close all; 
%%% ----------------- Parameter for finding ---------------------
% name_pat = {'csd', 'beta1', '97', 'p5', 'SDS', 'a212_37'};
name_pat = {'csd', 'alpha', '19', 'p5', 'SAS'};
% name_pat = {'csd','plv|wpli', 'alpha', '17', 'SAS','a212_37'};
% name_pat = {'coh|plv|wpli', 'beta1', '97', 'SDS','Putamen_Limbic'};
isFz = 0; % 1 for Fz, 0 for symptom
ifplot = 0; % 1 for plot, 0 for no plot. If you choose to plot, please enter the name_pat with a correct order
plotxlsx = 'pos'; % enter the filename when you choose to plot
ifclearhistoryfig = 0;% 1 for clear, 0 for keep
%%% -------------------------------------------------------------

main_path = 'D:\Longitudinal_EEG\datafromMac\result\ForMediation\For2idx';
xlsx_path = 'D:\Longitudinal_EEG\datafromMac\result\ForMediation\For2idx';
fig_path = fullfile( main_path, 'fig');
fmri_path = fullfile( main_path, 'result_Zheyi', 'sig_func_conn_beh2'); 
symptom_path = fullfile( main_path, 'longitudinal_prefinish', 'subjset_symptom.mat');
montage1 = {'csd', 'norm'};
Fre_file = {'delta','theta','alpha','beta1', 'beta2','gamma', 'all_band', 'all_band_50'};
% cnidx = {'coh','plv', 'wpli'};
Fz_year1 = {'2015', '2017', '2019'};
Fz_year2 = {'fz75', 'fz97', 'fz95'};
thresh = [0.05, 0.01]; % new change 2023/1/19
Tname = {'SAS', 'SAS_std', 'SDS', 'SDS_std'};
new_Fre = {'d', 't', 'a', 'b1', 'b2', 'g', 'a3', 'a5'};
new_Tname = {'A', 'As', 'D', 'Ds'};
new_thre = [5,1];
str_thre = {'p5','p1'};
addpath(genpath(main_path))
% what(fullfile(main_path))
cd(main_path)
if ifclearhistoryfig
    rmdir(fig_path, 's')
end
mkdir(fig_path)

result_dir = dir(fullfile(main_path));
result_folder = {result_dir(cellfun(@(x) x==1,{result_dir.isdir})).name};
result_folder(1:2) = []; % delete the '.' and '..'

%
if isFz
    result_folder = result_folder(contains(result_folder, 'fzresult'));
else
    result_folder = result_folder(~contains(result_folder, 'fzresult'));
end
% Searching
result_table = table();
sig_cell = {};
name_pat1 = strcat(name_pat(1:end-1),{'|'});
name_pat2 = strcat(name_pat1{:}, name_pat{end});


for rr = 1:numel(result_folder)
    tmp_folder = fullfile(main_path, result_folder{rr});
    folder_dir = dir(fullfile(tmp_folder, '*.mat'));
    if isempty(folder_dir)
        continue
    end
    matname = {folder_dir.name};
    searchidx = cellfun(@(x) numel(regexp(x, name_pat2)), matname, 'UniformOutput', true);
    sigrltidx = find(searchidx >= numel(name_pat));
    if isempty(sigrltidx)
        continue
    end
    % save the result to the cell and the table
    cellrowidx = size(sig_cell, 1);
    sig_cell(cellrowidx+1:cellrowidx+numel(sigrltidx),1) = result_folder(rr);
    % load the result
    for ss = 1:numel(sigrltidx)
        clear -regexp [^fig_path montage1,Fre_file,cnidx,thresh,Tname,new_Fre,new_Tname,new_thre,str_thre,Fz_year1,Fz_year2,isPlot,isFz cellrowidx sigrltidx searchidx matname folder_dir tmp_folder name_pat2 sig_cell name_pat1 name_pat result_folder main_path]
        sig_cell(cellrowidx+ss,2) = {matname{sigrltidx(ss)}(1:end-4)};
        load(fullfile(main_path, result_folder{rr}, matname{sigrltidx(ss)}))
        sig_cell(cellrowidx+ss,3) = {paths(end)};
        sig_cell(cellrowidx+ss,4) = {stats.p(end)};
        sig_cell(cellrowidx+ss,5) = {squeeze(stats.ci(:,end,:))'};
    end
end
if ~isempty(sig_cell)
    fprintf('\nSearching for the pattern: %s\n', name_pat2)
    fprintf('--------------------------------------------------------------\n')
    result_table = cell2table(sig_cell, 'VariableNames',{'FilesFolders', 'Filenames', 'Coeffciences', 'Pvalues', 'CIs'});
    result_table
else
    fprintf('\nSearching for the pattern: %s\n', name_pat2)
    fprintf('--------------------------------------------------------------\n')
    fprintf('No result in folders. TAT\n')
    for rr = 1:numel(result_folder)
        tmp_folder = fullfile(main_path, result_folder{rr});
        folder_dir = dir(fullfile(tmp_folder, '*.mat'));
        if isempty(folder_dir)
            continue
        end
        matname = {folder_dir.name};
        swname = matname(contains(matname, 'somethingwrong'));
        if isempty(swname)
           continue 
        end
        % save the result to the cell and the table
        cellrowidx = size(sig_cell, 1);
        sig_cell(cellrowidx+1:cellrowidx+numel(sigrltidx),1) = result_folder(rr);
        % load the result
        for ss = 1:numel(swname)
            clear -regexp [^fig_path montage1,Fre_file,cnidx,thresh,Tname,new_Fre,new_Tname,new_thre,str_thre,Fz_year1,Fz_year2,isPlot,isFz swname cellrowidx sigrltidx searchidx matname folder_dir tmp_folder name_pat2 sig_cell name_pat1 name_pat result_folder main_path]
            load(fullfile(main_path, result_folder{rr}, swname{ss}))
            % errorname, mismatching, missname
            if isempty(mismatching)
                continue
            end
            ciname = mismatching(cellfun(@(x) ~isempty(regexp(x, 'ci$', 'match')), mismatching, 'UniformOutput', true));
%             cnidx = {'coh','plv', 'wpli'};
            if isFz
                Fz_year = {'fz75', 'fz95', 'fz97'};
            else
                Fz_year = {'2015', '2017', '2019'};
            end
            Fre_file = {'delta','theta','alpha','beta1', 'beta2','gamma', 'all_band', 'all_band_50'};
            Tname = {'SAS', 'SAS_std', 'SDS', 'SDS_std'};
            new_Fre = {'d', 't', 'a', 'b1', 'b2', 'g', 'a3', 'a5'};
            new_Tname = {'A', 'As', 'D', 'Ds'};
            new_thre = [5,1];
            newname_spt = strsplit(name_pat2, '|');
            numidx = numel(find(cellfun(@(x) numel(regexp(x, 'coh|plv|wpli')), newname_spt, 'UniformOutput', true)));
            if numidx > 1
                cntidx(1,:) = newname_spt(1:numidx);
            else
                cntidx(1,:) = newname_spt(1);
            end
            % change the frequency band abbreviate
            newname_spt2 = newname_spt;
            newname_freband = cellfun(@(x) ~isempty(regexp(x, 'delta|theta|alpha|beta1|beta2|gamma|all_band|all_band_50', 'once')), newname_spt);
            newname_spt2(newname_freband) = new_Fre(contains(Fre_file, newname_spt2{newname_freband}));
            % change the Tname (symptom)
            newname_symptom = cellfun(@(x) ~isempty(regexp(x, 'SAS$|SDS$|SAS_std$|SDS_std$', 'once')), newname_spt);
            newname_spt2(newname_symptom) = new_Tname(cellfun(@(x) ~isempty(regexp(x, [newname_spt2{newname_symptom},'$'], 'once')), Tname));
            % chage the connectome idx
            newname_spt2(1:numidx) = cellfun(@(x) x(1), newname_spt2(1:numidx), 'UniformOutput', false);
            % change the threshold
            newname_spt2 = cellfun(@(x) regexprep(x, 'p+\d', x(end)), newname_spt2, 'UniformOutput', false);
            % combine for the regexp pattern
            newname_pat1 = strcat(newname_spt2(numidx:end-1),{'_|_'});
            newname_patidx_pre = strcat(newname_spt2{1:numidx-1},{'_|'});
            newname_pat2 = strcat(newname_patidx_pre, newname_pat1{:}, newname_spt2{end});
            % search
            searchidx = cellfun(@(x) numel(regexp(x(3:end-18), newname_pat2)), ciname, 'UniformOutput', true);
            sigrltidx = find(searchidx >= numel(newname_spt2)-numidx); % e.g. 7 - 1 (-2 + 1), exclude a212_37
            if isempty(sigrltidx)
                continue
            end
            % save the result to the cell and the table
            cellrowidx = size(sig_cell, 1);
            sig_cell(cellrowidx+1:cellrowidx+numel(sigrltidx),1) = result_folder(rr);
            % load the result
            for ss = 1:numel(sigrltidx)
                clear -regexp [^fig_path montage1,Fre_file,cnidx,thresh,Tname,new_Fre,new_Tname,new_thre,str_thre,Fz_year1,Fz_year2,isPlot,isFz ciname cellrowidx sigrltidx searchidx matname folder_dir tmp_folder name_pat2 sig_cell name_pat1 name_pat result_folder main_path]
                sig_cell(cellrowidx+ss,2) = {ciname{sigrltidx(ss)}};

            end
            
            
        end
        
    end
    if ~isempty(sig_cell)
        fprintf('\nSearching for the new pattern: %s\n', newname_pat2)
        fprintf('--------------------------------------------------------------\n')
        result_table = cell2table(sig_cell, 'VariableNames',{'FilesFolders', 'Filenames'});
        result_table
    else
        
        fprintf('No result in mismatching files as well. TAT\n')
    end

end

% for plot
if ~isempty(sig_cell)&& ifplot
    fprintf('\nLoading the data for plotting ......\n')
    data_dir = dir(fullfile(xlsx_path, '*.xlsx'));
    xlsxname = {data_dir.name};
    if isFz
        Fzname = xlsxname(contains(xlsxname, 'Fz'));
        Fz_year = Fz_year2;
    else
%         Fzname = xlsxname(~contains(xlsxname, {'Fz','1', '0', '2'})); % exclude Fz and 1 xlsx
        Fzname = xlsxname(~contains(xlsxname, {'Fz','1', '0'})); % exclude Fz and 1 xlsx
        Fz_year = Fz_year1;
    end
    % load the specific xlsx
    sp_xlsx = Fzname(~cellfun('isempty', regexp( Fzname,strcat({'^'},plotxlsx))));
    if numel(sp_xlsx)~=1
        error('Please check your xlsx name.')
    end
    T = readtable(fullfile(main_path, sp_xlsx{1}));
    % renew the vname 
    vname = T.Properties.VariableNames;
    spt_vname = cellfun(@(x) strsplit(x, '_'), vname, 'UniformOutput', false);
%     cnt_name = spt_vname(cellfun(@(x) numel(x)==6, spt_vname)); % extract the connectome names
    cnt_name = spt_vname(cellfun(@(x) numel(x)==5, spt_vname)); % extract the connectome names
    % check the dimension of each parameter
    m_ab = cellfun(@(x) x(1), montage1,'UniformOutput', false);
    y_ab = cellfun(@(x) x(end-1:end), Fz_year,'UniformOutput', false);
%     c_ab = cellfun(@(x) x(1), cnidx,'UniformOutput', false);
    t_ab = arrayfun(@(x) num2str(x), new_thre,'UniformOutput', false);
    all_ab = {m_ab, new_Fre, y_ab, t_ab, new_Tname};
%     all_uncover = {montage1, cnidx, Fre_file, Fz_year, str_thre, Tname}; 
    all_uncover = {montage1, Fre_file, Fz_year, str_thre, Tname}; 
    all_para_dim = {};
    for i = 1:numel(cnt_name{1})
        eval(['para',num2str(i),'_dim = cellfun(@(y) isequal(sort(y), sort(unique(cellfun(@(x) x{',num2str(i), '},cnt_name, ''UniformOutput'', false)))), all_ab);'])
        eval(['all_para_dim = [all_para_dim, para',num2str(i),'_dim];'])
        tmp_para_dim = eval(['para',num2str(i),'_dim']);
        if eval(['numel(find(para',num2str(i),'_dim)) == 1'])
            % find out the uncover name
            eval(['para',num2str(i),'_uncover = cellfun(@(x) all_uncover{tmp_para_dim}(~cellfun(''isempty'',',...
                'regexp(all_ab{tmp_para_dim},strcat(x(tmp_para_dim), ''$'')))),cnt_name,''UniformOutput'', true);']) 
        else
            error('Check your code and parameters')
        end
    end
    new_vname = [strcat(para1_uncover, '_',para2_uncover, '_',para3_uncover, '_',para4_uncover, '_',para5_uncover),...
        vname(cellfun(@(x) numel(x)~=5, spt_vname))];
%     new_vname = [strcat(para1_uncover, '_',para2_uncover, '_',para3_uncover, '_',para4_uncover, '_',para5_uncover, '_',para6_uncover),...
%         vname(cellfun(@(x) numel(x)~=6, spt_vname))];
    % search for the specific pattern
    name_pat2_spt = strsplit(name_pat2, '|');
    name_pat2_idx = contains(name_pat2_spt, Tname);
    name_pat2_spt(name_pat2_idx) = strcat(name_pat2_spt(name_pat2_idx),'$');
    name_pat2_add = strcat(name_pat2_spt(1:end-1),{'|'});
    name_pat3 = strcat(name_pat2_add{:}, name_pat{end});

    searchidx = cellfun(@(x) numel(regexp(x, name_pat3)), new_vname, 'UniformOutput', true);
    sigrltidx = find(searchidx >= numel(name_pat)-1);
    if isempty(sigrltidx)
        fprintf('No result can be visulized. TAT\n')
    end
    % prepare data for mediation
    symptom_idx = ~cellfun('isempty', regexp(new_vname,strcat('^',name_pat2_spt(name_pat2_idx),'$')));
    symptom_vname = name_pat2_spt(name_pat2_idx);
    symptom_data = table2array(T(:,symptom_idx));
    EEG_vname = new_vname(sigrltidx);
    EEG_data = table2array(T(:,sigrltidx));
    fmri_idx = ~cellfun('isempty', regexp(new_vname,name_pat2_spt(~cellfun('isempty', regexp(name_pat2_spt,'_')))));
    fmri_vname = new_vname(fmri_idx);
    fmri_data = table2array(T(:,fmri_idx));
    % run the mediation and plot
    cd(fig_path)
    for ss = 1:numel(symptom_vname)
        for ee = 1:numel(EEG_vname)
            for ff = 1:numel(fmri_vname)
                [paths, stats] = mediation(EEG_data(:,ee), symptom_data(:,ss), fmri_data(:,ff), 'boot', ...
                    'names', {regexprep(EEG_vname{ee}, '_', ' '), regexprep(symptom_vname{ss}(1:end-1), '_', ' '), ...
                    regexprep(fmri_vname{ff}, '_', ' ')}, 'doCIs','verbose', 'plots','dosave');% symptom need to delete the '$'
                %                 h = figure(2);
                %                 figure(h)
                fig_dir = dir(fullfile(fig_path, '*.png'));
                fig_dir2 = dir(fullfile(fig_path, '*.txt'));
                fig_name = {fig_dir.name};
                txt_name = {fig_dir2.name};
                spc_fig = fig_name(contains(fig_name, {'Path_Diagram','Histogram_Plot','Mediation_Scatterplots'}));
                spc_txt = txt_name(contains(txt_name, {'Mediation_Output'}));
                savetype = {'Path','Hist','Scat','Txt'};
                cellfun(@(x,y) eval(['!rename' , [',',x], [',',y,'_',EEG_vname{ee},'_',fmri_vname{ff},'.png' ]]),spc_fig, savetype(1:end-1));
                cellfun(@(x,y) eval(['!rename' , [',',x], [',',y,'_',EEG_vname{ee},'_',fmri_vname{ff},'.txt' ]]),spc_txt, savetype(end));
            end
        end
    end
    % clear the useless file
    close all;
    fig_dir3 = dir(fullfile(fig_path));
    filename = {fig_dir3.name};
    clearname = filename(contains(filename, {'Path_Diagram.','Histogram_Plot.','Mediation_Scatterplots.','Mediation_Output.'}));
    cellfun(@(x) delete(fullfile(fig_path, x)), clearname)
    
    fprintf('\nVisualization done.')
end
