%% Plot the result
clear; close all; clc;
m3_path = '/data/home/EEG001/tools/MediationToolbox-master';
m3_supt_path = '/data/home/EEG001/tools/CanlabCore-master';

restoredefaultpath;
addpath(genpath(m3_path));
addpath(genpath(m3_supt_path));
% addpath(genpath(spm12_path));

%% Mediation
clear; close all; 
%%% ----------------- Parameter for finding ---------------------
main_path = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data'; 
result_plot = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data/picresult';
montage1 = {'csd'};
Fre_file = {'alpha','beta1'};
Fz_year = {'2017', '2019', 'fz97'};
thresh = [0.05]; 
Tname = {'SAS',  'SDS'};
new_Fre = { 'a', 'b1'};
new_Tname = {'A',  'D'};
new_thre = [5];
str_thre = {'p5'};
pic_path = fullfile( result_plot, 'Mediation');
mkdir( pic_path)
data_path = fullfile( main_path, 'ForMediation', 'For2idx');
xlsx_path = data_path;
addpath(data_path)
fig_path = fullfile( result_plot,  'mediationfig');
fmri_path = fullfile( main_path, 'longitudinalData', 'FC'); 
symptom_path = fullfile( main_path,  'subjset_symptom.mat');
mkdir(fig_path)
%%% -------------------------------------------------------------


%% Run
result_dir = dir(fullfile(data_path));
result_folder = {result_dir(cellfun(@(x) x==1,{result_dir.isdir})).name};
result_folder(1:2) = []; % delete the '.' and '..'
for ff = 1:numel(Fre_file)
    for yy = 1:numel(Fz_year)
        xlsxname = {result_dir(contains({result_dir.name},'xlsx')).name};
        if contains(Fz_year{yy}, 'fz')
            tmp_path = fullfile( data_path, result_folder{contains(result_folder, 'fz')});
            tmp_xlsxname = xlsxname(contains(xlsxname, 'Fz'));
        else
            tmp_path = fullfile( data_path, result_folder{~contains(result_folder, {'fz','fig'})});
            tmp_xlsxname = xlsxname(~contains(xlsxname, 'Fz'));
        end
        
        for tt = 1:numel(new_thre)
            for ss = 1:numel(new_Tname)
                clear -except [^tmp_xlsxname, tmp_path,result_folder,fig_path,data_path,color_path,result_plot, montage1, Fre_file, Fz_year, new_thre, Tname, new_Fre, new_Tname, str_thre, pic_path]
                tmp_dir = dir( fullfile( tmp_path, '*.mat')); % load the name of the file
                tmp_foldername = {tmp_dir.name};
                tmp_x = [montage1{1}, '_', Fre_file{ff}, '_', Fz_year{yy}(end-1:end), '_p',...
                    num2str(new_thre(tt)), '_', regexprep(Tname{ss},'_', '')];% pattern name
                fprintf('\nLoading %s ...', tmp_x)
                tmp_subset = tmp_foldername(contains( tmp_foldername, tmp_x));
                spt_vname = cellfun(@(x) strsplit(x, '_'), tmp_subset, 'UniformOutput', false);
                for bb = 1:numel(tmp_subset)
                    tmp_filename = tmp_subset{bb};
                    load(fullfile( tmp_path, tmp_filename))
                    % 'paths', 'stats', 'mri_data', 'mri_colname', 'symptom_data', 'symptom_colname',...
                    %                 'eegconnetome', 'xlsx_save', 'colidx'
                    T = readtable(fullfile(data_path, xlsx_save));
                    vname = T.Properties.VariableNames;
                    spec_mri_idx = strfind(tmp_filename, '_');
                    spec_mri = tmp_filename(spec_mri_idx(5)+1:spec_mri_idx(7)-1);
                    % run mediation and plot
                    cd(fig_path)
                    [paths, stats] = mediation_dgz(table2array(T(:,colidx)), table2array(symptom_data), eval(['mri_data.', spec_mri]), 'boot', ...
                        'names', {regexprep(vname{colidx}, '_', ' '), regexprep(symptom_colname{:}, '_', ' '), ...
                        regexprep(spec_mri, '_', ' ')},  'doCIs','verbose', 'plots','dosave');% symptom need to delete the '$'
                    % change the filenames
                    fig_dir = dir(fullfile(fig_path, '*.png'));
                    fig_dir3 = dir(fullfile(fig_path, '*.svg'));
                    fig_dir2 = dir(fullfile(fig_path, '*.txt'));
                    fig_name = {fig_dir.name};
                    txt_name = {fig_dir2.name};
                    svg_name = {fig_dir3.name};
                    spc_fig = fig_name(contains(fig_name, {'Histogram_Plot','Path_Diagram','Mediation_Scatterplots'}));
                    spc_svg = svg_name(contains(svg_name, {'Histogram_Plot','Path_Diagram','Mediation_Scatterplots'}));
                    spc_txt = txt_name(contains(txt_name, {'Mediation_Output'}));
                    %                     savetype = {'Path','Hist','Scat','Txt'};
                    savetype = {'Hist','Scat','Path','Txt'};
                    if ispc
                        cellfun(@(x,y) eval(['!rename' , [',',x], [',',y,'_',tmp_filename(1:end-4),'.svg' ]]),spc_svg, savetype(1:end-1));
                        cellfun(@(x,y) eval(['!rename' , [',',x], [',',y,'_',tmp_filename(1:end-4),'.txt' ]]),spc_txt, savetype(end));
                        cellfun(@(x,y) eval(['!rename' , [',',x], [',',y,'_',tmp_filename(1:end-4),'.png' ]]),spc_fig, savetype(1:end-1));
                    elseif isunix
                        cellfun(@(x,y) eval(['!mv' , [' ',x], [' ',y,'_',tmp_filename(1:end-4),'.svg' ]]),spc_svg, savetype(1:end-1));
                        cellfun(@(x,y) eval(['!mv' , [' ',x], [' ',y,'_',tmp_filename(1:end-4),'.txt' ]]),spc_txt, savetype(end));
                        cellfun(@(x,y) eval(['!mv' , [' ',x], [' ',y,'_',tmp_filename(1:end-4),'.png' ]]),spc_fig, savetype(1:end-1));
                    end
                    % clear the useless file
                    close all;
                    fig_dir4 = dir(fullfile(fig_path));
                    filename = {fig_dir4.name};
                    clearname = filename(contains(filename, {'Path_Diagram.','Histogram_Plot.','Mediation_Scatterplots.','Mediation_Output.'}));
                    cellfun(@(x) delete(fullfile(fig_path, x)), clearname)
                    % save the stat data
                    stat_path = fullfile( fig_path, 'stat');mkdir(stat_path);
                    save(fullfile(stat_path, tmp_filename), 'stats', 'paths');
                end
            end
        end
    end
end
fprintf('\nVisualization saving done.')
        



