%% Plot the result
clear; close all; clc;
restoredefaultpath

%% Plot the permutation distribution
%%% ---------------------- parameters ------------------------
main_path = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes'; 
result_plot = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data/picresult';
eeglab_path = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/toolbox/eeglab14_0_0b';
pic_path = fullfile( result_plot, 'braintopo_plot');
mkdir( pic_path)
code_path = fullfile( result_plot, 'toolbox', 'connectome');
color_path = fullfile( result_plot, 'toolbox', 'slanCM');
addpath( genpath(eeglab_path))
addpath( genpath(code_path))
addpath( genpath( color_path))
eeglab; close all;

%% Specific parameters
symptom = {'SAS','SDS'};
year = {'2017', '2019', 'fz97'};
freq = {'alpha', 'beta1'};
% subedge = 'neg';
closepic = 0;
rgbcolor = [31,184,207];


%% load the Extractedge_perm
outpdir = fullfile( main_path, 'data','longitudinalData','csd', 'result','dataExtraction');
load(fullfile( outpdir, 'preLoo2idx_Extraction.mat'))
%%% perm_data perm_info corr_data corr_info
corr_forplot = struct();
for i = 1:numel(symptom)
    % extract the distribution
    for y = 1:numel(year)
        clear -except [^rgbcolor, corr_forplot, perm_data perm_info corr_data corr_info closepic, pic_path, symptom, year, freq, result_plot, main_path]
        catname = strcat(freq{i}, '_', year{y});
        idx1 = cellfun(@(x) contains(x, catname), perm_data(:,2), 'UniformOutput', true);
        tmp_corr = corr_data(idx1,:);
        tmp_perm = perm_data(idx1,:);
        idx2 = cellfun(@(x) strcmp(x, symptom{i}), perm_data(idx1,3), 'UniformOutput', true);
        subname = {'pos', 'neg', 'glm'};
        %%% find out the significant results of correlation
        spec_corr = tmp_corr(idx2,:);
        spec_corr_p = spec_corr{1,4};
        spec_corr_xy = spec_corr{1,7};
        idx3 = find( spec_corr_p <= 0.05);
        subname2 = subname(idx3);
        %%% find out the significant results of correlation
        spec_perm = tmp_perm(idx2,:);
        
        %%% Plot
        for s = 1:length(idx3)
            switch idx3(s)
                case 1
                    tmp_xy = spec_corr_xy{1,1};
                    tmp_matrix = zeros( 30,30);
                    tmp_matrix2 = spec_corr_xy{1,3};
                case 2
                    tmp_xy = spec_corr_xy{1,2};
                    tmp_matrix = zeros( 30,30);
                    tmp_matrix2 = spec_corr_xy{1,4};
                case 3
                    tmp_xy = unique([spec_corr_xy{1,1}; spec_corr_xy{1,2}], 'rows');
                    tmp_matrix = zeros( 30,30);
                    tmp_matrix2 = spec_corr_xy{1,3} + spec_corr_xy{1,4}; 
            end
            for z = 1:size(tmp_xy,1)
                eval(['tmp_matrix(tmp_xy(z, 1), tmp_xy(z, 2)) = tmp_matrix2(tmp_xy(z, 1), tmp_xy(z, 2));']);
            end
            
            % plot the connectome
            load(fullfile( main_path, 'data', 'chanlocs_demo.mat'))
%             gaf = figure;
%             scrsz = get(0, 'ScreenSize');
%             set( gaf, 'Position', scrsz);
            
            gaf = figure(1);
%             subplot(1,length(idx3),s)
            imagesc( tmp_matrix) % plot for each
            xlabel('channels')
            ylabel('channels')
            title( sprintf( [subname2{s}, ' connectome']), 'fontsize', 12)
            
            % save the picture
            picpath = fullfile(result_plot, 'connect_matrix');
            mkdir( picpath)
            saveas( gaf, fullfile( picpath, [freq{i},'_', year{y}, '_', symptom{i}, '_', subname2{s}, '.svg']), 'svg');
            save( fullfile( picpath, [freq{i},'_', year{y}, '_', symptom{i}, '_', subname2{s}, '.mat']), 'tmp_matrix', 'tmp_xy')
            
            
%             subplot(1,length(idx3),s)
            gbf = figure(2);
            sgt = sgtitle( sprintf( [freq{i},', ', year{y}, ', ', symptom{i}, ', ', subname2{s},...
                '\n p = ', num2str(spec_corr_p(idx3(s)))]), 'fontsize', 6, 'Color', 'Red');      
            ds_neg = struct();
            ds_neg.chanPairs = tmp_xy;
            for z = 1:size( tmp_xy, 1)
                ds_neg.connectStrength(z,1) = tmp_matrix(tmp_xy(z,1), tmp_xy(z,2));
            end
            %                     ds_neg.connectStrengthLimits = [1, 52];
            ds_neg.cx = [0 1];
%             ds_neg.EMARKERSIZE =  16;
            ds_neg.EMARKERSIZE =  25;
            ds_neg.add_size =  10;
            ds_neg.line_width = 5;
%             ds_neg.HEADRINGWIDTH = 0.03;
%             topoplot_connect_dgz( ds_neg, chanlocs_demo,slanCM('dense'))% dense, [133,100]
            topoplot_connect_dgz( ds_neg, chanlocs_demo,slanCM(133, 100))% dense, [133,100]
            clear ds.neg

            % save the picture
            picpath = fullfile(result_plot, 'pic', 'topoplot');
            mkdir( picpath)
            saveas( gbf, fullfile( picpath, [freq{i},'_', year{y}, '_', symptom{i}, '_', subname2{s}, '_new.tif']), 'tif');
            close all;
        end
    end
    
    if closepic == 1
        close all;
    end
end


