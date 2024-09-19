%% Plot the result
clear; close all; clc;

%% Plot the permutation distribution
%%% ---------------------- parameters ------------------------
main_path = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes'; 
result_plot = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data/picresult';
pic_path = fullfile( result_plot, 'perm_plot');
mkdir( pic_path)


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
        idx3 = find( spec_corr_p <= 0.05);
        subname2 = subname(idx3);
        %%% find out the significant results of correlation
        spec_perm = tmp_perm(idx2,:);
        %%% Plot
        for s = 1:length(idx3)
            tmp2 = tmp_perm{idx2,6}(:,idx3(s));
            truepred = tmp_corr{idx2, 5}(1,idx3(s));
            pvalue = tmp_perm{idx2, 4}(1,idx3(s));
            eval(['corr_forplot.',[freq{i},'_', year{y}, '_', symptom{i}, '_', subname2{s}],...
                ' = {[spec_corr{1,6}(:,1), spec_corr{1,6}(:,idx3(s)+1)],',...
                'spec_corr{1,5}(:,idx3(s)), spec_corr{1,4}(:,idx3(s))};'])
            % plot distribution
%             plot_rgb = [31,184,207]; % light blue
            % plot_rgb = [109,126,182];
            % plot_rgb = [182,231,249];
            plot_rgb = rgbcolor;
            figure; hold on;
            h = histogram(tmp2, 20);
            h.FaceColor = 'none';
            h.EdgeColor = plot_rgb./255;
            h.LineWidth = 2;
            ylim=get(gca,'Ylim');
            plot([truepred truepred],ylim,'r','linewidth',2, 'linestyle','--');
            title({[freq{i},' ', year{y}, ' ', symptom{i}, ' ', subname2{s}]; ['p = ',num2str(pvalue)]});
            set(gcf, 'renderer', 'painters');
            saveas(gcf, fullfile(pic_path, [freq{i},'_', year{y}, '_', symptom{i}, '_', subname2{s}, '.svg']))
            % print(gcf, fullfile(pic_path, [freq,'_', year, '_', symptom, '_', subedge]), '-dmeta', '-r600')
            if closepic == 1
                close all;
            end
        end
    end
end
