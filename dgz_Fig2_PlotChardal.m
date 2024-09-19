%% Plot the result
clear; close all; clc;

%% Plot the permutation distribution
%%% ---------------------- parameters ------------------------
main_path = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes'; 
result_plot = '/data/home/EEG001/Longitudinal/ForPlot/pub_codes/data/picresult';
pic_path = fullfile( result_plot, 'BubbleChord');
mkdir( pic_path)
code_path = fullfile( result_plot, 'toolbox', 'bubble');
color_path = fullfile( result_plot, 'toolbox', 'slanCM');
addpath( genpath(code_path))
addpath( genpath( color_path))



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
            ch_name = {chanlocs_demo.labels};
            Data1=tril(tmp_matrix)-tril(tmp_matrix)';
            pthresh = size(spec_corr{1,6},1) * 0.8;
            %
            nodeName = {'FP1', 'FP2', 'F3', 'F4',  'FZ', 'FCZ', 'CZ', 'FC3','FC4', 'FT7', 'FT8', 'T3', ...
                'T4', 'TP7','TP8', 'T5', 'T6', 'C3', 'CP3', 'P3', 'PZ', 'C4','CP4', 'P4', 'CPZ', 'O1', 'OZ', 'O2', 'F7', 'F8',};
            className={'Prefrontal',{'Dorsolateral'; 'prefrontal'},'Frontal',...
                'Temporal','Parietal','Occipital',{'Ventrolateral'; 'prefrontal'},''};
            nodeNum = [2,2,7,6,8,3, 2];
            nodeIdx = cellfun(@(x) find(strcmp(ch_name, x)), nodeName, 'UniformOutput', true);
            Class = [];
            for c = 1:numel(nodeNum)
                Class = [Class; c*ones(nodeNum(c), 1)];
            end
            Data = Data1(nodeIdx, nodeIdx);
            % add none bubbles
            oripos = 0;
            classmax = max(Class) + 1;
            Data2 = Data;
            nodeName2 = nodeName;
            Class2 = Class;
            addBubble = 1;
            for c = 1:numel(nodeNum)
                oripos  = oripos + nodeNum(c); 
                Data2 = [Data2(1:oripos,:); zeros(addBubble,size(Data2,2)); Data2(oripos+1:end,:)];
                Data2 = [Data2(:,1:oripos), zeros(size(Data2,1),addBubble), Data2(:,oripos+1:end)];
                nodeName2 = [nodeName2(1,1:oripos), repmat({''},1,addBubble), nodeName2(1,oripos+1:end)];
                Class2 = [Class2(1:oripos,:); classmax*ones(addBubble,size(Class2,2)); Class2(oripos+1:end,:)];
                oripos = oripos + addBubble;
            end
%             Class2(Class2==classmax,1) = 0; % set to print white
            
            % plot the bubble graph
            gaf = figure(1);
            BD=bubbleDigraph_dgz(Data2,Class2,'NodeName',nodeName2,'RClass',1.32,'RNode',1.1, 'ClassName',className,...
                'Pthresh', pthresh);
            BD=BD.draw();
            bubblesize([6,15])
            BD.setNodeLabel('FontName','Arial','Color',[0,0,0],'FontSize',4)
            BD.setClassLabel('FontName','Arial','Color',[181,84,137]./255,'FontSize',7)
%             BD.setBubble('MarkerFaceAlpha',.2,'MarkerEdgeColor',[0,0,0])
            BD.setBubbleColor(turbo(10))
            
            % save the picture
            picpath = fullfile(result_plot, 'pic', 'bubble');
            mkdir( picpath)
            saveas( gaf, fullfile( picpath, [freq{i},'_', year{y}, '_', symptom{i}, '_', subname2{s}, '.svg']), 'svg');
            save( fullfile( picpath, [freq{i},'_', year{y}, '_', symptom{i}, '_', subname2{s}, '.mat']), 'tmp_matrix', 'tmp_xy')
            close all;
        end
    end
    
    if closepic == 1
        close all;
    end
end


