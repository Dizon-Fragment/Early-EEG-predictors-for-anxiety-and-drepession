%% Initialization
clear; close all; clc;

%%% ---------------------- 1. Main Parameters ------------------------
main_path = '/data/home/EEG001/Longitudinal/'; 
result_plot_root = '/data/home/EEG001/Longitudinal/ForPlot/';
folders = {'csd','norm'}; % Data type for analysis

% --- Parameters for loading CONJUNCTION MASKS ---
Fre_file_masks = {'alpha','beta1'};
p_thresh_masks = [0.05,0.01];
year_data_masks = {'fz75', 'fz95', 'fz97'}; % These define WHICH mask to load
% opname_masks = 'Conjunction_result_groupedge_edge_resetyear';
opname_masks = 'Conjunction_result_multiplePlus_edge_resetyear';

% --- Parameters for loading LONGITUDINAL EEG DATA ---
longitudinal_years = {'2015', '2017', '2019'};
num_longitudinal_years = length(longitudinal_years);
idx_conn = {'coh','plv'}; % Connectivity measures to combine

% --- Parameters for ANALYSIS & SAVING ---
behaviors_to_analyze = {'SAS', 'SDS'};
output_foldername = 'Lateralization_Analysis_v5_mp_con_rainc';
perm_type = 'PERMconjunction';
% norm for PERMnorm, con for PERMconjunction
LME_Results = struct();

% --- Parameters for PLOTTING ---
close_figures_after_saving = true; % Set to false to keep figures open

%% Create result directories
output_path = fullfile(result_plot_root, output_foldername);
if ~exist(output_path, 'dir'), mkdir(output_path); end
for dd = 1:length(folders)
    output_dir_cond = fullfile(output_path, folders{dd});
    if ~exist(output_dir_cond, 'dir'), mkdir(output_dir_cond); end
end
warning('off');

%% ------------------- Main Analysis Loops -------------------

% Loop through data types ('csd', etc.)
for dd = 1:length(folders)
    current_folder = folders{dd};
    
    % Loop through frequency bands ('alpha', 'beta1')
    for ff = 1:length(Fre_file_masks)
        current_freq = Fre_file_masks{ff};
        
        fprintf('\n========================================================================\n');
        fprintf('Analyzing: FOLDER=%s | FREQ=%s\n', current_folder, current_freq);
        
        % --- PHASE 1: PRE-LOAD ALL LONGITUDINAL EEG DATA for this frequency band ---
        % This is efficient as we only load this large data once per frequency.
        fprintf('  Phase 1: Pre-loading longitudinal EEG data for %s...\n', current_freq);
        all_years_eeg_data = cell(1, num_longitudinal_years);
        try
            for yy_long = 1:num_longitudinal_years
                year_str = longitudinal_years{yy_long};
                temp_mats = cell(1, length(idx_conn));
                for ii = 1:length(idx_conn)
                    input_dir = fullfile(main_path, 'datafromMac', 'result', current_folder, 'CPM');
                    load(fullfile(input_dir, [idx_conn{ii},'_', current_freq, '_', year_str,'_','ForCPM.mat']));
                    var_name_to_eval = ['sub_', idx_conn{ii}, '_', current_freq, '_', year_str];
                    temp_mats{ii} = eval(var_name_to_eval);
                    clear(var_name_to_eval);
                end
                all_years_eeg_data{yy_long} = temp_mats;
            end
        catch ME
            fprintf('    ERROR: Failed to load longitudinal data for %s. Skipping this frequency band. Error: %s\n', current_freq, ME.message);
            continue; % Skip to the next frequency
        end
        
        % Now, loop through the conditions that define the conjunction mask
        for tt = 1:length(p_thresh_masks)
            current_p_thresh = p_thresh_masks(tt);
            p_thresh_name = ['p' strrep(num2str(current_p_thresh), '.', '')];

            for yy_mask = 1:length(year_data_masks)
                current_mask_year = year_data_masks{yy_mask};

                for behav_idx = 1:length(behaviors_to_analyze)
                    behavior_name = behaviors_to_analyze{behav_idx};
                    
                    analysis_condition_name = sprintf('%s_%s_%s_%s_%s', current_folder, current_freq, p_thresh_name, current_mask_year, behavior_name);
                    fprintf('\n--- Processing Condition: %s ---\n', analysis_condition_name);
                    
                    % --- PHASE 2: LOAD THE PRE-COMPUTED CONJUNCTION MASK ---
                    fprintf('  Phase 2: Loading conjunction mask...\n');
                    mask_input_dir = fullfile(main_path, 'datafromMac','result',  current_folder, opname_masks);
                    mask_file = fullfile(mask_input_dir, ['ind_',perm_type,'_area_', current_freq, '_', current_mask_year, '_', num2str(current_p_thresh), '_', opname_masks, '.mat']);
%                     mask_file = fullfile(mask_input_dir, ['ind_PERMconjunction_area_', current_freq, '_', current_mask_year, '_', num2str(current_p_thresh), '_', opname_masks, '.mat']);
                    
                    if ~exist(mask_file, 'file')
                        fprintf('    WARNING: Mask file not found, skipping. Path: %s\n', mask_file);
                        continue; % Skip to next behavior if mask file doesn't exist
                    end
                    
                    try
                        mask_data = load(mask_file);
                        res = mask_data.Extractedge_perm.(behavior_name);
                        
                        % According to your extraction script comments, these are the final masks
                        % res{29} is pos_matri_1 and res{30} is neg_matri_1
                        predictive_mask_pos = res{29};
                        predictive_mask_neg = res{30};
                        predictive_mask_glm = predictive_mask_pos | predictive_mask_neg;
                        
                        % Check if masks are empty
                        if ~any(predictive_mask_glm(:))
                           fprintf('    INFO: The loaded conjunction mask is empty. LI values will be zero.\n');
                        end

                    catch ME
                        fprintf('    WARNING: Could not load or extract data from mask file. Skipping. Error: %s\n', ME.message);
                        continue;
                    end
                    
                    % --- PHASE 3: CALCULATE LI TRAJECTORY USING THE LOADED MASK ---
                    fprintf('  Phase 3: Calculating LI trajectories...\n');
                    num_subjects = size(all_years_eeg_data{1}{1}, 3);
                    li_values_pos = zeros(num_subjects, num_longitudinal_years);
                    li_values_neg = zeros(num_subjects, num_longitudinal_years);
                    li_values_glm = zeros(num_subjects, num_longitudinal_years);
                    
                    left_nodes = [1,3,4,8,9,13,14,18,19,23,24,28];
                    right_nodes = [2,6,7,11,12,16,17,21,22,26,27,30];

                    for yy_long = 1:num_longitudinal_years
                        combined_conn_mat = abs(all_years_eeg_data{yy_long}{1}) + abs(all_years_eeg_data{yy_long}{2});
                        for subj = 1:num_subjects
                            subj_mat = combined_conn_mat(:, :, subj);
                            
                            % Apply masks and calculate LI
                            R_pos = sum(sum(subj_mat(right_nodes, right_nodes) .* predictive_mask_pos(right_nodes, right_nodes)));
                            L_pos = sum(sum(subj_mat(left_nodes, left_nodes) .* predictive_mask_pos(left_nodes, left_nodes)));
                            li_values_pos(subj, yy_long) = (R_pos - L_pos) / (R_pos + L_pos + 1e-9);
                            
                            R_neg = sum(sum(subj_mat(right_nodes, right_nodes) .* predictive_mask_neg(right_nodes, right_nodes)));
                            L_neg = sum(sum(subj_mat(left_nodes, left_nodes) .* predictive_mask_neg(left_nodes, left_nodes)));
                            li_values_neg(subj, yy_long) = (R_neg - L_neg) / (R_neg + L_neg + 1e-9);
                            
                            R_glm = sum(sum(subj_mat(right_nodes, right_nodes) .* predictive_mask_glm(right_nodes, right_nodes)));
                            L_glm = sum(sum(subj_mat(left_nodes, left_nodes) .* predictive_mask_glm(left_nodes, left_nodes)));
                            li_values_glm(subj, yy_long) = (R_glm - L_glm) / (R_glm + L_glm + 1e-9);
                        end
                    end
                    
                    % --- PHASE 4: OUTLIER REMOVAL & LME MODELING ---
                    analysis_types = {'Positive_Network_LI', 'Negative_Network_LI', 'GLM_Network_LI'};
                    all_li_data = {li_values_pos, li_values_neg, li_values_glm};
                    
                    for analysis_idx = 1:length(analysis_types)
                        analysis_name = analysis_types{analysis_idx};
                        current_li_data = all_li_data{analysis_idx};
                        
                        % Outlier removal using Z-score
                        z_thresh = 5;%2.5
                        outlier_mask = false(size(current_li_data));
                        for yy_long = 1:num_longitudinal_years
                            year_data = current_li_data(:, yy_long);
                            mean_val = mean(year_data, 'omitnan');
                            std_val = std(year_data, 'omitnan');
                            if std_val > 0
                                z_scores = abs((year_data - mean_val) / std_val);
                                outlier_mask(:, yy_long) = z_scores > z_thresh;
                            end
                        end
                        unique_outlier_indices = find(any(outlier_mask, 2));
                        
                        current_li_data_clean = current_li_data;
                        if ~isempty(unique_outlier_indices)
                            fprintf('    Found and removed %d outliers for %s\n', length(unique_outlier_indices), analysis_name);
                            current_li_data_clean(unique_outlier_indices, :) = NaN;
                        end

                        % Create table for LME
                        SubjectID = repelem((1:num_subjects)', num_longitudinal_years, 1);
                        Time = repmat([0; 2; 4], num_subjects, 1);
                        LI = reshape(current_li_data_clean, [], 1);
                        li_table = table(SubjectID, Time, LI);
                        li_table = rmmissing(li_table);
                        
                        % =================================================================
                        % --- FINAL "HYBRID" PLOTTING MODULE (Corrected for fliplr) ---
                        % This uses your original raincloud function and adds a mean line 
                        % with CORRECTLY MATCHED coordinates.
                        % =================================================================
                        
                        % --- 1. Prepare Data for your rainCloudsFor_longitudinal function ---
                        time_points_unique = unique(li_table.Time);
                        dataInput = cell(1, length(time_points_unique));
                        dataPlotName = cell(1, length(time_points_unique));
                        
                        for t_idx = 1:length(time_points_unique)
                            time_point = time_points_unique(t_idx);
                            dataInput{t_idx} = li_table.LI(li_table.Time == time_point);
                            dataPlotName{t_idx} = longitudinal_years{t_idx};
                        end
                        
                        % --- 2. Call YOUR rainCloudsFor_longitudinal function ---
                        % It will create the base plot with its own styling and reversed x-axis.
                        DB_fillrate = 0.068; DotSize = 60; markerEdgeAlpha = 0.2; 
                        markerFaceAlpha = 0.3; halfline = 0.05;
                        
                        fig = rainCloudsFor_longitudinal_v2(dataInput, dataPlotName, DB_fillrate, DotSize,...
                            markerEdgeAlpha, markerFaceAlpha, halfline, 0, 0, 0, 0, 0, ''); 
                        
                        % --- 3. Add the Group Mean Line with CORRECTED coordinates ---
                        ax = get(fig, 'CurrentAxes');
                        hold(ax, 'on');
                        
                        % Calculate mean values in the standard order [mean_2015, mean_2017, mean_2019]
                        mean_li = grpstats(li_table.LI, li_table.Time, 'mean');
                        
                        % !!! THE CRITICAL FIX !!!
                        % Your function plots using fliplr, so the x-axis locations are
                        % 2019 -> x=1, 2017 -> x=2, 2015 -> x=3.
                        % We MUST plot our mean line using these same reversed x-coordinates.
                        x_axis_for_mean_line = 1:length(time_points_unique);
                        
                        plot(ax, x_axis_for_mean_line, mean_li, '-o', ...
                            'Color', [0.2 0.2 0.2], 'LineWidth', 2, ...
                            'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.2 0.2 0.2], ...
                            'MarkerSize', 10, 'DisplayName', 'Group Mean');
                        
                        hold(ax, 'off');
                        
                        % --- 4. Final Touches (Update Title and Legend) ---
                        title_str = sprintf('%s Trajectory\n(%s, %s, mask: %s, %s)', ...
                                      strrep(analysis_name,'_',' '), current_folder, current_freq, current_mask_year, behavior_name);
                        title(ax, title_str, 'FontSize', 14, 'FontWeight', 'normal');
                        
%                         % Recreate the legend to include the new 'Group Mean'
                        ph = get(ax, 'children');
                        patch_handles = findobj(ph, 'Type', 'Patch');
                        mean_line_handle = findobj(ph, 'DisplayName', 'Group Mean');
                        
                        % The patch handles are also in reversed order, so we flip them back
%                         legend_handles = [flipud(patch_handles); mean_line_handle];
                        legend_handles = [(patch_handles); mean_line_handle];
%                         legend_names = [dataPlotName, {'Group Mean'}];
                        legend_names = [{'Group Mean'},dataPlotName];
%                         legend(ax, legend_handles, legend_names, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');
                        legend('Location', 'northeast', 'FontSize', 12);
                        
                        
                        % --- 5. Save in Multiple High-Quality Formats ---
                        plot_filename_base = fullfile(output_path, current_folder, sprintf('Plot_%s_%s', analysis_condition_name, analysis_name));
                        saveas(fig, [plot_filename_base, '.svg']);
                        print(fig, [plot_filename_base, '.png'], '-dpng', '-r300');
                        
                        if close_figures_after_saving, close(fig); end
                        
                        % LME Modeling
                        if height(li_table) < num_subjects || isempty(li_table) || all(li_table.LI == 0)
                            LME_Results.(current_folder).(current_freq).(p_thresh_name).(current_mask_year).(behavior_name).(analysis_name).model = 'Skipped_InsufficientData';
                            continue;
                        end
                        try
                            lme_model = fitlme(li_table, 'LI ~ 1 + Time + (1 + Time | SubjectID)');
                            LME_Results.(current_folder).(current_freq).(p_thresh_name).(current_mask_year).(behavior_name).(analysis_name).model = lme_model;
                        catch ME
                            fprintf('    LME model failed for %s. Error: %s\n', analysis_condition_name, ME.message);
                            LME_Results.(current_folder).(current_freq).(p_thresh_name).(current_mask_year).(behavior_name).(analysis_name).model = 'Failed_FitError';
                        end
                        
                        % two years
                        
                        % --- Exploratory Analysis: Trend from 2017-2019 ---
                        fprintf('    Running exploratory LME for 2017-2019 trend...\n');
                        li_table_late = li_table(li_table.Time >= 2, :); % Filter for Time points 2 (2017) and 4 (2019)
                        
                        % Recode Time for interpretability: 2017 becomes 0, 2019 becomes 2
                        li_table_late.Time = li_table_late.Time - 2;
                        
                        if height(li_table_late) > 0
                                lme_model_late = fitlme(li_table_late, 'LI ~ 1 + Time + (1 + Time | SubjectID)');
                                % ???????????????????
                                LME_Results.(current_folder).(current_freq).(p_thresh_name).(current_mask_year).(behavior_name).(analysis_name).exploratory_model_1719 = lme_model_late;
                                disp('    Exploratory LME (2017-2019) Results:');
                                disp(lme_model_late.Coefficients);

                        end
                    end % End of analysis_type loop
                end % End of behavior loop
            end % End of mask_year loop
        end % End of p_thresh loop
        
        % --- ROBUST SAVING: Save after each frequency band is fully processed ---
        fprintf('\n>>> CHECKPOINT SAVE after completing frequency band: %s <<<\n', current_freq);
        save(fullfile(output_path, 'LME_Conjunction_Results.mat'), 'LME_Results');

    end % End of frequency loop
end % End of folder loop

fprintf('\n\n SCRIPT HAS COMPLETED SUCCESSFULLY! \n The variable LME_Results has been saved.\n\n');

%% ------------------- Summary Script Section (Updated for Exploratory Analysis) -------------------
clear; clc; close all;

% --- Parameters for summary ---
result_plot_root = '/data/home/EEG001/Longitudinal/ForPlot/';
output_foldername = 'Lateralization_Analysis_v5_mp_con_rainc'; % Ensure this is correct
lme_results_file = fullfile(result_plot_root, output_foldername, 'LME_Conjunction_Results.mat');

% --- Load Data ---
fprintf('Loading LME results from: %s\n', lme_results_file);
if ~exist(lme_results_file, 'file'), error('LME results file not found.'); end
load(lme_results_file);

% --- Initialize Table with a new column for the analysis window ---
results_data = {};
header = {'DataType','Frequency','p_Threshold','Mask_Year','Behavior','Analysis_Type', ...
          'AnalysisWindow', ... % <<< NEW COLUMN HERE
          'Effect','Estimate','SE','tStat','DF','pValue','Lower_CI','Upper_CI'};
results_data(1, :) = header;

% --- Traverse the DEEP results structure ---
fprintf('Extracting and consolidating results...\n');
all_folders = fieldnames(LME_Results);
for d_idx = 1:length(all_folders)
    folder_name = all_folders{d_idx}; if isempty(LME_Results.(folder_name)), continue; end
    all_freqs = fieldnames(LME_Results.(folder_name));
    for f_idx = 1:length(all_freqs)
        freq_name = all_freqs{f_idx}; if isempty(LME_Results.(folder_name).(freq_name)), continue; end
        all_p_thresh = fieldnames(LME_Results.(folder_name).(freq_name));
        for p_idx = 1:length(all_p_thresh)
            p_name = all_p_thresh{p_idx}; if isempty(LME_Results.(folder_name).(freq_name).(p_name)), continue; end
            all_mask_years = fieldnames(LME_Results.(folder_name).(freq_name).(p_name));
            for my_idx = 1:length(all_mask_years)
                mask_year_name = all_mask_years{my_idx}; if isempty(LME_Results.(folder_name).(freq_name).(p_name).(mask_year_name)), continue; end
                all_behaviors = fieldnames(LME_Results.(folder_name).(freq_name).(p_name).(mask_year_name));
                for b_idx = 1:length(all_behaviors)
                    behav_name = all_behaviors{b_idx}; if isempty(LME_Results.(folder_name).(freq_name).(p_name).(mask_year_name).(behav_name)), continue; end
                    all_analysis_types = fieldnames(LME_Results.(folder_name).(freq_name).(p_name).(mask_year_name).(behav_name));
                    for a_idx = 1:length(all_analysis_types)
                        analysis_name = all_analysis_types{a_idx};
                        
                        % Get the parent struct that contains the models
                        parent_struct = LME_Results.(folder_name).(freq_name).(p_name).(mask_year_name).(behav_name).(analysis_name);
                        
                        % --- 1. Process the main model (2015-2019) ---
                        if isfield(parent_struct, 'model') && isstruct(parent_struct.model)
                            fixed_effects = parent_struct.model.Coefficients;
                            for effect_idx = 1:height(fixed_effects)
                                effect_row = fixed_effects(effect_idx, :);
                                new_row = {folder_name, freq_name, p_name, mask_year_name, behav_name, analysis_name, ...
                                           '2015-2019', ... % <<< Add window label
                                           effect_row.Name{1}, effect_row.Estimate, effect_row.SE, effect_row.tStat, ...
                                           effect_row.DF, effect_row.pValue, effect_row.Lower, effect_row.Upper};
                                results_data(end+1, :) = new_row;
                            end
                        end
                        
                        % --- 2. Process the exploratory model (2017-2019) ---
                        if isfield(parent_struct, 'exploratory_model_1719') && isstruct(parent_struct.exploratory_model_1719)
                            fixed_effects = parent_struct.exploratory_model_1719.Coefficients;
                            for effect_idx = 1:height(fixed_effects)
                                effect_row = fixed_effects(effect_idx, :);
                                new_row = {folder_name, freq_name, p_name, mask_year_name, behav_name, analysis_name, ...
                                           '2017-2019', ... % <<< Add window label
                                           effect_row.Name{1}, effect_row.Estimate, effect_row.SE, effect_row.tStat, ...
                                           effect_row.DF, effect_row.pValue, effect_row.Lower, effect_row.Upper};
                                results_data(end+1, :) = new_row;
                            end
                        end
                        
                    end
                end
            end
        end
    end
end

% --- Write to CSV ---
output_csv_path = fullfile(result_plot_root, output_foldername, 'lme_summary_conjunction_with_exploratory.csv');
if size(results_data, 1) < 2
    fprintf('\n\nNo valid LME model results were found to summarize.\n');
else
    results_table = cell2table(results_data(2:end,:), 'VariableNames', header);
    % Sort by the new column as well for clarity
    results_table = sortrows(results_table, {'DataType', 'Frequency', 'p_Threshold', 'Mask_Year', 'Behavior', 'Analysis_Type', 'AnalysisWindow', 'Effect'});
    writetable(results_table, output_csv_path);
    fprintf('\n\nSuccessfully consolidated all LME results into the file:\n%s\n', output_csv_path);
end
disp('--- LME summary script finished. ---');
