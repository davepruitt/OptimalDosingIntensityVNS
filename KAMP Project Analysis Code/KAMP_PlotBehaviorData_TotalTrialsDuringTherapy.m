function KAMP_PlotBehaviorData_TotalTrialsDuringTherapy ( kamp_project_data, parameter_name )

    if (nargin < 2)
        parameter_name = 'total_trials_per_day';
    end
    
    [rat_exclusion_list, ordered_rat_exclusion_list] = KAMP_GetExclusionList(kamp_project_data);
    
    %Black out the data for excluded rats so it doesn't get factored into
    %the calculations for this plot
%     for r = 1:length(rat_exclusion_list)
%         this_rat_name = rat_exclusion_list{r};
%         this_rat_row_idx = find(strcmpi(kamp_project_data.rat_names, this_rat_name), 1, 'first');
%         if (~isempty(this_rat_row_idx))
%             kamp_project_data.(parameter_name)(this_rat_row_idx, :) = NaN;
%         end
%     end
    

    % Grab the data and sum across days
    trials_data = nansum(kamp_project_data.(parameter_name)(:, 9:33), 2);
    
    figure;
    hold on;
    colors = colormap(lines);
    
    unique_groups = {'No VNS', 'VNS 0.4 mA', 'VNS 0.8 mA', 'VNS 1.6 mA'};
    
    all_means = [];
    all_sem = [];
    x_vals = [];
    y_vals = [];
    
    anova_data = struct('group', {}, 'data', {});
    
    %Iterate over each rat
    for g = 1:length(unique_groups)
        this_group_name = unique_groups{g};
        this_group_indices = find(strcmpi(kamp_project_data.groups, this_group_name));
        
        color_to_plot = colors(1, :);
        if (strcmpi(this_group_name, 'VNS 0.8 mA'))
            color_to_plot = colors(2, :);
        elseif (strcmpi(this_group_name, 'VNS 0.4 mA'))
            color_to_plot = colors(3, :);
        elseif (strcmpi(this_group_name, 'VNS 1.6 mA'))
            color_to_plot = colors(4, :);
        elseif (strcmpi(this_group_name, 'Unimpaired'))
            color_to_plot = colors(6, :);
        elseif (strcmpi(this_group_name, 'Not yet determined'))
            color_to_plot = colors(5, :);
        end
        
        this_group_exclusion_list = ordered_rat_exclusion_list(this_group_indices);
        this_group_indices = this_group_indices(~this_group_exclusion_list);
        this_group_data = trials_data(this_group_indices, :);
        this_group_mean = nanmean(this_group_data, 1);
        if (size(this_group_data, 1) > 1)
            this_group_sem = nanstd(this_group_data, 1) / sqrt(size(this_group_data, 1));
        else
            this_group_sem = zeros(size(this_group_data, 1), size(this_group_data, 2));
        end
        
        all_means = [all_means this_group_mean];
        all_sem = [all_sem this_group_sem];
        
        this_group_x_vals = generate_xvals_v4(this_group_data, g, 0.1, 400, 0.1);
        for r = 1:length(this_group_data)
            x_vals = [x_vals this_group_x_vals(r)];
            y_vals = [y_vals this_group_data(r)];
            anova_data(end+1) = struct('group', g, 'data', this_group_data(r));
        end
        
    end
    
    b = bar(all_means);
    b.FaceColor = 'flat';
    b.CData(1, :) = colors(1, :);
    b.CData(2, :) = colors(3, :);
    b.CData(3, :) = colors(2, :);
    b.CData(4, :) = colors(4, :);
    
    errorbar(1:4, all_means, all_sem, all_sem, 'Marker', 'none', 'LineStyle', 'none', 'Color', 'k', ...
        'MarkerFaceColor', 'k');
    
    plot(x_vals, y_vals, 'Marker', 'o', 'MarkerFaceColor', 'w', ...
        'Color', 'k', 'LineStyle', 'none', 'MarkerSize', 5);
    
    
    set(gca, 'XTick', 1:4);
    set(gca, 'XTickLabel', unique_groups);
    ylabel('Trials performed during therapy');
    
    
    dependent_variable = [anova_data.data];
    independent_variable = {[anova_data.group]};
    
    [p, tbl, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    
    disp('Percent recovery table:');
    disp(tbl);
    
    
end















































