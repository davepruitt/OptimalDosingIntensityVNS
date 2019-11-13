function KAMP_PlotLesionSizeData ( kamp_project_data )

    parameter_name = 'stroke_lesion_size';
    
    [rat_exclusion_list, ordered_rat_exclusion_list] = KAMP_GetExclusionList(kamp_project_data);
    
    %Black out the data for excluded rats so it doesn't get factored into
    %the calculations for this plot
    for r = 1:length(rat_exclusion_list)
        this_rat_name = rat_exclusion_list{r};
        this_rat_row_idx = find(strcmpi(kamp_project_data.rat_names, this_rat_name), 1, 'first');
        if (~isempty(this_rat_row_idx))
            kamp_project_data.(parameter_name)(this_rat_row_idx) = NaN;
        end
    end
    

    % Grab the data and sum across days
    lesion_size_data = kamp_project_data.(parameter_name);
    
    figure;
    hold on;
    colors = colormap(lines);
    
    unique_groups = {'No VNS', 'VNS 0.4 mA', 'VNS 0.8 mA', 'VNS 1.6 mA'};
    
    anova_data = struct('group', {}, 'data', {});
    all_means = [];
    all_sem = [];
    x_vals = [];
    y_vals = [];
    
    all_included_indices = [];
    
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
        this_group_data = lesion_size_data(this_group_indices);
        this_group_data = this_group_data(~isnan(this_group_data));
        this_group_mean = nanmean(this_group_data, 2);
        if (size(this_group_data, 2) > 1)
            this_group_sem = nanstd(this_group_data) / sqrt(size(this_group_data, 2));
        else
            this_group_sem = zeros(size(this_group_data, 1), size(this_group_data, 2));
        end
        
        all_means = [all_means this_group_mean];
        all_sem = [all_sem this_group_sem];
        all_included_indices = [all_included_indices this_group_indices];
        
        disp([this_group_name ', n = ' num2str(length(this_group_data))]);
        for r = 1:length(this_group_data)
            anova_data(end+1) = struct('group', g, 'data', this_group_data(r));
        end
        
        this_group_x_vals = generate_xvals_v4(this_group_data, g, 0.1, 1, 0.1);
        for r = 1:length(this_group_data)
            x_vals = [x_vals this_group_x_vals(r)];
            y_vals = [y_vals this_group_data(r)];
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
    ylabel('Lesion size (mm^3)');
    
    dependent_variable = [anova_data.data];
    independent_variable = {[anova_data.group]};
    [p, tbl, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', ...
        'varnames', {'Group'});
    [c, m] = multcompare(stats, 'CType', 'bonferroni', 'display', 'off');
    tbl
    c
    
    %% Figure 2:
    
    figure;
    hold on;
    
    for r = 1:length(kamp_project_data.rat_names)
        plot(kamp_project_data.stroke_lesion_progression(r, :), 'Color', 'b');
    end
    
    %% Now display the animal names with the smallest, average, and largest lesions:
    
    included_lesion_size_data = lesion_size_data(all_included_indices);
    
    [min_val, min_idx] = nanmin(included_lesion_size_data);
    min_rat_name = kamp_project_data.rat_names{all_included_indices(min_idx)};
    
    [max_val, max_idx] = nanmax(included_lesion_size_data);
    max_rat_name = kamp_project_data.rat_names{all_included_indices(max_idx)};
    
    median_lesion_size = nanmedian(included_lesion_size_data);
    lesions_minus_median = included_lesion_size_data - median_lesion_size;
    [~, med_idx] = nanmin(abs(lesions_minus_median));
    med_val = included_lesion_size_data(med_idx);
    med_rat_name = kamp_project_data.rat_names{all_included_indices(med_idx)};
    
    disp(['Smallest lesion: ' min_rat_name ', ' num2str(min_val) ' mm^3' ]);
    disp(['Largest lesion: ' max_rat_name ', ' num2str(max_val) ' mm^3' ]);
    disp(['Median lesion: ' med_rat_name ', ' num2str(med_val) ' mm^3' ]);
    disp(['True median: ' num2str(median_lesion_size)]);
    
    
end















































