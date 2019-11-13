function KAMP_PlotBehaviorData ( kamp_project_data, parameter_name )

    if (nargin < 2)
        parameter_name = 'hit_rate';
    end

    [rat_exclusion_list, ordered_rat_exclusion_list] = KAMP_GetExclusionList(kamp_project_data);
    
    %Black out the data for excluded rats so it doesn't get factored into
    %the calculations for this plot
    for r = 1:length(rat_exclusion_list)
        this_rat_name = rat_exclusion_list{r};
        this_rat_row_idx = find(strcmpi(kamp_project_data.rat_names, this_rat_name), 1, 'first');
        if (~isempty(this_rat_row_idx))
            kamp_project_data.(parameter_name)(this_rat_row_idx, :) = NaN;
        end
    end
    

    % Plot group data by week

    pre = nanmean(kamp_project_data.(parameter_name)(:, 1:3), 2);
    post = nanmean(kamp_project_data.(parameter_name)(:, 4:8), 2);
    wk1 = nanmean(kamp_project_data.(parameter_name)(:, 9:13), 2);
    wk2 = nanmean(kamp_project_data.(parameter_name)(:, 14:18), 2);
    wk3 = nanmean(kamp_project_data.(parameter_name)(:, 19:23), 2);
    wk4 = nanmean(kamp_project_data.(parameter_name)(:, 24:28), 2);
    wk5 = nanmean(kamp_project_data.(parameter_name)(:, 29:33), 2);
    wk6 = nanmean(kamp_project_data.(parameter_name)(:, 34:38), 2);
    plot_data_in_epochs = [pre post wk1 wk2 wk3 wk4 wk5 wk6];
    
    figure;
    hold on;
    colors = colormap(lines);
    
    %unique_groups = unique(kamp_project_data.groups);
    unique_groups = {'No VNS', 'VNS 0.4 mA', 'VNS 0.8 mA', 'VNS 1.6 mA'};
    
    legend_strings = {};
    legend_pieces = [];
    
    
    ptt_tail = 'right';
    ptt_ncomparisons = 6;
    ptt_alpha = 0.05 / ptt_ncomparisons;
    
    disp(' ');
    disp('EFFECT OF TIME IN EACH GROUP: POST THROUGH WEEK 6');
    disp('Repeated-measures ANOVA within a single group, looking only at time');
    disp('Followed up by paired t-tests, signed-rank tests, and sign-test within that group, comparing post to all therapy weeks');
    disp(['t-test parameters: tail = ' ptt_tail ', num comparisons = ' num2str(ptt_ncomparisons) ', alpha = ' num2str(ptt_alpha)]);
    
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
        this_group_data = plot_data_in_epochs(this_group_indices, :);
        this_group_mean = nanmean(this_group_data, 1);
        if (size(this_group_data, 1) > 1)
            this_group_sem = nanstd(this_group_data, 1) / sqrt(size(this_group_data, 1));
        else
            this_group_sem = zeros(size(this_group_data, 1), size(this_group_data, 2));
        end
        
        h = errorbar(1:length(this_group_mean), this_group_mean, this_group_sem, this_group_sem, ...
            'Marker', 'none', 'MarkerFaceColor', color_to_plot, 'Color', color_to_plot, ...
            'LineStyle', '-', 'LineWidth', 2);
        
        if (g == 3)
            legend_strings{end+1} = [this_group_name ' (n = 6)'];
        elseif (g == 1)
            legend_strings{end+1} = [this_group_name ' (n = 12)'];
        else
            legend_strings{end+1} = [this_group_name ' (n = ' num2str(sum(~this_group_exclusion_list)) ')'];
        end
        %legend_strings{end+1} = [this_group_name ' (n = ' num2str(sum(~this_group_exclusion_list)) ')'];
        legend_pieces(end+1) = h;
        
        
        t = table(this_group_data(:, 2), this_group_data(:, 3), this_group_data(:, 4), this_group_data(:, 5), ...
                  this_group_data(:, 6), this_group_data(:, 7), this_group_data(:, 8), ...
                  'VariableNames', {'Post', 'Wk1', 'Wk2', 'Wk3', 'Wk4', 'Wk5', 'Wk6'});
        timepoints_vector = 1:7;
        rm = fitrm(t, 'Post-Wk6~1', 'WithinDesign', timepoints_vector);
        ranova_tbl = ranova(rm);
        multcompare(rm, 'Time');
        ranova_pval = table2array(ranova_tbl(1, 5));
        
        disp(' ');
        disp(this_group_name);
        disp(['Repeated-measures ANOVA: p = ' num2str(ranova_pval)]);
        
        for t = 3:8
            post_data_cmp = this_group_data(:, 2);
            therapy_data_cmp = this_group_data(:, t);
            [p_signtest, ~] = signtest(therapy_data_cmp, post_data_cmp, 'tail', ptt_tail, 'alpha', ptt_alpha);
            [p_signrank, ~] = signrank(therapy_data_cmp, post_data_cmp, 'tail', ptt_tail, 'alpha', ptt_alpha);
            [h0, pval_ttest] = ttest(therapy_data_cmp, post_data_cmp, 'tail', ptt_tail, 'alpha', ptt_alpha);
            disp(['Post vs Wk' num2str(t - 2) ': (t-test) p = ' num2str(pval_ttest) ', (signrank) p = ' num2str(p_signrank) ', (signtest) p = ' num2str(p_signtest)]);
            
            if (pval_ttest < 0.05)
                plot(t, this_group_mean(t), 'Marker', 'o', 'MarkerFaceColor', color_to_plot, ...
                    'Color', color_to_plot, 'LineStyle', 'none', 'MarkerSize', 8);
            else
                plot(t, this_group_mean(t), 'Marker', 'o', 'MarkerFaceColor', [1 1 1], ...
                    'Color', color_to_plot, 'LineStyle', 'none', 'MarkerSize', 8);
            end
        end
        
    end
    
    %Set the ylimits
    ylim([0 100]);
    
    %Set the xlimits
    xlim([0.5 8.5]);
    
    epoch_names = {'Pre', 'Post', 'Wk1', 'Wk2', 'Wk3', 'Wk4', 'Wk5', 'Wk6'};
    set(gca, 'XTick', 1:length(epoch_names));
    set(gca, 'XTickLabel', epoch_names);
    set(gca, 'TickLength', [0 0]);
    
    if (strcmpi(parameter_name, 'hit_rate'))
        ylabel('Percent trials above 60 degrees');
    elseif (strcmpi(parameter_name, 'maximal_turn_angle'))
        ylabel('Maximal turn angle (degrees)');
    end
    
    %% Now let's do a lot of statistics
    
    %First, let's do a repeated-measures ANOVA
    rows_to_include = [];
    for r = 1:length(kamp_project_data.rat_names)
        this_rat_name = kamp_project_data.rat_names{r};
        this_rat_group = kamp_project_data.groups{r};
        is_rat_excluded = ~isempty(find(strcmpi(rat_exclusion_list, this_rat_name), 1, 'first'));
        is_group_included = ~isempty(find(strcmpi(unique_groups, this_rat_group), 1, 'first'));
        if (~is_rat_excluded && is_group_included)
            rows_to_include = [rows_to_include r];
        end
    end
    
    data_copy = kamp_project_data;
    data_copy.groups = data_copy.groups(rows_to_include);
    data_copy.(parameter_name) = data_copy.(parameter_name)(rows_to_include, :);
    
    groups = data_copy.groups';
    pre = nanmean(data_copy.(parameter_name)(:, 1:3), 2);
    post = nanmean(data_copy.(parameter_name)(:, 4:8), 2);
    wk1 = nanmean(data_copy.(parameter_name)(:, 9:13), 2);
    wk2 = nanmean(data_copy.(parameter_name)(:, 14:18), 2);
    wk3 = nanmean(data_copy.(parameter_name)(:, 19:23), 2);
    wk4 = nanmean(data_copy.(parameter_name)(:, 24:28), 2);
    wk5 = nanmean(data_copy.(parameter_name)(:, 29:33), 2);
    wk6 = nanmean(data_copy.(parameter_name)(:, 34:38), 2);
    
    t = table(groups, post, wk1, wk2, wk3, wk4, wk5, wk6, 'VariableNames', {'Group', 'Post', 'Wk1', 'Wk2', 'Wk3', 'Wk4', 'Wk5', 'Wk6'});
    timepoints_vector = 1:7;
    rm = fitrm(t, 'Post-Wk6~Group', 'WithinDesign', timepoints_vector);
    ranova_tbl = ranova(rm);
    anova_tbl = anova(rm);
    
    mauchly_result = mauchly(rm);
    ranova_multcompare_results = multcompare(rm, 'Group');
    
    disp(' ');
    disp('MAUCHLY TEST FOR SPHERICITY (Determines which kind of p-value we are allowed to use from repeated-measures ANOVA)');
    disp(['MAUCHLY, p = ' num2str(table2array(mauchly_result(1, 4)))]);
    
    disp(' ');
    disp('REPEATED-MEASURES ANOVA (Post through Week 6 timepoints):');
    disp('Table showing effects of time and interaction effects:');
    disp(ranova_tbl);
    disp('Table showing group effects:');
    disp(anova_tbl);
    disp('Post-hoc multiple comparisons: ');
    disp(ranova_multcompare_results);
    
    %Now let's do a one-way ANOVA at each timepoint
    disp(' ');
    disp('ONE-WAY ANOVA AT EACH TIMEPOINT: ');
    
    dependent_variable = pre;
    independent_variable = {groups};
    [p, t2, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    [c, ~, ~, agroups] = multcompare(stats, 'display', 'off');
    sig_rows = find(c(:, 6) < 0.05);
    sig_rows_str = '';
    if (~isempty(sig_rows))
        for r = 1:length(sig_rows)
            this_row_idx = sig_rows(r);
            g1_str = char(agroups(c(this_row_idx, 1)));
            g2_str = char(agroups(c(this_row_idx, 2)));;
            g1_str = erase(g1_str, 'Group=');
            g2_str = erase(g2_str, 'Group=');
            sig_rows_str = [sig_rows_str ' (' g1_str ' vs ' g2_str ', p = ' num2str(c(this_row_idx, 6)) ') '];
        end
    end
    disp(['Pre: p = ' num2str(p) sig_rows_str]);
    
    dependent_variable = post;
    independent_variable = {groups};
    [p, t2, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    c = multcompare(stats, 'display', 'off');
    sig_rows = find(c(:, 6) < 0.05);
    sig_rows_str = '';
    if (~isempty(sig_rows))
        for r = 1:length(sig_rows)
            this_row_idx = sig_rows(r);
            g1_str = char(agroups(c(this_row_idx, 1)));
            g2_str = char(agroups(c(this_row_idx, 2)));;
            g1_str = erase(g1_str, 'Group=');
            g2_str = erase(g2_str, 'Group=');
            sig_rows_str = [sig_rows_str ' (' g1_str ' vs ' g2_str ', p = ' num2str(c(this_row_idx, 6)) ') '];
        end
    end
    disp(['Post: p = ' num2str(p) sig_rows_str]);
    
    dependent_variable = wk1;
    independent_variable = {groups};
    [p, t2, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    c = multcompare(stats, 'display', 'off');
    sig_rows = find(c(:, 6) < 0.05);
    sig_rows_str = '';
    if (~isempty(sig_rows))
        for r = 1:length(sig_rows)
            this_row_idx = sig_rows(r);
            g1_str = char(agroups(c(this_row_idx, 1)));
            g2_str = char(agroups(c(this_row_idx, 2)));;
            g1_str = erase(g1_str, 'Group=');
            g2_str = erase(g2_str, 'Group=');
            sig_rows_str = [sig_rows_str ' (' g1_str ' vs ' g2_str ', p = ' num2str(c(this_row_idx, 6)) ') '];
        end
    end
    disp(['Wk1: p = ' num2str(p) sig_rows_str]);
    
    dependent_variable = wk2;
    [p, t2, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    c = multcompare(stats, 'display', 'off');
    sig_rows = find(c(:, 6) < 0.05);
    sig_rows_str = '';
    if (~isempty(sig_rows))
        for r = 1:length(sig_rows)
            this_row_idx = sig_rows(r);
            g1_str = char(agroups(c(this_row_idx, 1)));
            g2_str = char(agroups(c(this_row_idx, 2)));;
            g1_str = erase(g1_str, 'Group=');
            g2_str = erase(g2_str, 'Group=');
            sig_rows_str = [sig_rows_str ' (' g1_str ' vs ' g2_str ', p = ' num2str(c(this_row_idx, 6)) ') '];
        end
    end
    disp(['Wk2: p = ' num2str(p) sig_rows_str]);
    
    dependent_variable = wk3;
    [p, t2, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    c = multcompare(stats, 'display', 'off');
    sig_rows = find(c(:, 6) < 0.05);
    sig_rows_str = '';
    if (~isempty(sig_rows))
        for r = 1:length(sig_rows)
            this_row_idx = sig_rows(r);
            g1_str = char(agroups(c(this_row_idx, 1)));
            g2_str = char(agroups(c(this_row_idx, 2)));;
            g1_str = erase(g1_str, 'Group=');
            g2_str = erase(g2_str, 'Group=');
            sig_rows_str = [sig_rows_str ' (' g1_str ' vs ' g2_str ', p = ' num2str(c(this_row_idx, 6)) ') '];
        end
    end
    disp(['Wk3: p = ' num2str(p) sig_rows_str]);
    
    dependent_variable = wk4;
    [p, t2, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    c = multcompare(stats, 'display', 'off');
    sig_rows = find(c(:, 6) < 0.05);
    sig_rows_str = '';
    if (~isempty(sig_rows))
        for r = 1:length(sig_rows)
            this_row_idx = sig_rows(r);
            g1_str = char(agroups(c(this_row_idx, 1)));
            g2_str = char(agroups(c(this_row_idx, 2)));;
            g1_str = erase(g1_str, 'Group=');
            g2_str = erase(g2_str, 'Group=');
            sig_rows_str = [sig_rows_str ' (' g1_str ' vs ' g2_str ', p = ' num2str(c(this_row_idx, 6)) ') '];
        end
    end
    disp(['Wk4: p = ' num2str(p) sig_rows_str]);
    
    dependent_variable = wk5;
    [p, t2, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    c = multcompare(stats, 'display', 'off');
    sig_rows = find(c(:, 6) < 0.05);
    sig_rows_str = '';
    if (~isempty(sig_rows))
        for r = 1:length(sig_rows)
            this_row_idx = sig_rows(r);
            g1_str = char(agroups(c(this_row_idx, 1)));
            g2_str = char(agroups(c(this_row_idx, 2)));;
            g1_str = erase(g1_str, 'Group=');
            g2_str = erase(g2_str, 'Group=');
            sig_rows_str = [sig_rows_str ' (' g1_str ' vs ' g2_str ', p = ' num2str(c(this_row_idx, 6)) ') '];
        end
    end
    disp(['Wk5: p = ' num2str(p) sig_rows_str]);
    
    dependent_variable = wk6;
    [p, t2, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    c = multcompare(stats, 'display', 'off');
    sig_rows = find(c(:, 6) < 0.05);
    sig_rows_str = '';
    if (~isempty(sig_rows))
        for r = 1:length(sig_rows)
            this_row_idx = sig_rows(r);
            g1_str = char(agroups(c(this_row_idx, 1)));
            g2_str = char(agroups(c(this_row_idx, 2)));;
            g1_str = erase(g1_str, 'Group=');
            g2_str = erase(g2_str, 'Group=');
            sig_rows_str = [sig_rows_str ' (' g1_str ' vs ' g2_str ', p = ' num2str(c(this_row_idx, 6)) ') '];
        end
    end
    disp(['Wk6: p = ' num2str(p) sig_rows_str]);
    
    %Now let's do a Kruskal-Wallis test at each timepoint
    disp(' ');
    disp('KRUSKAL-WALLIS TEST ACROSS ALL GROUPS AT EACH TIMEPOINT (NON-PARAMETRIC FORM OF ONE-WAY ANOVA): ');
    [p, tbl, stats] = kruskalwallis(pre, groups, 'off');
    disp(['Pre: p = ' num2str(p)]);
    [p, tbl, stats] = kruskalwallis(post, groups, 'off');
    disp(['Post: p = ' num2str(p)]);
    [p, tbl, stats] = kruskalwallis(wk1, groups, 'off');
    disp(['Wk1: p = ' num2str(p)]);
    [p, tbl, stats] = kruskalwallis(wk2, groups, 'off');
    disp(['Wk2: p = ' num2str(p)]);
    [p, tbl, stats] = kruskalwallis(wk3, groups, 'off');
    disp(['Wk3: p = ' num2str(p)]);
    [p, tbl, stats] = kruskalwallis(wk4, groups, 'off');
    disp(['Wk4: p = ' num2str(p)]);
    [p, tbl, stats] = kruskalwallis(wk5, groups, 'off');
    disp(['Wk5: p = ' num2str(p)]);
    [p, tbl, stats] = kruskalwallis(wk6, groups, 'off');
    disp(['Wk6: p = ' num2str(p)]);
    
    vns_1point6_mA_indices = find(strcmpi(groups, 'VNS 1.6 mA'));
    vns_point8_mA_indices = find(strcmpi(groups, 'VNS 0.8 mA'));
    vns_point4_mA_indices = find(strcmpi(groups, 'VNS 0.4 mA'));
    no_vns_indices = find(strcmpi(groups, 'No VNS'));
    
    %Now let's do a 2-sample KS test at each timepoint (0.8 mA vs Rehab)
    disp(' ');
    disp('TWO-SAMPLE KS-TEST OF 0.8 mA vs No VNS AT EACH TIMEPOINT (TESTS IF GROUPS BELONG TO SAME UNDERLYING DISTRIBUTION):');
    kstest_tail = 'smaller';
    kstest_ncomparisons = 3;
    kstest_alpha = 0.05 / kstest_ncomparisons;
    disp(['ks-test parameters: tail = ' kstest_tail ', num comparisons = ' num2str(kstest_ncomparisons) ', alpha = ' num2str(kstest_alpha)]);
    
    [h, p] = kstest2(wk1(vns_point8_mA_indices), wk1(no_vns_indices), 'tail', kstest_tail, 'alpha', kstest_alpha);
    disp(['Week 1: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = kstest2(wk2(vns_point8_mA_indices), wk2(no_vns_indices), 'tail', kstest_tail, 'alpha', kstest_alpha);
    disp(['Week 2: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = kstest2(wk3(vns_point8_mA_indices), wk3(no_vns_indices), 'tail', kstest_tail, 'alpha', kstest_alpha);
    disp(['Week 3: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = kstest2(wk4(vns_point8_mA_indices), wk4(no_vns_indices), 'tail', kstest_tail, 'alpha', kstest_alpha);
    disp(['Week 4: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = kstest2(wk5(vns_point8_mA_indices), wk5(no_vns_indices), 'tail', kstest_tail, 'alpha', kstest_alpha);
    disp(['Week 5: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = kstest2(wk6(vns_point8_mA_indices), wk6(no_vns_indices), 'tail', kstest_tail, 'alpha', kstest_alpha);
    disp(['Week 6: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    
    
    %Now let's do a Wilcoxon rank-sum test at each timepoint (0.8 mA vs Rehab)
    disp(' ');
    disp('WILCOXON RANK-SUM OF 0.8 mA vs No VNS AT EACH TIMEPOINT (NON-PARAMETRIC FORM OF T-TEST):');
    rstest_tail = 'right';
    rstest_ncomparisons = 3;
    rstest_alpha = 0.05 / rstest_ncomparisons;
    disp(['ranksum parameters: tail = ' rstest_tail ', num comparisons = ' num2str(rstest_ncomparisons) ', alpha = ' num2str(rstest_alpha)]);
    
    [p, h] = ranksum(wk1(vns_point8_mA_indices), wk1(no_vns_indices), 'tail', rstest_tail, 'alpha', rstest_alpha);
    disp(['Week 1: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [p, h] = ranksum(wk2(vns_point8_mA_indices), wk2(no_vns_indices), 'tail', rstest_tail, 'alpha', rstest_alpha);
    disp(['Week 2: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [p, h] = ranksum(wk3(vns_point8_mA_indices), wk3(no_vns_indices), 'tail', rstest_tail, 'alpha', rstest_alpha);
    disp(['Week 3: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [p, h] = ranksum(wk4(vns_point8_mA_indices), wk4(no_vns_indices), 'tail', rstest_tail, 'alpha', rstest_alpha);
    disp(['Week 4: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [p, h] = ranksum(wk5(vns_point8_mA_indices), wk5(no_vns_indices), 'tail', rstest_tail, 'alpha', rstest_alpha);
    disp(['Week 5: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [p, h] = ranksum(wk6(vns_point8_mA_indices), wk6(no_vns_indices), 'tail', rstest_tail, 'alpha', rstest_alpha);
    disp(['Week 6: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    
    %Now let's do a test for equal variance at each timepoint
    vartest_type = 'Bartlett';
    
    disp(' ');
    disp('F-TEST FOR EQUAL VARIANCE ACROSS ALL GROUPS AT EACH TIMEPOINT (USED TO DETERMINE EITHER WELCH T-TEST OR STUDENT T-TEST):');
    disp('If there is unequal variance across groups, a Welch t-test may be used instead of a Student t-test');
    disp(['vartestn parameters: type = ' vartest_type]);
    
    [p, stats] = vartestn(pre, groups, 'display', 'off', 'testtype', vartest_type);
    disp(['Pre: p = ' num2str(p)]);
    [p, stats] = vartestn(post, groups, 'display', 'off', 'testtype', vartest_type);
    disp(['Post: p = ' num2str(p)]);
    [p, stats] = vartestn(wk1, groups, 'display', 'off', 'testtype', vartest_type);
    disp(['Wk1: p = ' num2str(p)]);
    [p, stats] = vartestn(wk2, groups, 'display', 'off', 'testtype', vartest_type);
    disp(['Wk2: p = ' num2str(p)]);
    [p, stats] = vartestn(wk3, groups, 'display', 'off', 'testtype', vartest_type);
    disp(['Wk3: p = ' num2str(p)]);
    [p, stats] = vartestn(wk4, groups, 'display', 'off', 'testtype', vartest_type);
    disp(['Wk4: p = ' num2str(p)]);
    [p, stats] = vartestn(wk5, groups, 'display', 'off', 'testtype', vartest_type);
    disp(['Wk5: p = ' num2str(p)]);
    [p, stats] = vartestn(wk6, groups, 'display', 'off', 'testtype', vartest_type);
    disp(['Wk6: p = ' num2str(p)]);
    
    %Now let's do a 2-sample test for equal variance at each timepoint
    vt2_tail = 'both';
    vt2_ncomparisons = 3;
    vt2_alpha = 0.05 / vt2_ncomparisons;
    
    disp(' ');
    disp('F-TEST FOR EQUAL VARIANCE OF 0.8 mA VNS vs No VNS AT EACH TIMEPOINT (USED TO DETERMINE EITHER WELCH T-TEST OR STUDENT T-TEST):');
    disp('If there is unequal variance across groups, a Welch t-test may be used instead of a Student t-test');
    disp(['vartest2 parameters: tail = ' vt2_tail ', num comparisons = ' num2str(vt2_ncomparisons) ', alpha = ' num2str(vt2_alpha)]);
    
    [h, p] = vartest2(wk1(vns_point8_mA_indices), wk1(no_vns_indices), 'tail', vt2_tail, 'alpha', vt2_alpha);
    disp(['Week 1: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = vartest2(wk2(vns_point8_mA_indices), wk2(no_vns_indices), 'tail', vt2_tail, 'alpha', vt2_alpha);
    disp(['Week 2: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = vartest2(wk3(vns_point8_mA_indices), wk3(no_vns_indices), 'tail', vt2_tail, 'alpha', vt2_alpha);
    disp(['Week 3: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = vartest2(wk4(vns_point8_mA_indices), wk4(no_vns_indices), 'tail', vt2_tail, 'alpha', vt2_alpha);
    disp(['Week 4: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = vartest2(wk5(vns_point8_mA_indices), wk5(no_vns_indices), 'tail', vt2_tail, 'alpha', vt2_alpha);
    disp(['Week 5: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    [h, p] = vartest2(wk6(vns_point8_mA_indices), wk6(no_vns_indices), 'tail', vt2_tail, 'alpha', vt2_alpha);
    disp(['Week 6: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    
    %Now let's do t-tests at each time point (0.8 mA vs Rehab)
    disp(' ');
    disp('UNPAIRED T-TEST OF 0.8 mA vs No VNS AT EACH TIMEPOINT:');
    ttest_tail = 'right';
    ttest_ncomparisons = 3;
    ttest_alpha = 0.05 / ttest_ncomparisons;
    ttest_vartype = 'equal';
    disp(['t-test parameters: tail = ' ttest_tail ', variances = ' ttest_vartype ', num comparisons = ' num2str(ttest_ncomparisons) ', alpha = ' num2str(ttest_alpha)]);
    
    [h, p] = ttest2(wk1(vns_point8_mA_indices), wk1(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype);
    disp(['Week 1: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    if (p < 0.05)
        plot(3, 80, 'Marker', '*', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
    end
    
    [h, p] = ttest2(wk2(vns_point8_mA_indices), wk2(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype);
    disp(['Week 2: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    if (p < 0.05)
        plot(4, 80, 'Marker', '*', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
    end
    
    [h, p] = ttest2(wk3(vns_point8_mA_indices), wk3(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype);
    disp(['Week 3: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    if (p < 0.05)
        plot(5, 80, 'Marker', '*', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
    end
    
    [h, p] = ttest2(wk4(vns_point8_mA_indices), wk4(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype);
    disp(['Week 4: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    if (p < 0.05)
        plot(6, 80, 'Marker', '*', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
    end
    
    [h, p] = ttest2(wk5(vns_point8_mA_indices), wk5(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype);
    disp(['Week 5: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    if (p < 0.05)
        plot(7, 80, 'Marker', '*', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
    end

    [h, p] = ttest2(wk6(vns_point8_mA_indices), wk6(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype);
    disp(['Week 6: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    if (p < 0.05)
        plot(8, 80, 'Marker', '*', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2);
    end
    
    %Now let's do t-test of 0.8 mA vs all other groups at Wk6 timepoint
    disp(' ');
    disp('UNPAIRED T-TEST OF 0.8 mA vs all other groups AT WEEK 6:');
    ttest_tail = 'right';
    ttest_ncomparisons = 3;
    ttest_alpha = 0.05 / ttest_ncomparisons;
    ttest_vartype = 'equal';
    disp(['t-test parameters: tail = ' ttest_tail ', variances = ' ttest_vartype ', num comparisons = ' num2str(ttest_ncomparisons) ', alpha = ' num2str(ttest_alpha)]);
    
    [h, p] = ttest2(wk6(vns_point8_mA_indices), wk6(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype);
    disp(['Week 1: 0.8 mA vs No VNS, p = ' num2str(p)]);
    [h, p] = ttest2(wk6(vns_point8_mA_indices), wk6(vns_point4_mA_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype);
    disp(['Week 1: 0.8 mA VNS vs 0.4 mA VNS, p = ' num2str(p)]);
    [h, p] = ttest2(wk6(vns_point8_mA_indices), wk6(vns_1point6_mA_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype);
    disp(['Week 1: 0.8 mA VNS vs 1.6 mA VNS, p = ' num2str(p)]);
    
    %Now let's use the bootstrapping technique
    bootstrapping_npermutations = 1000;
    bootstrapping_samplesize = 'equal';
    
    disp(' ');
    disp('BOOTSTRAPPING TECHNIQUE FOLLOWED BY UNPAIRED T-TEST OF 0.8 mA VNS vs No VNS AT EACH TIMEPOINT:');
    disp('Warning: these p-values are non-deterministic and are subject to change with each run of the code!');
    disp(['Bootstrapping parameters: num permutations = ' num2str(bootstrapping_npermutations)]);
    disp('t-test parameters are identical to previous t-test parameters');
    
    p = bootstrap_ttest2(wk1(vns_point8_mA_indices), wk1(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype, 'bootstrap_permutations', bootstrapping_npermutations, 'bootstrap_samplesize', bootstrapping_samplesize);
    disp(['Week 1: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    p = bootstrap_ttest2(wk2(vns_point8_mA_indices), wk2(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype, 'bootstrap_permutations', bootstrapping_npermutations, 'bootstrap_samplesize', bootstrapping_samplesize);
    disp(['Week 2: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    p = bootstrap_ttest2(wk3(vns_point8_mA_indices), wk3(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype, 'bootstrap_permutations', bootstrapping_npermutations, 'bootstrap_samplesize', bootstrapping_samplesize);
    disp(['Week 3: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    p = bootstrap_ttest2(wk4(vns_point8_mA_indices), wk4(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype, 'bootstrap_permutations', bootstrapping_npermutations, 'bootstrap_samplesize', bootstrapping_samplesize);
    disp(['Week 4: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    p = bootstrap_ttest2(wk5(vns_point8_mA_indices), wk5(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype, 'bootstrap_permutations', bootstrapping_npermutations, 'bootstrap_samplesize', bootstrapping_samplesize);
    disp(['Week 5: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    p = bootstrap_ttest2(wk6(vns_point8_mA_indices), wk6(no_vns_indices), 'tail', ttest_tail, 'alpha', ttest_alpha, 'vartype', ttest_vartype, 'bootstrap_permutations', bootstrapping_npermutations, 'bootstrap_samplesize', bootstrapping_samplesize);
    disp(['Week 6: No VNS vs 0.8 mA VNS, p = ' num2str(p)]);
    
    %% Display the legend very last!
    
    legend(legend_pieces, legend_strings);
    
    %% Create another figure as a bar plot for wk 6 data only
    
    figure;
    hold on;
    
    wk6_1point6_mean = nanmean(wk6(vns_1point6_mA_indices));
    wk6_1point6_err = nanstd(wk6(vns_1point6_mA_indices)) / sqrt(length(vns_1point6_mA_indices));
    wk6_point8_mean = nanmean(wk6(vns_point8_mA_indices));
    wk6_point8_err = nanstd(wk6(vns_point8_mA_indices)) / sqrt(length(vns_point8_mA_indices));
    wk6_point4_mean = nanmean(wk6(vns_point4_mA_indices));
    wk6_point4_err = nanstd(wk6(vns_point4_mA_indices)) / sqrt(length(vns_point4_mA_indices));
    wk6_novns_mean = nanmean(wk6(no_vns_indices));
    wk6_novns_err = nanstd(wk6(no_vns_indices)) / sqrt(length(no_vns_indices));
    
    wk6_bar_graph = bar([wk6_novns_mean wk6_point4_mean wk6_point8_mean wk6_1point6_mean]);
    wk6_bar_graph.FaceColor = 'flat';
    wk6_bar_graph.CData(1, :) = colors(1, :);
    wk6_bar_graph.CData(2, :) = colors(3, :);
    wk6_bar_graph.CData(3, :) = colors(2, :);
    wk6_bar_graph.CData(4, :) = colors(4, :);
    set(gca, 'XTick', 1:4);
%     set(gca, 'XTickLabel', unique_groups);
    set(gca, 'XTickLabel', {});
    set(gca, 'XTickLabelRotation', 45);
    
    all_means = [wk6_novns_mean wk6_point4_mean wk6_point8_mean wk6_1point6_mean]
    all_errs = [wk6_novns_err wk6_point4_err wk6_point8_err wk6_1point6_err];
    
    errorbar([1 2 3 4], all_means, all_errs, all_errs, 'Marker', 'none', 'LineStyle', 'none', ...
        'MarkerFaceColor', 'k', 'Color', 'k', 'LineWidth', 2);
    
    
%     if (strcmpi(parameter_name, 'hit_rate'))
%         ylabel('Percent trials above 60 degrees');
%     elseif (strcmpi(parameter_name, 'maximal_turn_angle'))
%         ylabel('Maximal turn angle (degrees)');
%     end
%     
    xlim([0.25 4.75]);
    
    %% Create another figure to plot the "% recovery" of each group
    
    figure;
    hold on;
    
    percent_recovery_means = [];
    percent_recovery_errs = [];
    percent_recovery_anova_data = struct('group', {}, 'data', {});
    
    for g = 1:4
        this_group_idx = [];
        switch (g)
            case 1
                this_group_idx = no_vns_indices;
            case 2
                this_group_idx = vns_point4_mA_indices;
            case 3
                this_group_idx = vns_point8_mA_indices;
            case 4
                this_group_idx = vns_1point6_mA_indices;
        end
        
        if (~isempty(this_group_idx))
            pre_data = pre(this_group_idx);
            post_data = post(this_group_idx);
            wk6_data = wk6(this_group_idx);
            
            percent_recovery = ((wk6_data - post_data) ./ (pre_data - post_data)) * 100;
            for r = 1:length(percent_recovery)
                percent_recovery_anova_data(end+1) = struct('group', g, 'data', percent_recovery(r));
            end
            
            mean_percent_recovery = nanmean(percent_recovery);
            err_percent_recovery = nanstd(percent_recovery) / sqrt(length(percent_recovery));
            
            percent_recovery_means = [percent_recovery_means mean_percent_recovery];
            percent_recovery_errs = [percent_recovery_errs err_percent_recovery];
        end
    end
    
    google_colors = [66 133 244; 244 180 0; 219 68 55; 15 157 88];
    google_colors_mat = google_colors ./ 256;
    
    b = bar(percent_recovery_means);
    b.FaceColor = 'flat';
    b.CData(1, :) = google_colors_mat(1, :);
    b.CData(2, :) = google_colors_mat(2, :);
    b.CData(3, :) = google_colors_mat(3, :);
    b.CData(4, :) = google_colors_mat(4, :);
    
    set(gca, 'XTick', 1:4);
    set(gca, 'XTickLabel', {});
    
    errorbar([1 2 3 4], percent_recovery_means, [], percent_recovery_errs, 'Marker', 'none', 'LineStyle', 'none', ...
        'MarkerFaceColor', 'k', 'Color', 'k', 'LineWidth', 2);
    
    xlim([0.25 4.75]);
    
    
    dependent_variable = [percent_recovery_anova_data.data];
    independent_variable = {[percent_recovery_anova_data.group]};
    
    [p, tbl, stats] = anovan(dependent_variable, independent_variable, 'display', 'off', 'model', 'full', 'varnames', {'Group'});
    
    disp('Percent recovery table:');
    disp(tbl);
    
end















































