function [rat_exclusion_list, ordered_rat_exclusion_list] = KAMP_GetExclusionList ( kamp_project_data )

    exclusion_criterion = 11.0;
    vns_groups_list = {'VNS 0.4 mA', 'VNS 0.8 mA', 'VNS 1.6 mA'};
    rat_exclusion_list = {};
    ordered_rat_exclusion_list = zeros(size(kamp_project_data.rat_names));

    %Decide which rats to exclude
    for r = 1:length(kamp_project_data.rat_names)
        %First check to see if this rat is in a VNS group
        if (any(strcmpi(vns_groups_list, kamp_project_data.groups{r})))
            %If the rat is in a VNS group, then calculate the average
            %impedance for this rat during therapy
            this_rat_impedance_measurements = kamp_project_data.vns_cuff_impedance(r, 9:33);
            this_rat_impedance = nanmean(this_rat_impedance_measurements);
            
            %If the average impedance exceeds the exclusion criterion, add
            %this rat to the list of excluded rats
            if (this_rat_impedance >= exclusion_criterion)
                rat_exclusion_list = [rat_exclusion_list kamp_project_data.rat_names{r}];
                ordered_rat_exclusion_list(r) = 1;
            end    
        end
    end
    
    %Create the "ordered" rat exclusion list
    for r = 1:length(rat_exclusion_list)
        idx = find(strcmpi(kamp_project_data.rat_names, rat_exclusion_list{r}), 1, 'first');
        ordered_rat_exclusion_list(idx) = 1;
    end
    
end