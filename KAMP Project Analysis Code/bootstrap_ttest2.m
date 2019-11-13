function p_result = bootstrap_ttest2 (s1, s2, varargin)

    %Parse input parameters
    defaultTail = 'both';
    defaultAlpha = 0.05;
    defaultVartype = 'equal';
    defaultBootstrapSampleSize = 'equal';
    defaultBootstrapPermutations = 100;
    p = inputParser;
    p.CaseSensitive = false;
    addOptional(p, 'tail', defaultTail);
    addOptional(p, 'alpha', defaultAlpha);
    addOptional(p, 'vartype', defaultVartype);
    addOptional(p, 'bootstrap_samplesize', defaultBootstrapSampleSize);
    addOptional(p, 'bootstrap_permutations', defaultBootstrapPermutations);
    parse(p, varargin{:});
    tail = p.Results.tail;
    alpha = p.Results.alpha;
    vartype = p.Results.vartype;
    bootstrap_samplesize = p.Results.bootstrap_samplesize;
    bootstrap_permutations = p.Results.bootstrap_permutations;
    
    %Determine new sample size for each group
    numeric_samplesize_s1 = length(s1);
    numeric_samplesize_s2 = length(s2);
    if (isnumeric(bootstrap_samplesize))
        numeric_samplesize_s1 = bootstrap_samplesize;
        numeric_samplesize_s2 = bootstrap_samplesize;
    elseif (ischar(bootstrap_samplesize))
        if (strcmpi(bootstrap_samplesize, 'max'))
            s = max([length(s1) length(s2)]);
            numeric_samplesize_s1 = s;
            numeric_samplesize_s2 = s;
        elseif (strcmpi(bootstrap_samplesize, 'min'))
            s = min([length(s1) length(s2)]);
            numeric_samplesize_s1 = s;
            numeric_samplesize_s2 = s;
        end
    end
    
    p_statistics = nan(1, bootstrap_permutations);
    for i = 1:bootstrap_permutations
        %Draw a random sample from each original dataset (with replacement)
        s1_indices = randsample(length(s1), numeric_samplesize_s1, true);
        s2_indices = randsample(length(s2), numeric_samplesize_s2, true);
        s1_new = s1(s1_indices);
        s2_new = s2(s2_indices);

        [~, p] = ttest2(s1_new, s2_new, 'tail', tail, 'alpha', alpha, 'vartype', vartype);
        p_statistics(i) = p;
    end
    
    p_result = nanmean(p_statistics);
end