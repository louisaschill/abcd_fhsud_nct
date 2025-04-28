function [all_results, stats, overall_eta_squared, overall_R_squared] = run_ANOVA(response_variables, ...
    anova_vars, continuous_vars, var_names, interaction_term_sets, display)
% Run ANOVA using the specified variables for multiple response variables
% Initialize output variables
all_results = cell(size(response_variables, 2), 1);
overall_eta_squared = zeros(size(response_variables, 2), 1);
overall_R_squared = zeros(size(response_variables, 2), 1); % Initialize R^2 output
rowNames = var_names;

% Extend rowNames to accommodate interaction terms
if ~isempty(interaction_term_sets)
    for int = 1:size(interaction_term_sets, 2)
        int_terms = interaction_term_sets{int};
        rowNames{end+1} = strjoin(int_terms, ':');  % Handles 2 or 3-way terms
    end 
end

% Iterate over each response variable
for i = 1:size(response_variables, 2)
    response_variable = response_variables(:, i);

    % Construct the initial model matrix (without interaction terms)
    model_matrix = eye(length(var_names));

    % Add interaction terms for each set
    for j = 1:numel(interaction_term_sets)
        interaction_terms = interaction_term_sets{j};
        indices = nan(length(interaction_terms), 1);  % Now flexible for 2 or 3 terms
        
        % Find the indices of the interaction terms
        for k = 1:length(interaction_terms)
            indices(k) = find(contains(var_names, interaction_terms{k}));
        end
        
        if ~isempty(indices)
            % Add a row to the model matrix representing the interaction
            model_matrix = [model_matrix; zeros(1, length(var_names))];
            model_matrix(end, indices) = 1;
        end
    end

    % Run ANOVA for the current response variable
    if size(response_variables, 2) == 1
        [pvals, terms, stats, ~] = anovan(response_variable, anova_vars, ...
            'continuous', continuous_vars, 'varnames', var_names, ...
            'model', model_matrix, 'display', display);
    else
        [pvals, terms, stats(i), ~] = anovan(response_variable, anova_vars, ...
            'continuous', continuous_vars, 'varnames', var_names, ...
            'model', model_matrix, 'display', display);
    end

    % Adjust p-values using Benjamini-Hochberg FDR
    p_corr = mafdr(pvals, 'bhfdr', 'true');

    % Extract F-statistics and p-values
    F_stats = cell2mat(terms(2:end-2, 6));

    % Calculate eta squared and partial eta squared for terms
    SS = cell2mat(terms(2:end-2, 2)); % Sum of squares for terms
    SS_total = terms{end, 2}; % Total sum of squares
    SS_error = terms{end-1, 2}; % Error sum of squares

    eta_squared = SS / SS_total; % Eta squared
    partial_eta_squared = SS ./ (SS + SS_error); % Partial eta squared

    % Calculate overall model eta squared and R-squared
    SS_model = sum(SS); % Sum of squares for all model terms
    overall_eta_squared(i) = SS_model / SS_total; % Overall eta squared
    overall_R_squared(i) = 1 - (SS_error / SS_total); % Overall R-squared

    % Create a table of results
    table_results = array2table([F_stats, pvals, p_corr, eta_squared, partial_eta_squared], ...
        'VariableNames', {'F_stats', 'p_value', 'p_value_corr', 'eta_squared', 'partial_eta_squared'}, ...
        'RowNames', rowNames);

    % Store the results for this response variable
    if size(response_variables, 2) == 1
        all_results = table_results;
    else
        all_results{i} = table_results;
    end
end
end
