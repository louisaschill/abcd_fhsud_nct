function [F_stats, p_values, p_values_corr, eta_squared, partial_eta_squared, ...
    significant_indices, significant_indices_corrected] = ...
    extract_anova_results(all_results, var_name)
% Extract F-stats, p-values, eta squared, and partial eta squared for a specific variable or interaction term

% Initialize output variables
F_stats = [];
p_values = [];
p_values_corr = [];
eta_squared = [];
partial_eta_squared = [];
significant_indices = [];
significant_indices_corrected = [];

if ~iscell(all_results)
    results = all_results;

    % Check if the variable or interaction term exists in the results
    if any(strcmp(results.Properties.RowNames, var_name))
        % Apply Benjamini-Hochberg correction for multiple comparisons
        results.p_values_corr = mafdr(results.p_value, 'BHFDR', 'true');

        % Extract F-stat, p-values, eta squared, and partial eta squared
        F_stat = results.F_stats(var_name);
        p_value = results.p_value(var_name);
        p_value_corr = results.p_values_corr(var_name);
        eta_sq = results.eta_squared(var_name);
        partial_eta_sq = results.partial_eta_squared(var_name);

        % Store the results
        F_stats = [F_stats; F_stat];
        p_values = [p_values; p_value];
        p_values_corr = [p_values_corr; p_value_corr];
        eta_squared = [eta_squared; eta_sq];
        partial_eta_squared = [partial_eta_squared; partial_eta_sq];

        % Check if the original p-value is significant
        if p_value < 0.05
            significant_indices = 1;
        end

        % Check if the corrected p-value is significant
        if p_value_corr < 0.05
            significant_indices_corrected = 1;
        end
    end
else
    % Iterate over multiple results (e.g., across regions)
    for i = 1:numel(all_results)
        results = all_results{i};

        % Check if the variable or interaction term exists in the results
        if any(strcmp(results.Properties.RowNames, var_name))
            % Apply Benjamini-Hochberg correction
            results.p_values_corr = mafdr(results.p_value, 'BHFDR', 'true');

            % Extract F-stat, p-values, eta squared, and partial eta squared
            F_stat = results.F_stats(var_name);
            p_value = results.p_value(var_name);
            p_value_corr = results.p_values_corr(var_name);
            eta_sq = results.eta_squared(var_name);
            partial_eta_sq = results.partial_eta_squared(var_name);

            % Store the results
            F_stats = [F_stats; F_stat];
            p_values = [p_values; p_value];
            p_values_corr = [p_values_corr; p_value_corr];
            eta_squared = [eta_squared; eta_sq];
            partial_eta_squared = [partial_eta_squared; partial_eta_sq];

            % Check if the original p-value is significant
            if p_value < 0.05
                significant_indices = [significant_indices; i];
            end

            % Check if the corrected p-value is significant
            if p_value_corr < 0.05
                significant_indices_corrected = [significant_indices_corrected; i];
            end
        end
    end
end
end
