function [T] = mean_x_plot(AIE_Means, t, params_char, SortedParams)
means = [];

% Organize the parameters
for s = 1:length(AIE_Means)
    sens = abs(AIE_Means{s});
    switch(s)
        case 1
            pc = params_char(1);
        case 2
            pc = params_char(2);
        case 3
            pc = params_char(3:16);
        case 4
            pc = params_char(17:30);
        case 5
            pc = params_char(31:44);
        case 6
            pc = params_char(45:92);
        case 7
            pc = params_char(93:107);
        case 8
            pc = params_char(108:120);
    end

    means = [means; sens.'];
end

% Create a table
mean_table = table(params_char.', means, 'VariableNames', {'Parameters', 'MeanSensitivity'});

% Find the parameters less than 0.1
p = [];
for val = 1:size(means, 1)
    par = means(val);

    if par < 0.1
        p = [p, par];
    end
end
disp(length(p)/110);

if size(SortedParams) > 1
    % Match the order of `Table1` based on `Table2`
    [~, idx] = ismember(SortedParams.Parameters, mean_table.Parameters);
    
    % Sort Table1 rows to match the order in Table2
    T = mean_table(idx, :);
else
    % Sort the table by 'MeanSensitivity' in ascending order
    T = sortrows(mean_table, 'MeanSensitivity');
end

% Extract the sorted data
ParametersSorted = T.Parameters;
MeanSensitivitySorted = T.MeanSensitivity;

% Create a heatmap plot
imagesc(MeanSensitivitySorted.');
set(gca, 'XTick', 1:numel(ParametersSorted), 'XTickLabel', ParametersSorted);
set(gca, 'YTick', 1, 'YTickLabel', "Mean Sensitivity")

% Set labels and title
xlabel('Parameters');
title(t);

% Adjust the colorbar
colormap("turbo");
colorbar;
end

