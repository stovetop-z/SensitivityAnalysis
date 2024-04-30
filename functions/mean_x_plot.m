function [] = mean_x_plot(x_means, t, params_char)
sie_means = [];
for s = 1:length(x_means)
    figure;
    sens = cell2mat(x_means{s});

    % Now, we can create a heatmap
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
            pc = params_char(108:110);
    end

    % Plot heat map of means
    heatmap(pc, "Mean Sensitivity", sens);
    sie_means = [sie_means; sens.'];
end

% Create a table
sie_table = table(params_char.', sie_means, 'VariableNames', {'Parameters', 'MeanSensitivity'});

% Sort the table by 'MeanSensitivity' in ascending order
T = sortrows(sie_table, 'MeanSensitivity');

% Extract the sorted data
ParametersSorted = T.Parameters;
MeanSensitivitySorted = T.MeanSensitivity;

% Create a heatmap plot
heatmap(1:numel(ParametersSorted), 1, MeanSensitivitySorted', 'XData', ParametersSorted);

% Set labels and title
xlabel('Parameters');
ylabel('Mean Sensitivity');
title(t);

% Adjust the colorbar
colorbar;
end

