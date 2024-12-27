function [] = mean_x_plot(x_means, t, params_char)
means = [];

for s = 1:length(x_means)
    figure;
    sens = abs(cell2mat(x_means{s}));

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

    heatmap(pc, "Mean Sensitivity", sens, 'Colormap', turbo);
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
disp (length(p)/110);

% Sort the table by 'MeanSensitivity' in ascending order
T = sortrows(mean_table, 'MeanSensitivity');

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

