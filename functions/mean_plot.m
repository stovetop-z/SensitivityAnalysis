function [] = mean_plot(Means, dof_loc, params_char, t)

all_means = [];
for m = 1:size(Means, 2)
    figure;
    sens = Means{m};
    means = [];

    for p = 1:size(sens, 2)
        p_means = sens{p};
        param = abs(p_means(dof_loc));

        means = [means, param];
    end
    switch(m)
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
    end

        % Plot heat map of means
    heatmap(pc, "Mean Sensitivity", means, 'Colormap', turbo);
    all_means = [all_means; means.'];
end

% Create a table
mean_table = table(params_char(1:107).', all_means, 'VariableNames', {'Parameters', 'MeanSensitivity'});

% Sort the table by 'MeanSensitivity' in ascending order
T = sortrows(mean_table, 'MeanSensitivity');

% Extract the sorted data
ParametersSorted = T.Parameters;
MeanSensitivitySorted = T.MeanSensitivity;

% Create a heatmap plot
heatmap(1:numel(ParametersSorted), 1, MeanSensitivitySorted', 'XData', ParametersSorted, 'Colormap', turbo);

% Set labels and title
xlabel('Parameters');
ylabel('Mean Sensitivity');
title(t);

% Adjust the colorbar
colorbar;

end