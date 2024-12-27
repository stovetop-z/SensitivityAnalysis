function [] = plot_imagesc(T, t)
% Extract the sorted data
ParametersSorted = T.Parameters;
MeanSensitivitySorted = T.MeanSensitivity;

% Create a heatmap plot
imagesc(MeanSensitivitySorted.');
set(gca, 'XTick', 1:numel(ParametersSorted), 'XTickLabel', ParametersSorted);
set(gca, 'YTick', 1, 'YTickLabel', '');

% Set labels and title
xlabel('Parameters');
title(t);

% Adjust the colorbar
colormap("turbo");
colorbar;
end

