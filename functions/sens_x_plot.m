function [] = sens_x_plot(sens, t, params_char)
figure;
frequency = (0.1:0.1:15);


% Plotting algorithm for a graph of each parameter's flow through 0-15Hz
for m = 1:length(sens)
    matrix_of_matrices = sens{m};

    for r = 1:size(matrix_of_matrices, 1)
        for c = 1:size(matrix_of_matrices, 2)
            param_matrix{c} = matrix_of_matrices{r, c};
            params(c) = param_matrix{c}(1);
        end

        plot(frequency, params); 
        hold on;
    end
end

legend(params_char);
title(t);
xlabel("Frequency (Hz)");
ylabel("Sensitivities");
hold off;
end

