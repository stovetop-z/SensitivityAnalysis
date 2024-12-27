function [] = sens_x_plot(sens, t, params_char, muscle)
figure;
frequency = (0.1:0.1:15) * 10;

% If we are using SIE sensitivity
if muscle > 0
    % Plotting algorithm for a graph of each parameter's flow through 0-15Hz
    for m = 1:length(sens)
        matrix_of_matrices = sens{m};
    
        for r = 1:size(matrix_of_matrices, 1)
            for c = 1:size(matrix_of_matrices, 2)
                param_matrix = matrix_of_matrices{r, c};
                params(c) = param_matrix(muscle);
            end
    
            plot(frequency, params); 
            hold on;
        end
    end
else % Using AIE sensitivity
    % Plotting algorithm for a graph of each parameter's flow through 0-15Hz
    for m = 1:length(sens)
        matrix_of_matrices = sens{m};
    
        for r = 1:size(matrix_of_matrices, 1)
            for c = 1:size(matrix_of_matrices, 2)
                param_matrix = matrix_of_matrices{r, c};
                params(c) = param_matrix;
            end
    
            plot(frequency, params); 
            hold on;
        end
    end
end

legend(params_char);
title(t);
xlabel("Frequency (Hz)");
ylabel("Sensitivities");
xregion(40, 80);
hold off;


end

