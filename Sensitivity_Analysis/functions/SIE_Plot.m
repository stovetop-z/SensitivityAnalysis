function [] = SIE_Plot(SIE_Sensitivities, dof_loc, muscle_loc, params_char, dof_to_muscles)
figure;
frequency = (0.1:0.1:15);

% Plotting algorithm for a graph of each parameter's flow through 0-15Hz
for m = 1:length(SIE_Sensitivities)
    matrix_of_matrices = SIE_Sensitivities{m};

    for p = 1:size(matrix_of_matrices, 1)
        for w = 1:size(matrix_of_matrices, 2)
            param_matrix = matrix_of_matrices{p, w};
            params(w) = param_matrix(dof_loc, muscle_loc);
        end

        plot(frequency, params);
        hold on;
    end
end

legend(params_char);
title("Single Input Excitation " + dof_to_muscles(dof_loc, muscle_loc));
xlabel("Frequency (Hz)");
ylabel("Sensitivities");
hold off;

end

