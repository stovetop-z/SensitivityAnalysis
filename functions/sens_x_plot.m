function [] = sens_x_plot(sens, t, params_char, muscle)
figure;
frequency = (0.1:0.1:15);

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

            switch(m)
                case 1
                    plot(frequency, params, 'Color','r');
                case 2
                    plot(frequency, params, 'Color','#D95319');
                case 3
                    plot(frequency, params, 'Color','g');
                case 4
                    plot(frequency, params, 'Color','b');
                case 5
                    plot(frequency, params, 'Color','m');
                case 6
                    plot(frequency, params, 'Color','c');
                case 7
                    plot(frequency, params, 'Color','y');
                case 8
                    plot(frequency, params, 'Color','#7E2F8E');

            end
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
    
            switch(m)
                case 1
                    plot(frequency, params, 'Color','r');
                case 2
                    plot(frequency, params, 'Color','#D95319');
                case 3
                    plot(frequency, params, 'Color','g');
                case 4
                    plot(frequency, params, 'Color','b');
                case 5
                    plot(frequency, params, 'Color','m');
                case 6
                    plot(frequency, params, 'Color','c');
                case 7
                    plot(frequency, params, 'Color','y');
                case 8
                    plot(frequency, params, 'Color','#7E2F8E');

            end
            hold on;
        end
    end
end

title(t);
xlabel("Frequency (Hz)");
ylabel("Sensitivities");
reg = xregion(4, 8);
reg.DisplayName = "Tremor Band";
legend(params_char);
hold off;


end

