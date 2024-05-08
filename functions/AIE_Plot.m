function [] = AIE_Plot(AIE_Sens, dof_loc, params_char, dof)

frequency = 0.1:0.1:15;
for m = 1:size(AIE_Sens, 2)
    sens = AIE_Sens{m};
    
    for r = 1:size(sens, 1)
        params = [];
        for c = 1:size(sens, 2)
            par = sens{r, c};
            params = [params, par(dof_loc)];
        end

        plot(frequency, params);
        hold on;
    end
end

legend(params_char);
title("All Input Excitation " + dof);
xlabel("Frequency (Hz)");
ylabel("Sensitivities");
hold off;
end

