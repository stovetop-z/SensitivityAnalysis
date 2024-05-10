function [] = mean_x_dof(x_sens, t, params_char)
% Take the sensitivities and find the mean of each parameter at each DOF
band_start = 40;
band_end = 80;
means = [];

% ps = parameter sensitivities
for ps = 1:length(x_sens) 
    sens = x_sens{ps};

    % p = parameter, which corresponds to rows
    for p = 1:size(sens, 1)
        param_sens = [];
        
        % w = the frequency column, do this between the tremor band
        % where the tremor band is between 4-8 Hz
        for w = band_start:band_end
            % Now we need to access individual sensitivity vectors
            param_sens = [param_sens; sens(p, w)];
        end
        params = cell2mat(param_sens);
        means = [means; mean(params(:))];
    end
end

% Let's plot the data where each column of each row of the `means` will be the y
% and each row corresponds to a parameters
imagesc(means);

title(t);
end