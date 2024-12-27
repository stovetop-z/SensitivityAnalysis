clear; 
clear all;

addpath("functions/");
Parameters;
load("mats/SIE_Means.mat", "SIE_X_Means");
load("mats/SIE_X_Sensitivities.mat", "SIE_X_Sens");
load("mats/AIE_Means.mat", "AIE_X_Means");
load("mats/AIE_X_Sensitivities.mat", "AIE_X_Sens");


muscles = {'DELT1', 'DELT2', 'DELT3', 'PECM2', 'BIClong', 'BICshort', 'TRIlong', 'TRIlat', 'BRA', 'BRD', 'PT', 'FCR', 'FCU', 'ECR', 'ECU'};
ms = {1, muscles{1}; 2, muscles{2}; 3, muscles{3}; 4, muscles{4}; 5, muscles{5}; 6, muscles{6}; 7, muscles{7}; 8, muscles{8}; 9, muscles{9}; 10, muscles{10}; 11, muscles{11}; 12, muscles{12}; 13, muscles{13}; 14, muscles{14}; 15, muscles{15}};

%% Find the highest sensitivities
all_aie = [];

for m = 1:length(AIE_X_Sens)
    matrices = AIE_X_Sens{m};
    for pm = 1:size(matrices, 1)
        for wm = 1:size(matrices, 2)
            all_aie = [all_aie, matrices{pm, wm}];
        end
    end
end

all_aie = sort(all_aie, "ascend");
all_aie(length(all_aie) - 19: length(all_aie))

%% SIE and AIE
sens_x_plot(SIE_X_Sens, "Single Input Excitation Sensitivities at " + ms{1, 2}, params_char, ms{1, 1});
sens_x_plot(AIE_X_Sens, "All Input Excitation Sensitivities", params_char, 0);

%% Go into individual sensitivities
mean_x_plot(SIE_X_Means, "Mean of Single Input Excitation Sensitivities", params_char);
mean_x_plot(AIE_X_Means, "Means of All Input Excitation Sensitivities", params_char);

% Take the sensitivities and find the mean of each parameter at each DOF
band_start = 40;
band_end = 80;
means = [];

% ps = parameter sensitivities
for ps = 1:length(SIE_X_Sens) 
    sens = SIE_X_Sens{ps};

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
        means = [means; mean(params)];
    end
end

% Let's plot the data where each column of each row of the `means` will be the y
% and each row corresponds to a parameters
imagesc(means.');
set(gca, 'XTick', 1:numel(params_char), 'XTickLabel', params_char);
set(gca, 'YTick', 1:numel(muscles), 'YTickLabel', muscles);
title("Mean Total Output Sensitivities for Muscles"); 
colormap("turbo");
colorbar;

band_start = 40;
band_end = 80;
stds = [];

% ps = parameter sensitivities
for ps = 1:length(SIE_X_Sens) 
    sens = SIE_X_Sens{ps};

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
        stds = [stds; std(params)];
    end
end

% Let's plot the data where each column of each row of the `means` will be the y
% and each row corresponds to a parameters
imagesc(stds.');
set(gca, 'XTick', 1:numel(params_char), 'XTickLabel', params_char);
set(gca, 'YTick', 1:numel(muscles), 'YTickLabel', muscles);
title("Standard Deviation of Total Output Sensitivities for Muscles"); 
colormap("turbo");
colorbar;

%% Data normalized to their respective parameter matrices
load("mats/AIE_Org_NOT_NORM.mat");

figure;
colorbar;
tiledlayout(2, 4);
title("AIE Normalized to Parameter Matrices");
for m = 1:length(AIE_Org_N)
    switch(m)
        case 1
            nexttile;
            T1_img = heatmap(cell2mat(AIE_Org_N{m}));
            title("T1");
        case 2
            nexttile;
            T2_img = heatmap(cell2mat(AIE_Org_N{m}));
            title("T2");
        case 3
            nexttile;
            I_img = heatmap(AIE_Org_N{m});
            title("I");
        case 4
            nexttile;
            D_img = heatmap(AIE_Org_N{m});
            title("D");
        case 5
            nexttile;
            K_img = heatmap(AIE_Org_N{m});
            title("K");
        case 6
            nexttile;
            M_img = heatmap(AIE_Org_N{m});
            title("M");
        case 7
            nexttile;
            C_img = heatmap(AIE_Org_N{m});
            title("C");
        case 8
            nexttile;
            LL_img = heatmap(cell2mat(AIE_Org_N{m}));
            title("Link-Lengths");
    end
    c = colormap("turbo");
    colormap(c);
end

%% Data normalized to their respective parameter matrices
load("mats/AIE_Org_Not_Norm.mat");

figure;
colorbar;
tiledlayout(2, 4);
title("AIE Normalized to Parameter Matrices");
for m = 1:length(AIE_Org_N)
    switch(m)
        case 1
            nexttile;
            T1_img = heatmap(cell2mat(AIE_Org_N{m}));
            title("T1");
        case 2
            nexttile;
            T2_img = heatmap(cell2mat(AIE_Org_N{m}));
            title("T2");
        case 3
            nexttile;
            I_img = heatmap(AIE_Org_N{m});
            title("I");
        case 4
            nexttile;
            D_img = heatmap(AIE_Org_N{m});
            title("D");
        case 5
            nexttile;
            K_img = heatmap(AIE_Org_N{m});
            title("K");
        case 6
            nexttile;
            M_img = heatmap(AIE_Org_N{m});
            title("M");
        case 7
            nexttile;
            C_img = heatmap(AIE_Org_N{m});
            title("C");
        case 8
            nexttile;
            LL_img = heatmap(cell2mat(AIE_Org_N{m}));
            title("Link-Lengths");
    end
    c = colormap("turbo");
    colormap(c);
end

%% Data normalized to transfer function matrix
load("mats/AIE_Org_Not_Norm_OGPosture_1.mat");
figure;
colorbar;
tiledlayout(2, 4);
title("AIE Normalized");
for m = 1:length(AIE_Org_N)
    switch(m)
        case 1
            nexttile;
            T1_img = heatmap(cell2mat(AIE_Org_N{m}));
            title("T1");
        case 2
            nexttile;
            T2_img = heatmap(cell2mat(AIE_Org_N{m}));
            title("T2");
        case 3
            nexttile;
            I_img = heatmap(AIE_Org_N{m});
            title("I");
        case 4
            nexttile;
            D_img = heatmap(AIE_Org_N{m});
            title("D");
        case 5
            nexttile;
            K_img = heatmap(AIE_Org_N{m});
            title("K");
        case 6
            nexttile;
            M_img = heatmap(AIE_Org_N{m});
            title("M");
        case 7
            nexttile;
            C_img = heatmap(AIE_Org_N{m});
            title("C");
        case 8
            nexttile;
            LL_img = heatmap(cell2mat(AIE_Org_N{m}));
            title("Link-Lengths");
    end
    c = colormap("turbo");
    colormap(c);
end