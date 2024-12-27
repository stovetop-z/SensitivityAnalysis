clear; 
clear all;

Parameters;
%% Load Sensitivities
load("mats/AIE_Postures_Means.mat", "AIE_Postures_Means");
load("mats/AIE_Postures_Sens.mat", "AIE_Postures_Sens");
load("mats/AIE_Postures_Params_Matrix.mat", "AIE_Postures_Params_Matrix");
load("mats/AIE_Postures_Stds.mat", "AIE_Postures_Stds");

muscles = {'DELT1', 'DELT2', 'DELT3', 'PECM2', 'BIClong', 'BICshort', 'TRIlong', 'TRIlat', 'BRA', 'BRD', 'PT', 'FCR', 'FCU', 'ECR', 'ECU'};
ms = {1, muscles{1}; 2, muscles{2}; 3, muscles{3}; 4, muscles{4}; 5, muscles{5}; 6, muscles{6}; 7, muscles{7}; 8, muscles{8}; 9, muscles{9}; 10, muscles{10}; 11, muscles{11}; 12, muscles{12}; 13, muscles{13}; 14, muscles{14}; 15, muscles{15}};

%% Means
poses = ["Posture_1" "Posture_2" "Posture_3" "Posture_4" "Posture_5" "Posture_6" "Posture_7"];
postures = [1 2 3 4 5 6 7];
% This will hold our mean and std tables
tables = {};
t = {"Mean of All Postures of the Means of AIE Sensitivities", "Mean Standard Deviations of All Postures"};

tables{1} = mean_x_plot(AIE_Postures_Means, "Mean of All Postures of the Means of AIE Sensitivities", params_char, table({1}));
T = table({1});
tables{2} = mean_x_plot(AIE_Postures_Stds, "Mean Standard Deviations of All Postures", params_char, T);
for p = 1:length(postures)
    posture = postures(p);

    % Get posture mat location
    aie_p_pos_loc = "mats/AIE_Means_Inputs_" + posture + ".mat";
    aie_pos_loc = "mats/AIE_Sens_Inputs_" + posture + ".mat";
    aie_m_pos_loc = "mats/AIE_Params_Matrix_Inputs_" + posture + ".mat";

    % Load postures sensitivities
    AIE_Pos_Means = load(aie_p_pos_loc);
    AIE_Pos_Sens = load(aie_pos_loc);
    AIE_Pos_Params_Matrix = load(aie_m_pos_loc);
    aie_mean_field = fieldnames(AIE_Pos_Means);
    aie_matrix_field = fieldnames(AIE_Pos_Params_Matrix);
    aie_field = fieldnames(AIE_Pos_Sens);
    Pos_AIE_Means = AIE_Pos_Means.(aie_mean_field{1});
    Pos_AIE_Params_Matrix = AIE_Pos_Params_Matrix.(aie_matrix_field{1});
    Pos_AIE_Sens = AIE_Pos_Sens.(aie_field{1});
    
    t{2 + p} = "Means of AIE Sensitivities at Posture " + num2str(postures(p));
    tables{p + 2} = mean_x_plot(Pos_AIE_Means, t, params_char, T);
end

% Create figure
fig = figure('Name', 'Postural AIE Senstivities', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);

% Create left panel for 7 tables
leftPanel = uipanel(fig, 'Title', 'Postures', 'FontSize', 10, 'Position', [0.05, 0.05, 0.4, 0.9], 'TitlePosition', 'centertop');

% Create a tiled layout inside the left panel (7 tiles, 1 column)
leftLayout = tiledlayout(leftPanel, 7, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Create right panel for 2 tables
rightPanel = uipanel(fig, 'Title', 'Mean of All Postures and Standard Deviation', 'FontSize', 10, 'Position', [0.55, 0.05, 0.4, 0.9], 'TitlePosition', 'centertop');

% Create a tiled layout inside the right panel (2 tiles, 1 column)
rightLayout = tiledlayout(rightPanel, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot 7 tables in the left panel using tiledlayout
for i = 1:7
    nexttile(leftLayout, i);  % Use nexttile for the left panel
    plot_imagesc(tables{2 + i}, t{2 + i});
end

% Plot 2 tables in the right panel using tiledlayout
for i = 1:2
    nexttile(rightLayout, i);  % Use nexttile for the right panel
    plot_imagesc(tables{i}, t{i});
end

%% AIE frequency plot
f2 = figure;
frequency = (0.1:0.1:15);
num = 1;
% Plotting algorithm for a graph of each parameter's flow through 0-15Hz
for m = 1:length(AIE_Postures_Sens)
    
    if m == 1 || m == 2
        params = AIE_Postures_Sens{m};
        plot(frequency, params, 'Color','r');
        m = 2;
        params = AIE_Postures_Sens{m};
        plot(frequency, params, 'Color','#D95319');
        m = 3;
    end

    matrix_of_matrices = AIE_Postures_Sens{m};

    for c = 1:size(matrix_of_matrices, 2)
        params = squeeze(matrix_of_matrices(:, c));

        switch(m)
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
        num = num + 1;
    end
    disp(num);
    hold on;
end
xlabel("Frequency (Hz)");
ylabel("Sensitivities");
reg = xregion(4, 8);
reg.DisplayName = "Tremor Band";
legend(params_char);
hold off;

%% Parameters in their respective matrices to show relative sensitivity to other parameters within same category (i.e. I matrix parameters)
load("mats/AIE_Postures_Params_Matrix.mat");

t1_P_Matrix = (AIE_Postures_Params_Matrix{1}) ./ t1(1);
t2_P_Matrix = (AIE_Postures_Params_Matrix{2}) ./ t2(1);
I_Sens = AIE_Postures_Params_Matrix{3};
D_Sens = AIE_Postures_Params_Matrix{4};
K_Sens = AIE_Postures_Params_Matrix{5};
M_Sens = AIE_Postures_Params_Matrix{6};
C_Sens = AIE_Postures_Params_Matrix{7};
L_P_Matrix = abs(AIE_Postures_Params_Matrix{8} ./ L);

I_Params =   [ I_Sens(1),  I_Sens(2),  I_Sens(3),  I_Sens(4),  I_Sens(5),  I_Sens(6),  I_Sens(7);...
               I_Sens(2),  I_Sens(8),  I_Sens(9), I_Sens(10), I_Sens(11), I_Sens(12), I_Sens(13);...
               I_Sens(3),  I_Sens(9), I_Sens(14), I_Sens(15), I_Sens(16), I_Sens(17), I_Sens(18);...
               I_Sens(4), I_Sens(10), I_Sens(15), I_Sens(19),          0, I_Sens(20), I_Sens(21);...
               I_Sens(5), I_Sens(11), I_Sens(16),          0, I_Sens(22),          0,          0;...
               I_Sens(6), I_Sens(12), I_Sens(17), I_Sens(20),          0, I_Sens(23),          0;...
               I_Sens(7), I_Sens(13), I_Sens(18), I_Sens(21),          0,          0, I_Sens(24)];
I_P_Matrix = {abs(I_Params./I{1}), abs(I_Params./I{2}), abs(I_Params./I{3}), abs(I_Params./I{4}), abs(I_Params./I{5}), abs(I_Params./I{6}), abs(I_Params./I{7})};

D_Params =    [ D_Sens(1), D_Sens(2), D_Sens(3), D_Sens(4),     0,      0,      0;... 
                D_Sens(1), D_Sens(5), D_Sens(6), 0,             0,      0,      0;...
                D_Sens(2), D_Sens(6), D_Sens(7), 0,             0,      0,      0;...
                D_Sens(3), 0,         0,         D_Sens(8),     0,      0,      0;...
                0,         0,         0,         0,         D_Sens(9),  D_Sens(10), D_Sens(11);...
                0,         0,         0,         0,         D_Sens(10), D_Sens(12), D_Sens(13);...
                0,         0,         0,         0,         D_Sens(11), D_Sens(13), D_Sens(14) ];
D_P_Matrix = abs(D_Params./D);

K_Params =    [ K_Sens(1), K_Sens(2), K_Sens(3), K_Sens(4),     0,      0,      0;...
                K_Sens(1), K_Sens(5), K_Sens(6), 0,             0,      0,      0;...
                K_Sens(2), K_Sens(6), K_Sens(7), 0,             0,      0,      0;...
                K_Sens(3), 0,         0,         K_Sens(8),     0,      0,      0;...
                0,         0,         0,         0,         K_Sens(9),  K_Sens(10), K_Sens(11);...
                0,         0,         0,         0,         K_Sens(10), K_Sens(12), K_Sens(13);...
                0,         0,         0,         0,         K_Sens(11), K_Sens(13), K_Sens(14) ];
K_P_Matrix = abs(K_Params./K);

M_Params =    1e-3*[M_Sens(1), M_Sens(2), M_Sens(3), M_Sens(4), M_Sens(5), M_Sens(6), M_Sens(7),     0,    0,    0,    0,     0,    0,     0,     0;...
                    M_Sens(8),     M_Sens(9),      M_Sens(10),      M_Sens(11),     M_Sens(12),     M_Sens(13),      M_Sens(14),     0,    0,    0,    0,     0,    0,     0,     0;...
                    M_Sens(15),      M_Sens(16),      M_Sens(17),      M_Sens(18),      M_Sens(19),     M_Sens(20),     M_Sens(21),     0,    0,    0,    0,     0,    0,     0,     0;...
                       0,         0,          0,         0,      M_Sens(22),     M_Sens(23),     M_Sens(24), M_Sens(25), M_Sens(26), M_Sens(27), M_Sens(28),  M_Sens(29), M_Sens(30),  M_Sens(31), M_Sens(32);...
                       0,         0,          0,         0,     M_Sens(33),    M_Sens(34),         0,     0,    0, M_Sens(35), M_Sens(36),  M_Sens(37), M_Sens(38), M_Sens(39), M_Sens(40);...
                       0,         0,          0,         0,         0,        0,         0,     0,    0,    0,    0,  M_Sens(41), M_Sens(42), M_Sens(43), M_Sens(44);...
                       0,         0,          0,         0,         0,        0,         0,     0,    0,    0,    0, M_Sens(45), M_Sens(46), M_Sens(47),  M_Sens(48)];
M_P_Matrix = {abs(M_Params./M{1}), abs(M_Params./M{2}), abs(M_Params./M{3}), abs(M_Params./M{4}), abs(M_Params./M{5}), abs(M_Params./M{6}), abs(M_Params./M{7})};

C_Params = diag([C_Sens(1), C_Sens(2), C_Sens(3), C_Sens(4), C_Sens(5), C_Sens(6), C_Sens(7), C_Sens(8), C_Sens(9), C_Sens(10), C_Sens(11), C_Sens(12), C_Sens(13), C_Sens(14), C_Sens(15)]);
C_P_Matrix = abs(C_Params./C);

f3 = figure;
t = tiledlayout(5, 4);
title(t, 'Mean Sensitivities Non-Normalized and Compared Within Parameter Matrices');
nexttile(t);
heatmap(t1_P_Matrix, 'Title', 'T1', 'Colormap', turbo);
nexttile(t);
heatmap(t2_P_Matrix, 'Title', 'T2', 'Colormap', turbo);
nexttile(t);

I_P_Matrix{1}(isnan(I_P_Matrix{1})) = 0;
heatmap(I_P_Matrix{1}, 'Title', 'I Posture 1', 'Colormap', turbo);
nexttile(t);
I_P_Matrix{2}(isnan(I_P_Matrix{2})) = 0;
heatmap(I_P_Matrix{2}, 'Title', 'I Posture 2', 'Colormap', turbo);
nexttile(t);
I_P_Matrix{3}(isnan(I_P_Matrix{3})) = 0;
heatmap(I_P_Matrix{3}, 'Title', 'I Posture 3', 'Colormap', turbo);
nexttile(t);
I_P_Matrix{4}(isnan(I_P_Matrix{4})) = 0;
heatmap(I_P_Matrix{4}, 'Title', 'I Posture 4', 'Colormap', turbo);
nexttile(t);
I_P_Matrix{4}(isnan(I_P_Matrix{4})) = 0;
heatmap(I_P_Matrix{5}, 'Title', 'I Posture 5', 'Colormap', turbo);
nexttile(t);
I_P_Matrix{6}(isnan(I_P_Matrix{6})) = 0;
heatmap(I_P_Matrix{6}, 'Title', 'I Posture 6', 'Colormap', turbo);
nexttile(t);
I_P_Matrix{7}(isnan(I_P_Matrix{7})) = 0;
heatmap(I_P_Matrix{7}, 'Title', 'I Posture 7', 'Colormap', turbo);
nexttile(t);

D_P_Matrix(isnan(D_P_Matrix)) = 0;
heatmap(D_P_Matrix, 'Title', 'D', 'Colormap', turbo);
nexttile(t);

K_P_Matrix(isnan(K_P_Matrix)) = 0;
heatmap(K_P_Matrix, 'Title', 'K', 'Colormap', turbo);
nexttile(t);

M_P_Matrix{1}(isnan(M_P_Matrix{1})) = 0;
heatmap(M_P_Matrix{1}, 'Title', 'M Posture 1', 'Colormap', turbo);
nexttile(t);
M_P_Matrix{2}(isnan(M_P_Matrix{2})) = 0;
heatmap(M_P_Matrix{2}, 'Title', 'M Posture 2', 'Colormap', turbo);
nexttile(t);
M_P_Matrix{3}(isnan(M_P_Matrix{3})) = 0;
heatmap(M_P_Matrix{3}, 'Title', 'M Posture 3', 'Colormap', turbo);
nexttile(t);
M_P_Matrix{4}(isnan(M_P_Matrix{4})) = 0;
heatmap(M_P_Matrix{4}, 'Title', 'M Posture 4', 'Colormap', turbo);
nexttile(t);
M_P_Matrix{5}(isnan(M_P_Matrix{5})) = 0;
heatmap(M_P_Matrix{5}, 'Title', 'M Posture 5', 'Colormap', turbo);
nexttile(t);
M_P_Matrix{6}(isnan(M_P_Matrix{6})) = 0;
heatmap(M_P_Matrix{6}, 'Title', 'M Posture 6', 'Colormap', turbo);
nexttile(t);
M_P_Matrix{7}(isnan(M_P_Matrix{7})) = 0;
heatmap(M_P_Matrix{7}, 'Title', 'M Posture 7', 'Colormap', turbo);
nexttile(t);

C_P_Matrix(isnan(C_P_Matrix)) = 0;
heatmap(C_P_Matrix, 'Title', 'C', 'Colormap', turbo);
nexttile(t);

heatmap(L_P_Matrix, 'Title', 'Link-Lengths', 'Colormap', turbo);