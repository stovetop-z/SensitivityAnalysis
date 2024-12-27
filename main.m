%% Sensitivity Analysis %%
% Steps to run:
% 1. Create robotic arm to compute jacobian (if applicable)
% 2. Change specs of analysis
%  - 'postures', there are 7 total postures to analyze
% 3. Run Analysis
% 4. Plot data
clear;
clear all;


%% Amplitude and theta inputs
A = rand(100, 15, 1, 1);
theta = rand(100, 15, 1, 1) * 2 * pi;

% Create matrix of inputs pre calculated to save computation time later
INPUT = A .* exp(1).^(-1i.*theta);

%% Tremor band
frequency_band = (0.1:0.1:15);
w2 = 2 * pi * frequency_band;

%% Step 1 - Compute Jacobian - SKIP IF DONE
% create_robot_arm;

%% Step 2 - Change specs of analysis
poses = ["Posture_1" "Posture_2" "Posture_3" "Posture_4" "Posture_5" "Posture_6" "Posture_7"];
postures = [1 2 3 4 5 6 7];

%% Step 3 - Run Analysis

% Load parameters and other necessary data
disp("Loading parameters...")
Parameters;

%% Loop through analysis with given postures
disp("Welcome! Running All Input Excitation Sensitivity Analysis...")
for p = 1:length(postures)    
    disp(poses(postures(p)));
    pos = postures(p);

    tic
    AIE_Sensitivity_Analysis(I{pos}, D, K, C, M{pos}, t1, t2, L, J_sp{pos}, J_p{pos}, INPUT, w2, poses(postures(p)));
    toc
end

%% Analysis continued : computing mean of all postures
disp("Computing mean of all postures...");
tic;

AIE_I_Means = zeros(length(postures), 14);
AIE_D_Means = zeros(length(postures), 14);
AIE_K_Means = zeros(length(postures), 14);
AIE_M_Means = zeros(length(postures), 48);
AIE_C_Means = zeros(length(postures), 15);
AIE_t1_Means = zeros(length(postures), 1);
AIE_t2_Means = zeros(length(postures), 1);
AIE_L_Means = zeros(length(postures), 3);

for p = 1:length(postures)
    posture = poses(postures(p));

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
    
    for param = 1:8
        AIE_P_Means = Pos_AIE_Means{param};
        AIE_P_Sens = Pos_AIE_Sens{param};
        AIE_P_Matrix = Pos_AIE_Params_Matrix{param};
        if param == 1
            AIE_t1_Means(p, :) = AIE_P_Means;
            AIE_t1_Sens(p, :) = AIE_P_Sens;
            AIE_t1_Matrix(p, :) = AIE_P_Matrix;
        end
        if param == 2
            AIE_t2_Means(p, :) = AIE_P_Means;
            AIE_t2_Sens(p, :) = AIE_P_Sens;
            AIE_t2_Matrix(p, :) = AIE_P_Matrix;
        end
        if param == 3
            AIE_I_Means(p, :) = AIE_P_Means;
            AIE_I_Sens(p, :, :) = AIE_P_Sens;
            AIE_I_Matrix(p, :, :) = AIE_P_Matrix;
        end
        if param == 4
            AIE_D_Means(p, :) = AIE_P_Means;
            AIE_D_Sens(p, :, :) = AIE_P_Sens;
            AIE_D_Matrix(p, :, :) = AIE_P_Matrix;
        end
        if param == 5
            AIE_K_Means(p, :) = AIE_P_Means;
            AIE_K_Sens(p, :, :) = AIE_P_Sens;
            AIE_K_Matrix(p, :, :) = AIE_P_Matrix;
        end
        if param == 6
            AIE_M_Means(p, :) = AIE_P_Means;
            AIE_M_Sens(p, :, :) = AIE_P_Sens;
            AIE_M_Matrix(p, :, :) = AIE_P_Matrix;
        end
        if param == 7
            AIE_C_Means(p, :) = AIE_P_Means;
            AIE_C_Sens(p, :, :) = AIE_P_Sens;
            AIE_C_Matrix(p, :, :) = AIE_P_Matrix;
        end
        if param == 8
            AIE_L_Means(p, :) = AIE_P_Means;
            AIE_L_Sens(p, :, :) = AIE_P_Matrix;
            AIE_L_Matrix(p, :, :) = AIE_P_Matrix;
        end
    end
end

% Take the mean of all the above postures
AIE_I_Std = squeeze(std(AIE_I_Means));
AIE_I_Means = squeeze(mean(AIE_I_Means));
AIE_D_Std = squeeze(std(AIE_D_Means));
AIE_D_Means = squeeze(mean(AIE_D_Means));
AIE_K_Std = squeeze(std(AIE_K_Means));
AIE_K_Means = squeeze(mean(AIE_K_Means));
AIE_M_Std = squeeze(std(AIE_M_Means));
AIE_M_Means = squeeze(mean(AIE_M_Means));
AIE_C_Std = squeeze(std(AIE_C_Means));
AIE_C_Means = squeeze(mean(AIE_C_Means));
AIE_t1_Std = squeeze(std(AIE_t1_Means));
AIE_t1_Means = squeeze(mean(AIE_t1_Means));
AIE_t2_Std = squeeze(std(AIE_t2_Means));
AIE_t2_Means = squeeze(mean(AIE_t2_Means));
AIE_L_Std = squeeze(std(AIE_L_Means));
AIE_L_Means = squeeze(mean(AIE_L_Means));

AIE_I_Sens = squeeze(mean(AIE_I_Sens));
AIE_D_Sens = squeeze(mean(AIE_D_Sens));
AIE_K_Sens = squeeze(mean(AIE_K_Sens));
AIE_M_Sens = squeeze(mean(AIE_M_Sens));
AIE_C_Sens = squeeze(mean(AIE_C_Sens));
AIE_t1_Sens = squeeze(mean(AIE_t1_Sens));
AIE_t2_Sens = squeeze(mean(AIE_t2_Sens));
AIE_L_Sens = squeeze(mean(AIE_L_Sens));

AIE_I_Matrix = squeeze(mean(AIE_I_Matrix));
AIE_D_Matrix = squeeze(mean(AIE_D_Matrix));
AIE_K_Matrix = squeeze(mean(AIE_K_Matrix));
AIE_M_Matrix = squeeze(mean(AIE_M_Matrix));
AIE_C_Matrix = squeeze(mean(AIE_C_Matrix));
AIE_t1_Matrix = squeeze(mean(AIE_t1_Matrix));
AIE_t2_Matrix = squeeze(mean(AIE_t2_Matrix));
AIE_L_Matrix = squeeze(mean(AIE_L_Matrix));
toc;

AIE_Postures_Stds = {AIE_t1_Std, AIE_t2_Std, AIE_I_Std, AIE_D_Std, AIE_K_Std, AIE_M_Std, AIE_C_Std, AIE_L_Std};
AIE_Postures_Means = {AIE_t1_Means, AIE_t2_Means, AIE_I_Means, AIE_D_Means, AIE_K_Means, AIE_M_Means, AIE_C_Means, AIE_L_Means};
AIE_Postures_Sens = {AIE_t1_Sens, AIE_t2_Sens, AIE_I_Sens, AIE_D_Sens, AIE_K_Sens, AIE_M_Sens, AIE_C_Sens, AIE_L_Sens};
AIE_Postures_Params_Matrix = {AIE_t1_Matrix, AIE_t2_Matrix, AIE_I_Matrix, AIE_D_Matrix, AIE_K_Matrix, AIE_M_Matrix, AIE_C_Matrix, AIE_L_Matrix};

save("mats/AIE_Postures_Stds.mat", "AIE_Postures_Stds");
save("mats/AIE_Postures_Means.mat", "AIE_Postures_Means");
save("mats/AIE_Postures_Sens.mat", "AIE_Postures_Sens");
save("mats/AIE_Postures_Params_Matrix.mat", "AIE_Postures_Params_Matrix");
