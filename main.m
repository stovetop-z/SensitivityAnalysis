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
sz = 1500;
A_rp = randperm(sz);
A_rp = (A_rp - 1) / (sz - 1);
A = reshape(A_rp, 100, 15);

theta_rp = randperm(sz);
theta_rp = (theta_rp - 1) / (sz - 1);
theta = reshape(theta_rp, 100, 15) * 2 * pi;

% Create matrix of inputs pre calculated to save computation time later
INPUT = A .* exp(1).^(-1i.*theta);

%% Frequency/Tremor band
frequency_band = (0.1:0.1:15);
w2 = 2 * pi * frequency_band;
tremor_band = [4, 8] * 10;

%% Step 1 - Compute Jacobian - SKIP IF DONE
% create_robot_arm;

%% Step 2 - Change specs of analysis
poses = ["Posture_1" "Posture_2" "Posture_3" "Posture_4" "Posture_5" "Posture_6" "Posture_7"];
postures = [1 2 3 4 5 6 7];

% Let's create constant's for the number of paramteres we are calculating
% per matrix
I_params = 24;
D_params = 14;
K_params = 14;
C_params = 15;
M_params = 48;
L_params = 3;

% the input_length is simply the number of inputs that are going into the
% system, while delta is a small value to create a small perturbation for
% a given parameter
input_length = size(INPUT, 1);
delta = sqrt(eps);

%% Step 3 - Run Analysis

% Load parameters and other necessary data
disp("Loading parameters...")
Parameters;

%% Loop through analysis with given postures
disp("Welcome! Running All Input Excitation Sensitivity Analysis...\n")
for p = 1:length(postures)    
    disp(poses(postures(p)));
    pos = postures(p);

    tic
    % Give all the values that change with posture a specific variable for consistency and ease of programming, and load all the matrices that have the parameters of interest
    j_p = J_p{pos};
    j_sym = J_sp{pos};
    I_pos = I{pos};
    M_pos = M{pos};
    Matrix = {I_pos, D, K, C, M_pos, t1, t2, L};
    
    % We need to first allocate memory for the sensitivities
    I_Sensitivity_3x1 = zeros(3, input_length, 150, I_params);
    I_Mean_Matrix = zeros(3, input_length, 40, I_params);
    
    D_Sensitivity_3x1 = zeros(3, input_length, 150, D_params);
    D_Mean_Matrix = zeros(3, input_length, 40, D_params);
    
    K_Sensitivity_3x1 = zeros(3, input_length, 150, K_params);
    K_Mean_Matrix = zeros(3, input_length, 40, K_params);
    
    C_Sensitivity_3x1 = zeros(3, input_length, 150, C_params);
    C_Mean_Matrix = zeros(3, input_length, 40, C_params);
    
    M_Sensitivity_3x1 = zeros(3, input_length, 150, M_params);
    M_Mean_Matrix = zeros(3, input_length, 40, M_params);
    
    t1_Sensitivity_3x1 = zeros(3, input_length, 150, 1);
    t1_Mean_Matrix = zeros(3, input_length, 40, 1);
    
    t2_Sensitivity_3x1 = zeros(3, input_length, 150, 1);
    t2_Mean_Matrix = zeros(3, input_length, 40, 1);
    
    L_Sensitivity_3x1 = zeros(3, input_length, 150, L_params);
    L_Mean_Matrix = zeros(3, input_length, 40, L_params);
    
    % This index will point to which parameter we are on in each matrix
    p = 1;
    
    % Counts for t1 and t2 so that we only iterate over them once
    t1_count = 0;
    t2_count = 0;
    
    %% Multiparameter sensitivity
    for i=1:length(Matrix)
        % Set the current matrix
        matrix = Matrix{i};
    
        for j=1:size(matrix, 1) % Row
            for k=1:size(matrix, 2) % Column
                % Set up key variables such as the row and column of the matrix
                row = j;
                col = k;
    
                % Check first if the parameter is a zero value for all matrices
                % except the I matrix, since we are computing some zero
                % parameters in I that are significant depending on posture.
                if i == 5 && row == 1 && col < 8
                elseif matrix(row, col) == 0 && i ~= 1
                    continue
                end
                
                % Now let's check if we are going on the left side of the I, D,
                % and K. If we are, let's skip since we assume these are symmetric
                % matrices.
                if (i == 1 || i == 2 || i == 3) && row > col
                    continue
                end
    
                % Finally, we can check if we are calculating on an appropriate
                % parameter while in I. This is due the fact that each posture
                % has it's own values of I parameters. And so, there are common
                % zero values for all postures that we want to skip. I45, I56, I57, I67
                if i == 1 
                    if (row == 4 && col == 5) || (row == 5 && col == 6) || (row == 5 && col == 7) || (row == 6 && col == 7)
                        continue
                    end

                    if matrix(row, col) == 0
                        for in = 1:size(INPUT, 2)
                            for n = 1:length(w2)
                                if n >= tremor_band(1) && n <= tremor_band(2)
                                    I_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = 0;
                                end

                                I_Sensitivity_3x1(:, in, n, p) = 0;
                            end
                        end

                        continue
                    end
                end
    
                % The parameter of interest is going to be set to the current
                % row and column of the matrix
                parameter = matrix(row, col);
    
                % Make x plus and minus dx for secant differentiation. Where
                X_plus_dx = matrix;
                X_min_dx = matrix;
                dx = delta * parameter; % Make dx as small as possible to increase accuracy with machine epsilon as delta
                x_plus_dx = parameter + dx;
                x_min_dx = parameter - dx;
                h = dx * 2; % Secant differentiation requries h to be multiplied by 2
                X_plus_dx(row, col) = x_plus_dx;
                X_min_dx(row, col) = x_min_dx;
                
                for in = 1:size(INPUT, 2)
                    input = INPUT(in, :);
    
                    for n = 1:length(w2)
                        w = w2(n);
    
                        if i == 1 % We are at I
                            X_plus_dx(col, row) = x_plus_dx;
                            X_min_dx(col, row) = x_min_dx;
                            I_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, X_plus_dx, D, K, M_pos, t1, t2, C, input, w)), 2) - sum(abs(Gx(j_p, X_min_dx, D, K, M_pos, t1, t2, C, input, w)), 2))./h;
                            
                            if n >= tremor_band(1) && n <= tremor_band(2)
                                I_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = I_Sensitivity_3x1(:, in, n, p);
                            end
                            
                            I_Sensitivity_3x1(:, in, n, p) = I_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);
                        end
            
                        if i == 2 % We are at D
                            X_plus_dx(col, row) = x_plus_dx;
                            X_min_dx(col, row) = x_min_dx;
                            D_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, X_plus_dx, K, M_pos, t1, t2, C, input, w)), 2) - sum(abs(Gx(j_p, I_pos, X_min_dx, K, M_pos, t1, t2, C, input, w)), 2))./h;    
        
                            if n >= tremor_band(1) && n <= tremor_band(2)
                                D_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = D_Sensitivity_3x1(:, in, n, p);
                            end
    
                            D_Sensitivity_3x1(:, in, n, p) = D_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);            
                        end
        
                        if i == 3 % We are at K
                            X_plus_dx(col, row) = x_plus_dx;
                            X_min_dx(col, row) = x_min_dx;
                            K_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, D, X_plus_dx, M_pos, t1, t2, C, input, w)), 2) - sum(abs(Gx(j_p, I_pos, D, X_min_dx, M_pos, t1, t2, C, input, w)), 2))./h;                
        
                            if n >= tremor_band(1) && n <= tremor_band(2)
                                K_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = K_Sensitivity_3x1(:, in, n, p);
                            end
    
                            K_Sensitivity_3x1(:, in, n, p) = K_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);   
                        end
            
                        if i == 4 % We are at C
                            C_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, X_plus_dx, input, w)), 2) - sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, X_min_dx, input, w)), 2))./h;
        
                            if n >= tremor_band(1) && n <= tremor_band(2)
                                C_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = C_Sensitivity_3x1(:, in, n, p);
                            end
    
                            C_Sensitivity_3x1(:, in, n, p) = C_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);   
                        end
            
                        if i == 5 % We are at M
                            M_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, D, K, X_plus_dx, t1, t2, C, input, w)), 2) - sum(abs(Gx(j_p, I_pos, D, K, X_min_dx, t1, t2, C, input, w)), 2))./h;
        
                            if n >= tremor_band(1) && n <= tremor_band(2)
                                M_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = M_Sensitivity_3x1(:, in, n, p);
                            end
    
                            M_Sensitivity_3x1(:, in, n, p) = M_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);   
                        end
            
                        if i == 6 && t1_count == 0 % We are at t1
                            t1_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1 + t1*delta, t2, C, input, w)), 2) - sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1 - t1*delta, t2, C, input, w)), 2))./h;
        
                            if n == 150
                                t1_count = 1;
                            end
    
                            if n >= tremor_band(1) && n <= tremor_band(2)
                                t1_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = t1_Sensitivity_3x1(:, in, n, p);
                            end
    
                            t1_Sensitivity_3x1(:, in, n, p) = t1_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);   
                        end
            
                        if i == 7 && t2_count == 0 % We are at t2
                            t2_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2 + t2*delta, C, input, w)), 2) - sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2 - t2*delta, C, input, w)), 2))./h;
        
                            if n == 150
                                t2_count = 1;
                            end
    
                            if n >= tremor_band(1) && n <= tremor_band(2)
                                t2_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = t2_Sensitivity_3x1(:, in, n, p);
                            end
    
                            t2_Sensitivity_3x1(:, in, n, p) = t2_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);   
                        end
        
                        if i == 8 % We are at L
                            L_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx_L(j_sym, X_plus_dx, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2) - sum(abs(Gx_L(j_sym, X_min_dx, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2))./h;    
    
                            if n >= tremor_band(1) && n <= tremor_band(2)
                                L_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = L_Sensitivity_3x1(:, in, n, p);
                            end
    
                            L_Sensitivity_3x1(:, in, n, p) = L_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);   
                        end
                    end
                end     
    
                % Increment 'p' so that we can proceed to next parameter in the
                % matrix
                p = p + 1;
            end
        end
    
    
        % New matrix, new set of parameters. Go back to first iteration
        p = 1;
    end
    
    % Calculate the magnitude of the 1st dimension
    % 3 x 100 (input_length) x 150 x parameter -> 1 x 100
    % (input_length) x 150 x parameter
    % Mean 3 x 100 (input_length) x 40 x parameter - > 1 x
    % 100 (input_length) x 40 x parameter
    I_Means_Matrix = sqrt(sum(I_Mean_Matrix).^2);
    I_Sensitivity_1x1 = sqrt(sum(I_Sensitivity_3x1).^2);
    
    D_Means_Matrix = sqrt(sum(D_Mean_Matrix).^2);
    D_Sensitivity_1x1 = sqrt(sum(D_Sensitivity_3x1).^2);
    
    K_Means_Matrix = sqrt(sum(K_Mean_Matrix).^2);
    K_Sensitivity_1x1 = sqrt(sum(K_Sensitivity_3x1).^2);
    
    C_Means_Matrix = sqrt(sum(C_Mean_Matrix).^2);
    C_Sensitivity_1x1 = sqrt(sum(C_Sensitivity_3x1).^2);
    
    M_Means_Matrix = sqrt(sum(M_Mean_Matrix).^2);
    M_Sensitivity_1x1 = sqrt(sum(M_Sensitivity_3x1).^2);
    
    t1_Means_Matrix = sqrt(sum(t1_Mean_Matrix).^2);
    t1_Sensitivity_1x1 = sqrt(sum(t1_Sensitivity_3x1).^2);
    
    t2_Means_Matrix = sqrt(sum(t2_Mean_Matrix).^2);
    t2_Sensitivity_1x1 = sqrt(sum(t2_Sensitivity_3x1).^2);
    
    L_Means_Matrix = sqrt(sum(L_Mean_Matrix).^2);
    L_Sensitivity_1x1 = sqrt(sum(L_Sensitivity_3x1).^2);
    
    % We have a 1 x 1 x 1 x 100 (input_length) x 150 double matrix. We .are going to take
    % the mean of the 4th dimension (so, all the VARIATIONS of
    % the input Amplitude and Theta), and leave the 3rd
    % dimension (all the inputs Amplitude and Theta, which are
    % 15) later for analysis.
    % Sensitivity goes from 1 x 100 (input_length) x 150 -> 150 x parameter
    % Means goes from 1 x 100 (input_length) x 40 x parameter -> 1 x parameter
    I_Sensitivity_1x1 = squeeze(squeeze(mean(I_Sensitivity_1x1, 2)));
    I_Means = mean(I_Sensitivity_1x1(40:80, :));
    I_Param_Matrix = squeeze(squeeze(mean(mean(I_Means_Matrix, 2), 3)));
    
    D_Sensitivity_1x1 = squeeze(squeeze(mean(D_Sensitivity_1x1, 2)));
    D_Means = mean(D_Sensitivity_1x1(40:80, :));
    D_Param_Matrix = squeeze(squeeze(mean(mean(D_Means_Matrix, 2), 3)));
    
    K_Sensitivity_1x1 = squeeze(squeeze(mean(K_Sensitivity_1x1, 2)));
    K_Means = mean(K_Sensitivity_1x1(40:80, :));
    K_Param_Matrix = squeeze(squeeze(mean(mean(K_Means_Matrix, 2), 3)));
    
    C_Sensitivity_1x1 = squeeze(squeeze(mean(C_Sensitivity_1x1, 2)));
    C_Means = mean(C_Sensitivity_1x1(40:80, :));
    C_Param_Matrix = squeeze(squeeze(mean(mean(C_Means_Matrix, 2), 3)));
    
    M_Sensitivity_1x1 = squeeze(squeeze(mean(M_Sensitivity_1x1, 2)));
    M_Means = mean(M_Sensitivity_1x1(40:80, :));
    M_Param_Matrix = squeeze(squeeze(mean(mean(M_Means_Matrix, 2), 3)));
    
    L_Sensitivity_1x1 = squeeze(squeeze(mean(L_Sensitivity_1x1, 2)));
    L_Means = mean(L_Sensitivity_1x1(40:80, :));
    L_Param_Matrix = squeeze(squeeze(mean(mean(L_Means_Matrix, 2), 3)));

    t1_Sensitivity_1x1 = squeeze(mean(t1_Sensitivity_1x1, 2));
    t1_Means = mean(t1_Sensitivity_1x1(40:80));
    t1_Param_Matrix = squeeze(squeeze(mean(mean(t1_Means_Matrix, 2), 3)));
    
    t2_Sensitivity_1x1 = squeeze(mean(t2_Sensitivity_1x1, 2));
    t2_Means = mean(t2_Sensitivity_1x1(40:80));
    t2_Param_Matrix = squeeze(squeeze(mean(mean(t2_Means_Matrix, 2), 3)));
    
    %% Combine all sensitivities real quick
    AIE_Sens = {t1_Sensitivity_1x1, t2_Sensitivity_1x1, I_Sensitivity_1x1, D_Sensitivity_1x1, K_Sensitivity_1x1, M_Sensitivity_1x1, C_Sensitivity_1x1, L_Sensitivity_1x1};
    AIE_Means = {t1_Means, t2_Means, I_Means, D_Means, K_Means, M_Means, C_Means, L_Means};
    AIE_Params_Matrix = {t1_Param_Matrix  t2_Param_Matrix  I_Param_Matrix  D_Param_Matrix  K_Param_Matrix  M_Param_Matrix  C_Param_Matrix  L_Param_Matrix};
    
    aie_m_pos_title = "mats/AIE_Means_Inputs_" + pos + ".mat";
    aie_p_pos_title = "mats/AIE_Params_Matrix_Inputs_" + pos + ".mat";
    aie_s_pos_title = "mats/AIE_Sens_Inputs_" + pos + ".mat";
    
    save(aie_m_pos_title, "AIE_Means");
    save(aie_p_pos_title, 'AIE_Params_Matrix');
    save(aie_s_pos_title, "AIE_Sens");
    toc
end

%% Analysis continued : computing mean of all postures
disp("Computing mean of all postures...");
tic;

AIE_I_Means = zeros(length(postures), 24);
AIE_D_Means = zeros(length(postures), 14);
AIE_K_Means = zeros(length(postures), 14);
AIE_M_Means = zeros(length(postures), 48);
AIE_C_Means = zeros(length(postures), 15);
AIE_t1_Means = zeros(length(postures), 1);
AIE_t2_Means = zeros(length(postures), 1);
AIE_L_Means = zeros(length(postures), 3);

for p = 1:length(postures)
    pos = postures(p);

    % Get pos mat location
    aie_p_pos_loc = "mats/AIE_Means_Inputs_" + pos + ".mat";
    aie_pos_loc = "mats/AIE_Sens_Inputs_" + pos + ".mat";
    aie_m_pos_loc = "mats/AIE_Params_Matrix_Inputs_" + pos + ".mat";

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
            AIE_L_Sens(p, :, :) = AIE_P_Sens;
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
AIE_M_Std(isnan(AIE_M_Std)) = 0;
AIE_M_Means = squeeze(mean(AIE_M_Means));
AIE_M_Means(isnan(AIE_M_Means)) = 0;
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
AIE_M_Sens(isnan(AIE_M_Sens)) = 0;
AIE_C_Sens = squeeze(mean(AIE_C_Sens));
AIE_t1_Sens = squeeze(mean(AIE_t1_Sens));
AIE_t2_Sens = squeeze(mean(AIE_t2_Sens));
AIE_L_Sens = squeeze(mean(AIE_L_Sens));

AIE_I_Matrix = squeeze(mean(AIE_I_Matrix));
AIE_D_Matrix = squeeze(mean(AIE_D_Matrix));
AIE_K_Matrix = squeeze(mean(AIE_K_Matrix));
AIE_M_Matrix = squeeze(mean(AIE_M_Matrix));
AIE_M_Matrix(isnan(AIE_M_Matrix)) = 0;
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
