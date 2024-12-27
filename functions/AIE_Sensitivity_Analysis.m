function [] = AIE_Sensitivity_Analysis(I_pos, D, K, C, M_pos, t1, t2, L, j_sym, j_p, INPUT, w2, posture)

addpath("functions/");
% Load all the matrices that have the parameters of interest
Matrix = {I_pos, D, K, C, M_pos, t1, t2, L};

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
tremor_band = [4, 8] * 10;
delta = eps;

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

% counts for t1 and t2 so that we only iterate over them once
t1_count = 0;
t2_count = 0;

%% Multiparameter sensitivity
for i=1:length(Matrix)
    for j=1:size(Matrix{i},1)
        for k=1:size(Matrix{i},2)
            % Set the current matrix
            matrix = Matrix{i};

            % Set up key variables such as the row and column of the matrix
            row = j;
            col = k;

            % Check first if the inue is a zero value for all matrices
            % except the I matrix, since we are computing some zero
            % parameters.
            if matrix(row, col) == 0 && i ~= 1
                continue
            end
            
            % Now let's check if we are going on the left side of the I, D,
            % and K. If we are, let's skip.
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
            end

             % The parameter of interest is going to be set to the current
            % row and column of the matrix
            parameter = matrix(row, col);

            % Now we are going to create our forward increment for the
            % forward numerical differentiation where Xph = the matrix with the parameter(x)
            % plus h
            Xph = matrix;
            Xmh = matrix;
            dp = delta * parameter; % Make h as small as possible to increase accuracy with machine epsilon as delta
            xph = parameter + dp;
            xmh = parameter - dp;
            h = dp * 2;
            Xph(row, col) = xph;
            Xmh(row, col) = xmh;
            
            for in = 1:size(INPUT, 2)
                input = INPUT(in, :);

                for n = 1:length(w2)
                    w = w2(n);

                    if i == 1 % We are at I
                        Xph(col, row) = xph;
                        Xmh(col, row) = xmh;
                        I_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, Xph, D, K, M_pos, t1, t2, C, input, w)), 2) - sum(abs(Gx(j_p, Xmh, D, K, M_pos, t1, t2, C, input, w)), 2))./h;
                        
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            I_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = I_Sensitivity_3x1(:, in, n, p);
                        end
                        
                        I_Sensitivity_3x1(:, in, n, p) = I_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);
                    end
        
                    if i == 2 % We are at D
                        Xph(col, row) = xph;
                        Xmh(col, row) = xmh;
                        D_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, Xph, K, M_pos, t1, t2, C, input, w)), 2) - sum(abs(Gx(j_p, I_pos, Xmh, K, M_pos, t1, t2, C, input, w)), 2))./h;    
    
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            D_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = D_Sensitivity_3x1(:, in, n, p);
                        end

                        D_Sensitivity_3x1(:, in, n, p) = D_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);            
                    end
    
                    if i == 3 % We are at K
                        Xph(col, row) = xph;
                        Xmh(col, row) = xmh;
                        K_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, D, Xph, M_pos, t1, t2, C, input, w)), 2) - sum(abs(Gx(j_p, I_pos, D, Xmh, M_pos, t1, t2, C, input, w)), 2))./h;                
    
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            K_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = K_Sensitivity_3x1(:, in, n, p);
                        end

                        K_Sensitivity_3x1(:, in, n, p) = K_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);   
                    end
        
                    if i == 4 % We are at C
                        C_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, Xph, input, w)), 2) - sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, Xmh, input, w)), 2))./h;
    
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            C_Mean_Matrix(:, in, n - tremor_band(1) + 1, p) = C_Sensitivity_3x1(:, in, n, p);
                        end

                        C_Sensitivity_3x1(:, in, n, p) = C_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx(j_p, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2);   
                    end
        
                    if i == 5 % We are at M
                        M_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(j_p, I_pos, D, K, Xph, t1, t2, C, input, w)), 2) - sum(abs(Gx(j_p, I_pos, D, K, Xmh, t1, t2, C, input, w)), 2))./h;
    
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
                        L_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx_L(j_sym, Xph, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2) - sum(abs(Gx_L(j_sym, Xmh, I_pos, D, K, M_pos, t1, t2, C, input, w)), 2))./h;    

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

t1_Sensitivity_1x1 = squeeze(mean(t1_Sensitivity_1x1, 2));
t1_Means = mean(t1_Sensitivity_1x1(40:80));
t1_Param_Matrix = squeeze(squeeze(mean(mean(t1_Means_Matrix, 2), 3)));

t2_Sensitivity_1x1 = squeeze(mean(t2_Sensitivity_1x1, 2));
t2_Means = mean(t2_Sensitivity_1x1(40:80));
t2_Param_Matrix = squeeze(squeeze(mean(mean(t2_Means_Matrix, 2), 3)));

L_Sensitivity_1x1 = squeeze(squeeze(mean(L_Sensitivity_1x1, 2)));
L_Means = mean(L_Sensitivity_1x1(40:80, :));
L_Param_Matrix = squeeze(squeeze(mean(mean(L_Means_Matrix, 2), 3)));

%% Combine all sensitivities real quick
AIE_Sens = {t1_Sensitivity_1x1, t2_Sensitivity_1x1, I_Sensitivity_1x1, D_Sensitivity_1x1, K_Sensitivity_1x1, M_Sensitivity_1x1, C_Sensitivity_1x1, L_Sensitivity_1x1};
AIE_Means = {t1_Means, t2_Means, I_Means, D_Means, K_Means, M_Means, C_Means, L_Means};
AIE_Params_Matrix = {t1_Param_Matrix  t2_Param_Matrix  I_Param_Matrix  D_Param_Matrix  K_Param_Matrix  M_Param_Matrix  C_Param_Matrix  L_Param_Matrix};

aie_m_pos_title = "mats/AIE_Means_Inputs_" + posture + ".mat";
aie_p_pos_title = "mats/AIE_Params_Matrix_Inputs_" + posture + ".mat";
aie_s_pos_title = "mats/AIE_Sens_Inputs_" + posture + ".mat";

save(aie_m_pos_title, "AIE_Means");
save(aie_p_pos_title, 'AIE_Params_Matrix');
save(aie_s_pos_title, "AIE_Sens");
end

