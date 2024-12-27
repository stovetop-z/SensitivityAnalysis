Parameters;

Matrix = {I, D, K, C, M, t1, t2, L};
input_length = size(INPUT, 1);
tremor_band = [4, 8] * 10;
tb_length = tremor_band(2) - tremor_band(1);
delta = sqrt(eps);

% We need to first allocate memory for the sensitivities
I_Org = zeros(7, 7);
I_Sensitivity_3x1 = zeros(3, input_length, 150, 14);
I_Mean = zeros(3, input_length, 40, 14);

D_Org = zeros(7 ,7);
D_Sensitivity_3x1 = zeros(3, input_length, 150, 14);
D_Mean = zeros(3, input_length, 40, 14);

K_Org = zeros(7 ,7);
K_Sensitivity_3x1 = zeros(3, input_length, 150, 14);
K_Mean = zeros(3, input_length, 40, 14);

C_Org = zeros(15);
C_Sensitivity_3x1 = zeros(3, input_length, 150, 15);
C_Mean = zeros(3, input_length, 40, 15);

M_Org = zeros(7, 15);
M_Sensitivity_3x1 = zeros(3, input_length, 150, 48);
M_Mean = zeros(3, input_length, 40, 48);

t1_Sensitivity_3x1 = zeros(3, input_length, 150, 1);
t1_Mean = zeros(3, input_length, 40, 1);

t2_Sensitivity_3x1 = zeros(3, input_length, 150, 1);
t2_Mean = zeros(3, input_length, 40, 1);

L_Sensitivity_3x1 = zeros(3, input_length, 150, 3);
L_Mean = zeros(3, input_length, 40, 3);

% This index will point to which parameter we are on in each matrix
p = 1;
num = 1;

% counts for t1 and t2 so that we only iterate over them once
t1_count = 0;
t2_count = 0;

%% Multiparameter sensitivity
tic
for i=1:length(Matrix)
    for j=1:size(Matrix{i},1)
        for k=1:size(Matrix{i},2)
            % Set the current matrix
            matrix = Matrix{i};
            disp(num/110 * 100);

            % Set up key variables such as the row and column of the matrix
            row = j;
            col = k;

            % Check first if the inue is significant and do not include
            % situations where the row is greater than the column for I, D,
            % and K since they are symmetrical
            if matrix(row, col) == 0 || ((i == 1 || i == 2 || i == 3) && row > col)
                continue
            end

             % The parameter of interest is going to be set to the current
            % row and column of the matrix
            parameter = matrix(row, col);

            % Now we are going to create our forward increment for the
            % forward numerical differentiation where Xph = the matrix with the parameter(x)
            % plus h
            Xph = matrix;
            Xmh = matrix;
            dp = delta * parameter; % Make h as small as possible to increase accuracy
            xph = dp + parameter;
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
                        I_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(J_p1, L, Xph, D, K, M, t1, t2, C, input, w)), 2) - sum(abs(Gx(J_p1, L, Xmh, D, K, M, t1, t2, C, input, w)), 2))./h;
                        
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            I_Mean(:, in, n - tremor_band(1) + 1, p) = I_Sensitivity_3x1(:, in, n, p);
                        end
                        
                        I_Sensitivity_3x1(:, in, n, p) = I_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx_L(J_sp1, L, I, D, K, M, t1, t2, C, input, w)), 2);
                    end
        
                    if i == 2 % We are at D
                        Xph(col, row) = xph;
                        Xmh(col, row) = xmh;
                        D_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(J_p1, L, I, Xph, K, M, t1, t2, C, input, w)), 2) - sum(abs(Gx(J_p1, L, I, Xmh, K, M, t1, t2, C, input, w)), 2))./h;    
    
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            D_Mean(:, in, n - tremor_band(1) + 1, p) = D_Sensitivity_3x1(:, in, n, p);
                        end

                        D_Sensitivity_3x1(:, in, n, p) = D_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx_L(J_sp1, L, I, D, K, M, t1, t2, C, input, w)), 2);            
                    end
    
                    if i == 3 % We are at K
                        Xph(col, row) = xph;
                        Xmh(col, row) = xmh;
                        K_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(J_p1, L, I, D, Xph, M, t1, t2, C, input, w)), 2) - sum(abs(Gx(J_p1, L, I, D, Xmh, M, t1, t2, C, input, w)), 2))./h;                
    
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            K_Mean(:, in, n - tremor_band(1) + 1, p) = K_Sensitivity_3x1(:, in, n, p);
                        end

                        K_Sensitivity_3x1(:, in, n, p) = K_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx_L(J_sp1, L, I, D, K, M, t1, t2, C, input, w)), 2);   
                    end
        
                    if i == 4 % We are at C
                        C_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(J_p1, L, I, D, K, M, t1, t2, Xph, input, w)), 2) - sum(abs(Gx(J_p1, L, I, D, K, M, t1, t2, Xmh, input, w)), 2))./h;
    
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            C_Mean(:, in, n - tremor_band(1) + 1, p) = C_Sensitivity_3x1(:, in, n, p);
                        end

                        C_Sensitivity_3x1(:, in, n, p) = C_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx_L(J_sp1, L, I, D, K, M, t1, t2, C, input, w)), 2);   
                    end
        
                    if i == 5 % We are at M
                        M_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(J_p1, L, I, D, K, Xph, t1, t2, C, input, w)), 2) - sum(abs(Gx(J_p1, L, I, D, K, Xmh, t1, t2, C, input, w)), 2))./h;
    
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            M_Mean(:, in, n - tremor_band(1) + 1, p) = M_Sensitivity_3x1(:, in, n, p);
                        end

                        M_Sensitivity_3x1(:, in, n, p) = M_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx_L(J_sp1, L, I, D, K, M, t1, t2, C, input, w)), 2);   
                    end
        
                    if i == 6 && t1_count == 0 % We are at t1
                        t1_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(J_p1, L, I, D, K, M, t1 + t1*delta, t2, C, input, w)), 2) - sum(abs(Gx(J_p1, L, I, D, K, M, t1 - t1*delta, t2, C, input, w)), 2))./h;
    
                        if n == 150
                            t1_count = 1;
                        end

                        if n >= tremor_band(1) && n <= tremor_band(2)
                            t1_Mean(:, in, n - tremor_band(1) + 1, p) = t1_Sensitivity_3x1(:, in, n, p);
                        end

                        t1_Sensitivity_3x1(:, in, n, p) = t1_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx_L(J_sp1, L, I, D, K, M, t1, t2, C, input, w)), 2);   
                    end
        
                    if i == 7 && t2_count == 0 % We are at t2
                        t2_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(J_p1, L, I, D, K, M, t1, t2 + t2*delta, C, input, w)), 2) - sum(abs(Gx(J_p1, L, I, D, K, M, t1, t2 - t2*delta, C, input, w)), 2))./h;
    
                        if n == 150
                            t2_count = 1;
                        end

                        if n >= tremor_band(1) && n <= tremor_band(2)
                            t2_Mean(:, in, n - tremor_band(1) + 1, p) = t2_Sensitivity_3x1(:, in, n, p);
                        end

                        t2_Sensitivity_3x1(:, in, n, p) = t2_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx_L(J_sp1, L, I, D, K, M, t1, t2, C, input, w)), 2);   
                    end
    
                    if i == 8 % We are at L
                        L_Sensitivity_3x1(:, in, n, p) = (sum(abs(Gx(J_p1, Xph, I, D, K, M, t1, t2, C, input, w)), 2) - sum(abs(Gx(J_p1, Xmh, I, D, K, M, t1, t2, C, input, w)), 2))./h;
                        
                        if n >= tremor_band(1) && n <= tremor_band(2)
                            L_Mean(:, in, n - tremor_band(1) + 1, p) = L_Sensitivity_3x1(:, in, n, p);
                        end

                        L_Sensitivity_3x1(:, in, n, p) = L_Sensitivity_3x1(:, in, n, p) .* parameter./sum(abs(Gx_L(J_sp1, L, I, D, K, M, t1, t2, C, input, w)), 2);   
                    end
                end
            end     

            % Increment 'm' so that we can proceed to next parameter in the
            % matrix
            p = p + 1;
            num = num + 1;
        end
    end


    % New matrix, new set of parameters. Go back to first iteration
    num = num - 1;
    p = 1;
end

% Calculate the magnitude of the 1st dimension
% 3 x 100 (input_length) x 150 x parameter -> 1 x 100
% (input_length) x 150 x parameter
% Mean 3 x 100 (input_length) x 40 x parameter - > 1 x
% 100 (input_length) x 40 x parameter
I_Means = sqrt(sum(I_Mean).^2);
I_Sensitivity_1x1 = sqrt(sum(I_Sensitivity_3x1).^2);

D_Means = sqrt(sum(D_Mean).^2);
D_Sensitivity_1x1 = sqrt(sum(D_Sensitivity_3x1).^2);

K_Means = sqrt(sum(K_Mean).^2);
K_Sensitivity_1x1 = sqrt(sum(K_Sensitivity_3x1).^2);

C_Means = sqrt(sum(C_Mean).^2);
C_Sensitivity_1x1 = sqrt(sum(C_Sensitivity_3x1).^2);

M_Means = sqrt(sum(M_Mean).^2);
M_Sensitivity_1x1 = sqrt(sum(M_Sensitivity_3x1).^2);

t1_Means = sqrt(sum(t1_Mean).^2);
t1_Sensitivity_1x1 = sqrt(sum(t1_Sensitivity_3x1).^2);

t2_Means = sqrt(sum(t2_Mean).^2);
t2_Sensitivity_1x1 = sqrt(sum(t2_Sensitivity_3x1).^2);

L_Means = sqrt(sum(L_Mean).^2);
L_Sensitivity_1x1 = sqrt(sum(L_Sensitivity_3x1).^2);

% We have a 1 x 1 x 1 x 100 (input_length) x 150 double matrix. We .are going to take
% the mean of the 4th dimension (so, all the VARIATIONS of
% the input Amplitude and Theta), and leave the 3rd
% dimension (all the inputs Amplitude and Theta, which are
% 15) later for analysis.
% Sensitivity goes from 1 x 100 (input_length) x 150 -> 150 x parameter
% Means goes from 1 x 100 (input_length) x 40 x parameter -> 1 x parameter
I_Sensitivity_1x1 = squeeze(squeeze(mean(I_Sensitivity_1x1, 2)));
I_Mean = squeeze(squeeze(mean(mean(I_Means, 2), 3)));

D_Sensitivity_1x1 = squeeze(squeeze(mean(D_Sensitivity_1x1, 2)));
D_Mean = squeeze(squeeze(mean(mean(D_Means, 2), 3)));

K_Sensitivity_1x1 = squeeze(squeeze(mean(K_Sensitivity_1x1, 2)));
K_Mean = squeeze(squeeze(mean(mean(K_Means, 2), 3)));

C_Sensitivity_1x1 = squeeze(squeeze(mean(C_Sensitivity_1x1, 2)));
C_Mean = squeeze(squeeze(mean(mean(C_Means, 2), 3)));

M_Sensitivity_1x1 = squeeze(squeeze(mean(M_Sensitivity_1x1, 2)));
M_Mean = squeeze(squeeze(mean(mean(M_Means, 2), 3)));

t1_Sensitivity_1x1 = squeeze(mean(t1_Sensitivity_1x1, 2));
t1_Mean = squeeze(squeeze(mean(mean(t1_Means, 2), 3)));

t2_Sensitivity_1x1 = squeeze(mean(t2_Sensitivity_1x1, 2));
t2_Mean = squeeze(squeeze(mean(mean(t2_Means, 2), 3)));

L_Sensitivity_1x1 = squeeze(squeeze(mean(L_Sensitivity_1x1, 2)));
L_Mean = squeeze(squeeze(mean(mean(L_Means, 2), 3)));
toc

%% Combine all sensitivities real quick
AIE_X_Sens = {t1_Sensitivity_1x1, t2_Sensitivity_1x1, I_Sensitivity_1x1, D_Sensitivity_1x1, K_Sensitivity_1x1, M_Sensitivity_1x1, C_Sensitivity_1x1, L_Sensitivity_1x1};
AIE_X_Means = {t1_Mean, t2_Mean, I_Mean, D_Mean, K_Mean, M_Mean, C_Mean, L_Mean};

aie_m_pos_title = "mats/AIE_X_Means_15_Inputs.mat";
aie_s_pos_title = "mats/AIE_X_Sens_15_Inputs.mat";

save(aie_m_pos_title, 'AIE_X_Means');
save(aie_s_pos_title, "AIE_X_Sens");

