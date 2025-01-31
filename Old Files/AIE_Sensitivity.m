clear;
clear all;

%% Load the parameters
Parameters;

%% Important values to load
% Load all the matrices that have the parameters of interest
Matrix = {I, D, K, C, M, t1, t2};

% This indices will remember where we are in the sensitivity matrices
m = 1;

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

            % Set up key variables such as the row and column of the matrix
            row = j;
            col = k;

            % Check first if the value is significant and do not include
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
            dp = sqrt(eps) * parameter; % Make h as small as possible to increase accuracy
            xph = dp + parameter;
            xmh = parameter - dp;
            h = 2 * dp;
            Xph(row, col) = xph;
            Xmh(row, col) = xmh;
            
            for n = 1:length(w2)
                if i == 1 % We are at I
                    Xph(col, row) = xph;
                    Xmh(col, row) = xmh;
                    I_Sensitivity{m, n} = (abs(sum(G_Paul(Xph, D, K, M, t1, t2, C, w2(n)), 02)) - abs(sum(G_Paul(Xmh, D, K, M, t1, t2, C, w2(n)), 02)))./h .* parameter./abs(sum(G_Paul(I, D, K, M, t1, t2, C, w2(n)), 02));
                
                end
    
                if i == 2 % We are at D
                    Xph(col, row) = xph;
                    Xmh(col, row) = xmh;
                    D_Sensitivity{m, n} = (abs(sum(G_Paul(I, Xph, K, M, t1, t2, C, w2(n)), 02)) - abs(sum(G_Paul(I, Xmh, K, M, t1, t2, C, w2(n)), 02)))./h .* parameter./abs(sum(G_Paul(I, D, K, M, t1, t2, C, w2(n)), 02));                
                end

                if i == 3 % We are at K
                    Xph(col, row) = xph;
                    Xmh(col, row) = xmh;
                    K_Sensitivity{m, n} = (abs(sum(G_Paul(I, D, Xph, M, t1, t2, C, w2(n)), 02)) - abs(sum(G_Paul(I, D, Xmh, M, t1, t2, C, w2(n)), 02)))./h .* parameter./abs(sum(G_Paul(I, D, K, M, t1, t2, C, w2(n)), 02));                
                end
    
                if i == 4 % We are at C
                    C_Sensitivity{m, n} = (abs(sum(G_Paul(I, D, K, M, t1, t2, Xph, w2(n)), 02)) - abs(sum(G_Paul(I, D, K, M, t1, t2, Xmh, w2(n)), 02)))./h .* parameter./abs(sum(G_Paul(I, D, K, M, t1, t2, C, w2(n)), 02));                
                end
    
                if i == 5 % We are at M
                    M_Sensitivity{m, n} = (abs(sum(G_Paul(I, D, K, Xph, t1, t2, C, w2(n)), 02)) - abs(sum(G_Paul(I, D, K, Xmh, t1, t2, C, w2(n)), 02)))./h .* parameter./abs(sum(G_Paul(I, D, K, M, t1, t2, C, w2(n)), 02));                
                end
    
                if i == 6 && t1_count == 0 % We are at t1
                    t1_Sensitivity{m, n} = (abs(sum(G_Paul(I, D, K, M, t1 + t1*sqrt(eps), t2, C, w2(n)), 02)) - abs(sum(G_Paul(I, D, K, M, t1 - t1*sqrt(eps), t2, C, w2(n)), 02)))./h .* parameter./abs(sum(G_Paul(I, D, K, M, t1, t2, C, w2(n)), 02));                
                    if n == 150
                        t1_count = 1;
                    end
                end
    
                if i == 7 && t2_count == 0 % We are at t2
                    t2_Sensitivity{m, n} = (abs(sum(G_Paul(I, D, K, M, t1, t2 + t2*sqrt(eps), C, w2(n)), 02)) - abs(sum(G_Paul(I, D, K, M, t1, t2 - t2*sqrt(eps), C, w2(n)), 02)))./h .* parameter./abs(sum(G_Paul(I, D, K, M, t1, t2, C, w2(n)), 02));                
                    if n == 150
                        t2_count = 1;
                    end
                end
            end

            % Get the means within the specified tremor_band in the
            % Parameters.m
            switch i
                case 1
                    I_S = I_Sensitivity{m, tremor_band(1):tremor_band(2)};
                    I_Means{m} = mean(I_S, 2);
                case 2
                    D_S = D_Sensitivity{m, tremor_band(1):tremor_band(2)};
                    D_Means{m} = mean(D_S, 2);
                case 3
                    K_S = K_Sensitivity{m, tremor_band(1):tremor_band(2)};
                    K_Means{m} = mean(K_S, 2);
                case 4
                    C_S = C_Sensitivity{m, tremor_band(1):tremor_band(2)};
                    C_Means{m} = mean(C_S, 2);
                case 5
                    M_S = M_Sensitivity{m, tremor_band(1):tremor_band(2)};
                    M_Means{m} = mean(M_S, 2);
            end

            m = m + 1;
        end
    end
    m = 1;
end

% Get the means of t1 and t2 separately
t1_S = t1_Sensitivity{1, tremor_band(1):tremor_band(2)};
t2_S = t2_Sensitivity{1, tremor_band(1):tremor_band(2)};
t1_Means = {mean(t1_S, 2)};
t2_Means = {mean(t2_S, 2)};
toc

%% Combine all sensitivities real quick
AIE_Sensitivities = {t1_Sensitivity, t2_Sensitivity, I_Sensitivity, D_Sensitivity, K_Sensitivity, M_Sensitivity, C_Sensitivity};
AIE_Means = {t1_Means, t2_Means, I_Means, D_Means, K_Means, M_Means, C_Means};

save("mats/AIE_Means.mat", "AIE_Means");
save("mats/AIE_Sensitivities.mat", "AIE_Sensitivities");

%% TODO
% - Change to PICTURE
% -