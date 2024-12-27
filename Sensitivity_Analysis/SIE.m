clear;
clear all;
disp("Welcome to the Single Input Excitation Sensitivity Analysis with Means");
disp("-----------------------------------------------------------");
disp("");
disp("Loading Parameters...");

%% Load the parameters
tic
Parameters;

%% Important values to load
% Load all the matrices that have the parameters of interest
CMT = {C, M, T1, T2};
IDK = {Ivar, Dvar, Kvar};
toc

% This index will remember where we are in the parameter matrices
p = 1;

% counts for T1 and T2 so that we only iterate over them once
% since we treat T1 and T2 as scalars instead of matrices
T1_count = 0;
T2_count = 0;

%% IDK params
disp("Performing analysis on parameters part of the I, D, and K matrices...");
tic
for i = 1:length(IDK)
    switch i
        case 1
            disp("I...");
        case 2
            disp("D...");
        case 3
            disp("K...");
    end
    for j = 1:size(IDK{i}, 2)
        % Define temporary vectors of parameters and their values, then
        % sub the values into the transfer function for all parameters
        % except the parameter of interest.
        tempmyvars = [myvars(1:p-1) myvars(p+1:end)];
        tempmyvals = [myvals(1:p-1) myvals(p+1:end)];

        % Get parameter value
        param = myvals(p);

        % substitute all params except that of param of interest
        diffG3 = subs(G3inv, tempmyvars, tempmyvals);
        diffG =  subs(diff(G(diffG3, M, C, T1, T2, w2), myvars(p)), myvars(p), param);
        for n = 1:length(w3)
            w = w3(n);
            % Substitute values into G3inv to calculate the transfer
            % function at that freqeuncy
            g3inv = double(subs(g3, w2, w));
            g = abs(G(g3inv, M, C, T1, T2, w));
    
            % Calculate diff and subsitute the differentiated parameter
            diffG = abs(double(subs(diffG, w2, w)));
            switch i
                case 1
                    I_Sensitivity{j, n} = diffG .* (param./g);
                case 2
                    D_Sensitivity{j, n} = diffG .* (param./g);
                case 3
                    K_Sensitivity{j, n} = diffG .* (param./g);
            end
        end
        
        switch i
            case 1
                I_S = I_Sensitivity{j, tremor_band(1):tremor_band(2)};
                I_Means{j} = mean(I_S, 2);
            case 2
                D_S = D_Sensitivity{j, tremor_band(1):tremor_band(2)};
                D_Means{j} = mean(D_S, 2);
            case 3
                K_S = K_Sensitivity{j, tremor_band(1):tremor_band(2)};
                K_Means{j} = mean(K_S, 2);
        end

        % Iterate over parameters
        p = p + 1;
    end
end
toc

%% Multiparameter sensitivity
p = 1;
disp("Performing analysis of C, M, T1, and T2 matrices...");
tic
for i=1:length(CMT)
    switch i
        case 1
            disp("C...");
        case 2
            disp("M...");
        case 3
            disp("T1...");
        case 4
            disp("T2...");
    end

    for j=1:size(CMT{i},1)
        for k=1:size(CMT{i},2)
            % Set the current matrix
            matrix = CMT{i};

            % Set up key variables such as the row and column of the matrix
            row = j;
            col = k;

            % Check first if the value is significant and do not include
            % situations where the row is greater than the column for I, D,
            % and K since they are symmetrical
            if matrix(row, col) == 0
                continue
            end

             % The parameter of interest is going to be set to the current
            % row and column of the matrix
            parameter = matrix(row, col);

            % Now we are going to create our forward increment for the
            % forward numerical differentiation where Xph = the matrix with the parameter(x)
            % plus h
            Xph = matrix;
            Xmh = Xph;
            dp = sqrt(eps) * parameter; % Make h as small as possible to increase accuracy

            % Make small change to parameter
            xph = parameter + dp;
            xmh = parameter - dp;
            h = dp * 2; % Define h

            % Put changed parameter into matrix
            Xph(row, col) = xph; 
            Xmh(row, col) = xmh;
            
            for n = 1:length(w3)
                w = w3(n); 
                g3inv = double(subs(g3, w2, w));
                if i == 1 % We are at C
                    C_Sensitivity{p, n} = (abs(G(g3inv, M, Xph, T1, T2, w)) - abs(G(g3inv, M, Xmh, T1, T2, w)))./h .* parameter./abs(G(g3inv, M, C, T1, T2, w));                
                end
    
                if i == 2 % We are at M
                    M_Sensitivity{p, n} = (abs(G(g3inv, Xph, C, T1, T2, w)) - abs(G(g3inv, Xmh, C, T1, T2, w)))./h .* parameter./abs(G(g3inv, M, C, T1, T2, w));                
                end
    
                if i == 3 && T1_count == 0 % We are at T1
                    T1_Sensitivity{p, n} = (abs(G(g3inv, M, C, T1 + T1*sqrt(eps), T2, w)) - abs(G(g3inv, M, C, T1 - T1*sqrt(eps), T2, w)))./h .* parameter./abs(G(g3inv, M, C, T1, T2, w));                
                    
                    % Make sure we do not iterate over T1 again
                    if n == 150
                        T1_count = 1;
                    end
                end
    
                if i == 4 && T2_count == 0 % We are at T2
                    T2_Sensitivity{p, n} = (abs(G(g3inv, M, C, T1, T2 + T2*sqrt(eps), w)) - abs(G(g3inv, M, C, T1, T2 - T2*sqrt(eps), w)))./h .* parameter./abs(G(g3inv, M, C, T1, T2, w));                
                    
                    % Make sure we do not iterate over T2 again
                    if n == 150
                        T2_count = 1;
                    end
                end
            end
            
            switch i
                case 1
                    C_S = C_Sensitivity{p, tremor_band(1):tremor_band(2)};
                    C_Means{p} = mean(C_S, 2);
                case 2
                    M_S = M_Sensitivity{p, tremor_band(1):tremor_band(2)};
                    M_Means{p} = mean(M_S, 2);
            end

            p = p + 1;
        end
    end
    p = 1;
end

% Get the means of T1 and T2 separately
T1_S = T1_Sensitivity{1, tremor_band(1):tremor_band(2)};
T2_S = T2_Sensitivity{1, tremor_band(1):tremor_band(2)};
T1_Means = {mean(T1_S, 2)};
T2_Means = {mean(T2_S, 2)};
toc

%% Combine all sensitivities real quick
disp("Saving matrices in mats folder:");
disp("- SIE_Sensitivities.mat");
disp("- SIE_Means.mat");
SIE_Sensitivities = {T1_Sensitivity, T2_Sensitivity, I_Sensitivity, D_Sensitivity, K_Sensitivity, M_Sensitivity, C_Sensitivity};
SIE_Means = {T1_Means, T2_Means, I_Means, D_Means, K_Means, M_Means, C_Means};

% Save the sensitvities
save("mats/SIE_Means.mat", "SIE_Means");
save("mats/SIE_Sensitivities.mat", "SIE_Sensitivities");