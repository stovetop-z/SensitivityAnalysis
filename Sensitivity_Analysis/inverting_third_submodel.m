tic
m = 7;

%% SYMBOLIC INVERSION
% Create symbolic matrices
syms I11 I14 I17 I22 I23 I25 I26 I33 I36 I44 I47 I55 I66 I77

I_sym = [I11 0 0 I14 0 0 I17;
    0 I22 I23 0 I25 I26 0;
    0 I23 I33 0 0 I36 0;
    I14 0 0 I44 0 0 I47;
    0 I25 0 0 I55 0 0;
    0 I26 I36 0 0 I66 0;
    I17 0 0 I47 0 0 I77];

syms D11 D12 D13 D14 D22 D23 D33 D44 D55 D56 D57 D66 D67 D77

D_sym = [D11 D12 D13 D14 0 0 0;
    D12 D22 D23 0 0 0 0;
    D13 D23 D33 0 0 0 0;
    D14 0 0 D44 0 0 0;
    0 0 0 0 D55 D56 D57;
    0 0 0 0 D56 D66 D67;
    0 0 0 0 D57 D67 D77];

syms K11 K12 K13 K14 K22 K23 K33 K44 K55 K56 K57 K66 K67 K77

K_sym = [K11 K12 K13 K14 0 0 0;
    K12 K22 K23 0 0 0 0;
    K13 K23 K33 0 0 0 0;
    K14 0 0 K44 0 0 0;
    0 0 0 0 K55 K56 K57;
    0 0 0 0 K56 K66 K67;
    0 0 0 0 K57 K67 K77];

% Create symbolic G3
syms w_sym
G3_sym = I_sym(1:m,1:m)*(j*w_sym)^2 + D_sym(1:m,1:m)*(j*w_sym) + K_sym(1:m,1:m);

%% Invert symbolic G3
G3_sym_inv = inv(G3_sym);

% Evaluate quality of inverse
disp('G3_sym * G3_sym_inv = ')
G3_I_sym = G3_sym*G3_sym_inv;

%% NUMERICAL INVERSION
% Create numerical matrices
[I_num, D_num, K_num] = param();

% Create numerical G3
w_num = 5*(2*pi);
%%
G3_num = I_num(1:m,1:m)*(1i*w_num)^2 + D_num(1:m,1:m)*(1i*w_num) + K_num(1:m,1:m);

% Invert numerical G3
G3_num_inv = inv(G3_num);

% Evaluate quality of inverse
disp('G3_num * G3_num_inv = ')
disp(G3_num*G3_num_inv)

% COMPARISON B/W SYMBOLIC AND NUMERICAL INVERSES
% Assign values to elements of, and evaluate, G3_sym_inv
G3_sym_inv_num = subs(G3_sym_inv,...
    [I11 I14 I17 I22 I23 I25 I26 I33 I36 I44 I47 I55 I66 I77 ...
     D11 D12 D13 D14 D22 D23 D33 D44 D55 D56 D57 D66 D67 D77 ...
     K11 K12 K13 K14 K22 K23 K33 K44 K55 K56 K57 K66 K67 K77 ...
     w_sym],...
    [I_num(1,1) I_num(1,4) I_num(1,7) I_num(2,2) I_num(2,3) I_num(2,5) I_num(2,6) ...
     I_num(3,3) I_num(3,6) I_num(4,4) I_num(4,7) I_num(5,5) I_num(6,6) I_num(7,7) ...
     D_num(1,1) D_num(1,2) D_num(1,3) D_num(1,4) D_num(2,2) D_num(2,3) D_num(3,3) ...
     D_num(4,4) D_num(5,5) D_num(5,6) D_num(5,7) D_num(6,6) D_num(6,7) D_num(7,7) ...
     K_num(1,1) K_num(1,2) K_num(1,3) K_num(1,4) K_num(2,2) K_num(2,3) K_num(3,3) ...
     K_num(4,4) K_num(5,5) K_num(5,6) K_num(5,7) K_num(6,6) K_num(6,7) K_num(7,7) ...
     w_num]);
G3_sym_inv_num = eval(G3_sym_inv_num);

disp('Ratio of G3_sym_inv_num to G3_num_inv = ')
disp(G3_sym_inv_num./G3_num_inv)

toc

% load handel
% sound(y,Fs)