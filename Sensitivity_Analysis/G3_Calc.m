Parameters;

% Define symbolic variables
syms I11 I22 I32 I33 I41 I44 I52 I55 I62 I63 I66 I71 I74 I77 real
syms I23 I14 I25 I26 I36 I17 I47 real
I2 = [I11,   0,   0, I41,   0,   0, I71;
        0, I22, I32,   0, I52, I62,   0;
        0, I32, I33,   0,   0, I63,   0;
      I41,   0,   0, I44,   0,   0, I74;
        0, I52,   0,   0, I55,   0,   0;
        0, I62, I63,   0,   0, I66,   0;
      I71,   0,   0, I74,   0,   0, I77];

syms K11 K21 K22 K31 K32 K33 K41 K44 K55 K65 K66 K75 K76 K77 real
syms K12 K13 K23 K14 K56 K57 K67 real
K2 = [K11, K21, K31, K41,   0,   0,   0;
      K21, K22, K32,   0,   0,   0,   0;
      K31, K32, K33,   0,   0,   0,   0;
      K41,   0,   0, K44,   0,   0,   0;
        0,   0,   0,   0, K55, K65, K75;
        0,   0,   0,   0, K65, K66, K76;
        0,   0,   0,   0, K75, K76, K77];

syms D11 D21 D22 D31 D32 D33 D41 D44 D55 D65 D66 D75 D76 D77 real
syms D12 D13 D23 D14 D56 D57 D67 real
D2 = [D11, D21, D31, D41,   0,   0,   0;
      D21, D22, D32,   0,   0,   0,   0;
      D31, D32, D33,   0,   0,   0,   0;
      D41,   0,   0, D44,   0,   0,   0;
        0,   0,   0,   0, D55, D65, D75;
        0,   0,   0,   0, D65, D66, D76;
        0,   0,   0,   0, D75, D76, D77];

I_2 = [I11,   0,   0, I41,   0,   0, I71;
         0, I22, I32,   0, I52, I62,   0;
         0, I23, I33,   0,   0, I63,   0;
       I14,   0,   0, I44,   0,   0, I74;
         0, I25,   0,   0, I55,   0,   0;
         0, I26, I36,   0,   0, I66,   0;
       I17,   0,   0, I47,   0,   0, I77];

K_2 = [K11, K21, K31, K41,   0,   0,   0;
       K12, K22, K32,   0,   0,   0,   0;
       K13, K23, K33,   0,   0,   0,   0;
       K14,   0,   0, K44,   0,   0,   0;
         0,   0,   0,   0, K55, K65, K75;
         0,   0,   0,   0, K56, K66, K76;
         0,   0,   0,   0, K57, K67, K77];

D_2 = [D11, D21, D31, D41,   0,   0,   0;
       D12, D22, D32,   0,   0,   0,   0;
       D13, D23, D33,   0,   0,   0,   0;
       D14,   0,   0, D44,   0,   0,   0;
         0,   0,   0,   0, D55, D65, D75;
         0,   0,   0,   0, D56, D66, D76;
         0,   0,   0,   0, D57, D67, D77];

w2 = sym('w2','real');
%% Calc G3
disp("Calculating G3...");
tic
% G3 = -I2*w2^2 + D2*1i*w2 + K2;
G_3 = -I_2*w2^2 + D_2*1i*w2 + K_2;
toc

%% Calc inv of G3
disp("Calculating inv(G3)...");
tic
% invG3 = inv(-I2*w2^2 + D2*1i*w2 + K2);
invG_3 = inv(simplify(G_3)); % The simplify function should make things easier for matlab.
toc
save("invG_3.mat", "invG_3");
% g_3 = invG_3 * G_3;