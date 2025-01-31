addpath("functions/");
addpath("Per50AverageGenderParams/");

% Load the parameters that are from each of the 7 postures
% I (Inertia), M (Moment arms)
load("P1_Per50_average_Params.mat");
I_P1 = I;
M_P1 = M(:, [1 2 3 10 21 22 16 17 23 24 31 28 29 26 27]);
q_p1 = q_curr;
load("P2_Per50_average_Params.mat");
I_P2 = I;
M_P2 = M(:, [1 2 3 10 21 22 16 17 23 24 31 28 29 26 27]);
q_p2 = q_curr;
load("P3_Per50_average_Params.mat");
I_P3 = I;
M_P3 = M(:, [1 2 3 10 21 22 16 17 23 24 31 28 29 26 27]);
q_p3 = q_curr;
load("P4_Per50_average_Params.mat");
I_P4 = I;
M_P4 = M(:, [1 2 3 10 21 22 16 17 23 24 31 28 29 26 27]);
q_p4 = q_curr;
load("P5_Per50_average_Params.mat");
I_P5 = I;
M_P5 = M(:, [1 2 3 10 21 22 16 17 23 24 31 28 29 26 27]);
q_p5 = q_curr;
load("P6_Per50_average_Params.mat");
I_P6 = I;
M_P6 = M(:, [1 2 3 10 21 22 16 17 23 24 31 28 29 26 27]);
q_p6 = q_curr;
load("P7_Per50_average_Params.mat");
I_P7 = I;
M_P7 = M(:, [1 2 3 10 21 22 16 17 23 24 31 28 29 26 27]);
q_p7 = q_curr;

% The I parameter matrices have zeroes that are in common to ALL of them.
% They are listed below in row, column relation. 
% I45, I56, I57, I67

I = {I_P1, I_P2, I_P3, I_P4, I_P5, I_P6, I_P7};
M = {M_P1, M_P2, M_P3, M_P4, M_P5, M_P6, M_P7};

%% t1 (Muscle time constant 1), t2 (Muscle time constant 2), C (Max muscle force)
t1 = eye(15).*30e-3; %s
t2 = eye(15).*40e-3; %s

C = diag([1218.9, 1103.5, 201.6, 658.3, 525.1, 316.8, 771.8, 717.5,...
                       1177.4,   276, 557.2, 407.9, 479.8, 589.8, 192.9]);  % Newtons

%% D (Damping), K (Stiffness)

D =     [0.756, 0.184, 0.020, 0.187,     0,      0,      0;... % Nms/rad
          0.184, 0.383, 0.267,     0,     0,      0,      0;...
          0.020, 0.267, 0.524,     0,     0,      0,      0;...
          0.187,     0,     0, 0.607,     0,      0,      0;...
              0,     0,     0,     0, 0.021,  0.001,  0.008;...
              0,     0,     0,     0, 0.001,  0.028, -0.003;...
              0,     0,     0,     0, 0.008, -0.003,  0.082];

K =     [10.80, 2.626, 0.279, 2.670,     0,      0,      0;... % Nm/rad
          2.626, 5.468, 3.821,     0,     0,      0,      0;...
          0.279, 3.821, 7.486,     0,     0,      0,      0;...
          2.670,     0,     0, 8.670,     0,      0,      0;...
              0,     0,     0,     0, 0.756,  0.018,  0.291;
              0,     0,     0,     0, 0.018,  0.992, -0.099;...
              0,     0,     0,     0, 0.291, -0.099,  2.920];

%% Jacobian parameters and variables here
load("mats/J.mat", "J");

%% LINK LENGTHS
% - Length upper arm (Lua)
% - Length forearm (Lfa)
% - Distance from axis (wrist) to center of mass of hand (Lh)
male_hand = 0.069;
male_forearm = 0.2723;
male_upperarm = 0.2853;
female_upperarm = 0.2558;
female_forearm = .2457;
female_hand = .0725;

h = mean(male_hand + female_hand);
f = mean(male_forearm + female_forearm);
u = mean(male_upperarm + female_upperarm);

L = [h, f, u];

% Position variable for posture 1
% Check: Fundamental Principles of Tremor Propagation in the Upper Limb by
% Davidson and Charles at Fig 1 for values for different postures
syms q1 q2 q3 q4 q5 q6 q7 Lh Lfa Lua

% Sym Jacobians from postures 1-7
J_sp1 = subs(J, [q1 q2 q3 q4 q5 q6 q7], q_p1);
J_sp2 = subs(J, [q1 q2 q3 q4 q5 q6 q7], q_p2);
J_sp3 = subs(J, [q1 q2 q3 q4 q5 q6 q7], q_p3);
J_sp4 = subs(J, [q1 q2 q3 q4 q5 q6 q7], q_p4);
J_sp5 = subs(J, [q1 q2 q3 q4 q5 q6 q7], q_p5);
J_sp6 = subs(J, [q1 q2 q3 q4 q5 q6 q7], q_p6);
J_sp7 = subs(J, [q1 q2 q3 q4 q5 q6 q7], q_p7);
J_sp = {J_sp1 J_sp2 J_sp3 J_sp4 J_sp5 J_sp6 J_sp7};

% Double Jacobians from postures 1-7
J_p1 = double(subs(J_sp1, [Lh, Lfa, Lua], L));
J_p2 = double(subs(J_sp2, [Lh, Lfa, Lua], L));
J_p3 = double(subs(J_sp3, [Lh, Lfa, Lua], L));
J_p4 = double(subs(J_sp4, [Lh, Lfa, Lua], L));
J_p5 = double(subs(J_sp5, [Lh, Lfa, Lua], L));
J_p6 = double(subs(J_sp6, [Lh, Lfa, Lua], L));
J_p7 = double(subs(J_sp7, [Lh, Lfa, Lua], L));
J_p = {J_p1 J_p2 J_p3 J_p4 J_p5 J_p6 J_p7};

%% array to hold the values as chars symbolically for plotting
params_char = {'t1', 't2', 'I11', 'I21', 'I31', 'I41', 'I51', 'I61', 'I71', 'I22', 'I32', 'I42', 'I52', 'I62', 'I72', 'I33', 'I43', 'I53', 'I63', 'I73', 'I44', 'I64', 'I74', 'I55', 'I66', 'I77', 'D11', 'D21', 'D31', 'D41', 'D22', 'D32', 'D33', 'D44', 'D55', 'D65', 'D75', 'D66', 'D76', 'D77', 'K11', 'K21', 'K31', 'K41', 'K22', 'K32', 'K33', 'K44', 'K55', 'K65', 'K75', 'K66', 'K76', 'K77', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M31', 'M32', 'M33', 'M34', 'M35', 'M36', 'M37', 'M45', 'M46', 'M47', 'M48', 'M49', 'M410', 'M411', 'M412', 'M413', 'M414', 'M415', 'M55', 'M56', 'M510', 'M511', 'M512', 'M513', 'M514', 'M515', 'M612', 'M613', 'M614', 'M615', 'M712', 'M713', 'M714', 'M715', 'C11', 'C22', 'C33', 'C44', 'C55', 'C66', 'C77', 'C88', 'C99', 'C1010', 'C1111', 'C1212', 'C1313', 'C1414', 'C1515', 'Lh', 'Lfa', 'Lua'};  

