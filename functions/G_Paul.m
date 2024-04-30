function [T] = G_Paul(I, D, K, M, t1, t2, C, w2)
% G Summary of this function goes here
% Detailed explanation goes here
% identity matrix was 1 in original code of Paul Curtis

T = inv(-I.*w2.^2 + D*(1i.*w2) + K)*(M*inv(-(t1*t2).*w2.^2 + (t1+t2)*(1i.*w2) + eye(15))*C);
end

