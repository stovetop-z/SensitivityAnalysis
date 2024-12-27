function [T] = Orig_G(I, D, K, M, t1, t2, C, w)
% G Summary of this function goes here
% Detailed explanation goes here
% identity matrix was 1 in original code of Paul Curtis

T = inv(-I.*(w^2) + D.*(1i*w) + K)*M*inv(-((t1*t2).*(w^2)) + ((t1+t2).*(1i*w)) + 1)*C;
end

