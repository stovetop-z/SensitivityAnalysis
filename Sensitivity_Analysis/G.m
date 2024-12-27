function [T] = G(G3inv, M, C, t1, t2, w)
% G Summary of this function goes here
% Detailed explanation goes here
% identity matrix was 1 in original code of Paul Curtis

T = G3inv*M*inv(-(t1*t2).*w.^2 + (t1+t2).*1i.*w + eye(15))*C;
end
