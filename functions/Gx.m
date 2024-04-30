function [X] = Gx(J, L, I, D, K, M, t1, t2, C, w2)
% Set jacobian
syms Lh Lfa Lua

j = subs(J, [Lh Lfa Lua], L);
j = double(j);

X = j*inv(-I.*w2.^2 + D.*1i.*w2 + K)*M*inv(-(t1*t2).*w2.^2 + (t1+t2).*1i.*w2 + eye(15))*C;
end