function [X] = Gx(J, L, I, D, K, M, t1, t2, C, w)
% Set jacobian
syms Lh Lfa Lua

j = subs(J, [Lh Lfa Lua], L);
j = double(j);

X = j*inv(-I.*(w^2) + D.*1i.*w + K)*M*inv(-((t1*t2).*(w^2)) + ((t1+t2).*1i.*w) + eye(15))*C;
end