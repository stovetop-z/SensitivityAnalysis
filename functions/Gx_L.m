function [X] = Gx_L(Js, L_p, I_p, D_p, K_p, M_p, t1_p, t2_p, C_p, input, w)
syms Lh Lfa Lua

js = double(subs(Js, [Lh Lfa Lua], L_p));
X = (js*inv(-I_p.*(w^2) + D_p.*1i.*w + K_p)*M_p*inv(-((t1_p*t2_p).*(w^2)) + ((t1_p+t2_p).*1i.*w) + eye(15))*C_p) .* input;
end

