function [c, ceq] = constraint(x)
% nonlinear constraints
c = x(1)-x(2);  

% No equality constraints
ceq = [];