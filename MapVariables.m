function x = MapVariables(x)
global xupper xdelt xlower;
% map integer variables to a discrete set

% The possible values for x(1) and x(2)
allX1 = xlower:xdelt:xupper;

% Map x(1) and x(2) from the integer values used by GA to the discrete values required.
x(1) = allX1(x(1));
x(2) = allX1(x(2));

