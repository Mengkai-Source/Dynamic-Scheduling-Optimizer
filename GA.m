%%%%%----- Genetic Algorithm for discrete optimization (DMDII)----%%%%%
%%%%%----- Anqi He (20180722) %%%%%
%%%%% Main Program %%%%%
close all
clear all
clc
tic
% Global 
global X L Xf C1 C2 C3 Cws Cwc N M C0 delta_t xlower xupper xdelt;
% Import the data of the future predicted degradation path
load('path.mat');
X = MS;
%% feasible region
Xf = 1;  % Provided by PDX
L = 50;  % leading time(settle down regarding the real case )
xdelt= 0.005;
% Pre-check the validity of the predicted paths  
% Remove the columns in which the last obsevation is less than the value of failure threshold
if find(X(end,:)<Xf)==[]
    X = X;
else
    X(:,find(X(end,:)<Xf))=[];
end 
[row,col] = size(X);
N = row; % The time length of predicted degradation paths( provided by PDX )
M = col; % Number of predicted paths

% Plot the data of predicted paths
figure(1)
for i= 1:M
    plot(X(:,i));hold on;
end

% Initial Parameters
C1 = 15; % Maintenance cost for type 1
C2 = 20; % Maintenance cost for type 2
C3 = 40; % Maintenance cost for type 3
Cws = 1; % Cost of waiting time of suppliers per unit time  %%% original 1
Cwc = 1000; % Cost of waiting time of suppliers per unit time 10 %%% original 1000
C0 = 0.5; % Initial cost can be removed
delta_t = 1 ; % Unit time (provided by PDX)
%% GA Algorithm
%% Create vectors containing the lower bound (|lb|) and upper bound constraints(|ub|).
% Transform the bounds on the discrete variables.
xlower = X(1,1); %%% 
xupper = Xf;
% [0,3] deltaX = 0.01
lb = [1 1];   
ub = [length(xlower:xdelt:xupper) length(xlower:xdelt:xupper)];
% Now we can call |ga| to solve the problem with discrete variables.
rng(0, 'twister');
[xbestDisc, fbestDisc, exitflagDisc] = ga(@CostRate, ...
    2, [], [], [], [], lb, ub, @constraint, 1:2);  %%% need to be modified
toc
% _Analyze the Results_
xbestDisc = MapVariables(xbestDisc);
display(xbestDisc)
fprintf('\nCost rate return by ga = %g\n',fbestDisc);

