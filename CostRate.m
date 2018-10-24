function C = CostRate(x)
%  This function is the objective function to calculate the 'cost' metric.
%  (CostRate(x) calculates the cost rate based on the predicted degradation
%  path)

%  Inputs:
%         A vector x contains health index values, x(1) and x(2) map the
%         discrete variates 
%
%  Outputs:
%        The optimized result of the 'cost' metric for optimizing the two threshold.  
%         
%
%  NU_Version1.0
%  by Anqi He, Northeastern University, Boston

x = MapVariables(x);

% Calculate the cost rate
%% initial parameters
global X L Xf C1 C2 C3 Cws Cwc C0 N M delta_t;

%% Create variables and constraints
Xs = x(1);
Xm = x(2); 

%% Generate Objective Function
CR1 = [];CR2 = [];
P1_save = [];P2_save = [];P3_save = [];P4_save = [];P5_save = [];
n1_save = [];n2_save = [];n3_save = [];n4_save = [];n5_save = [];
Lm_save = [];
Lf_save = [];
% %% Calculate the probability of each type  when the current time i=1 
n_1 = 0;n_2 = 0;n_3 = 0;n_4 = 0;n_5 = 0; 
L_m = 0;
L_f = 0;
 for j = 1:M
        if X(1,j)>=Xs && X(1+L,j)<Xm
           n_1 = n_1+1; % Number of paths belong to type 1 at time point 1
           L_m = L_m+(find(X(1+L:end,j)>=Xm,1)-2); % sum of the waiting time of supplier of all paths in type 1
        elseif X(1,j)>=Xs && X(1+L,j)>=Xm && X(1+L,j)<Xf
            n_2 = n_2+1;  % Number of paths belong to type 2 at time point 1
        elseif X(1,j)>=Xs && X(1+L,j)>=Xf
            n_3 = n_3+1;   % Number of paths belong to type 3 at time point 1
            L_f = L_f+(find(X(1+L:-1:1,j)<Xf,1)-2); % sum of the waiting time of customers of all paths in type 3
        else n_4=n_4+1;   % Number of paths belong to type 4 at time point 1
        end
 end
%  P1_0 = n_1/M;P2_0 = n_2/M;P3_0 = n_3/M;P4_0 = n_4/M;   % the probability of each type
%%  Calculate the probability of each type  when the current time i>1
%%%% when i=6 
for i = 2:N-L
    n1 = 0; n2 = 0; n3 = 0; n4 = 0; n5 = 0; 
    Lm = 0;
    Lf = 0;
    for j = 1:M
        if X(i-1,j)<Xs && X(i,j)>=Xs && X(i+L,j)<Xm
           n1 = n1+1;
           Lm = Lm+(find(X(i+L:end,j)>=Xm,1)-2);   %%%% unsolved/degradation trajectory should be better monotone 
        elseif X(i-1,j)<Xs && X(i,j)>=Xs && X(i+L,j)>=Xm&&X(i+L,j)<Xf
            n2 = n2+1;
        elseif X(i-1,j)<Xs && X(i,j)>=Xs && X(i+L,j)>=Xf 
            n3 = n3+1;
            Lf = Lf+(find(X(i+L:-1:1,j)<Xf,1)-2);    %%%% unsolved/degradation trajectory should be better monotone 
        elseif X(i-1,j)<Xs && X(i,j)<Xs
            n4 = n4+1;
        else n5 = n5+1;
        end
    end  
    n1_save = [n1_save,n1];
    n2_save = [n2_save,n2];
    n3_save = [n3_save,n3];
    n4_save = [n4_save,n4];
    n5_save = [n5_save,n5];
    % Save the probability of each type in all time points;save Lm and Lf.
    Lm_save = [Lm_save,Lm];
    Lf_save = [Lf_save,Lf];
end
%% Save the probability of each type in all time points
n1_save = [n_1,n1_save];
n2_save = [n_2,n2_save];
n3_save = [n_3,n3_save];
n4_save = [n_4,n4_save];
n5_save = [0,n5_save];
P1_save = n1_save/M;  %%%%%%%%%%????????
P2_save = n2_save/M;
P3_save = n3_save/M;
P4_save = n4_save/M;
P5_save = n5_save/M;
Lm_save = [L_m,Lm_save];
Lf_save = [L_f,Lf_save];
%% Calculate the value of objective function
%% calculate CR1
 CR1 = sum(P1_save*C1+P2_save*C2+P3_save*C3); 
%% Calculate CR2
 a=1:N-L;
 a=(a+L*ones(1,N-L))*delta_t;
 CR2 = sum((P1_save+P2_save+P3_save).*a);
%% Calculate CR3 and CR4
% waiting time of suppliar(type 1)
% index1 = find(n1_save(1,:)==0);
% for i=1:length(index1)
%     llm=Lm_save(Lm_save~=Lm_save(index1(1,i)));
% end
% nn1=n1_save(n1_save~=0);
% pp1=P1_save(P1_save~=0);
% E_Tws=llm./nn1.*pp1;
index1 = find(n1_save(1,:)==0);
  llm=Lm_save;
for i=length(index1):-1:1
   llm(index1(1,i))=[];
end                      
nn1=n1_save(n1_save~=0);
pp1=P1_save(P1_save~=0);
E_Tws=llm./nn1.*pp1;   %%%%%%%%%%%%% problem 
% waiting time of customer(type 3)
index2 = find(n3_save(1,:)==0);
  llf=Lf_save;
for i=length(index2):-1:1
   llf(index2(1,i))=[];
end                      
nn3=n3_save(n3_save~=0);
pp3=P3_save(P3_save~=0);
E_Twc=llf./nn3.*pp3;  
 %% Calculate CR:judge if E_Tws or E_Twc=[]
 if isempty(E_Twc)
     CR3 = sum(E_Tws)*Cws;
     CR4 = 0;
     CR5 = sum(E_Tws);
     CR=C0+(CR1+CR3+CR4)/(CR2+CR5);
 elseif isempty(E_Tws) 
     CR3 = 0;
     CR4 = sum(E_Twc)*Cwc;
     CR5 = sum(-E_Twc);
     CR=C0+(CR1+CR3+CR4)/(CR2+CR5);
 else
     CR3 = sum(E_Tws)*Cws;
     CR4 = sum(E_Twc)*Cwc;
     CR5 = sum(E_Tws)-sum(E_Twc);
     CR=C0+(CR1+CR3+CR4)/(CR2+CR5);
 end    
 C = CR;
