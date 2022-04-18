%% ---- Swing up task of 2 dof arm using iLQR or DDP ----

clear; close all; clc;

%% ---compute symbolic  dynamics---
symbolic_2_link_dynamics



%% --- initial state + optimization parameters ---
X0s = {{[-pi/6;0;0;0],200},{[-pi;0;0;0],250},{[-pi/2;0;0;0],450},{[-pi/2;0;0;0],800}};

% --- uncomment for different initial states and time horizonts --- 
x0=X0s{1}{1}; N=X0s{1}{2};
% x0=X0s{2}{1}; N=X0s{2}{2};  
% x0=X0s{3}{1}; N=X0s{3}{2};  
% x0=X0s{4}{1}; N=X0s{4}{2};  
%-------------------------------------------

lambda_init=1000; %initial regularization
max_iter=100; % max number of iterations
lambda_min=1e-6;
lambda_max=1e4;
tolVal=1e-6; % cost tolerance 
alpha=10.^linspace(0,-3,16); %line search vector
lambda_factor=4;  % growth/reduction of lambda
n=4;m=2; % size of state and input variables
tolGrad=1e-4; % gradient tolerance
FullDDP=false; % if true full DDP with second order dynamics is calculated
SimulateZeroTorque = false; %simulate zero input torque 

params=struct('N',N,'x0',x0,'max_iter',max_iter,'lambda_init',lambda_init,'lambda_min',lambda_min,'lambda_max',lambda_max,'lambda_factor',lambda_factor,'n',n,'m',m,'tolVal',tolVal,'tolGrad',tolGrad,'alpha',alpha,'FullDDP',FullDDP,'SimulateZeroTorque',SimulateZeroTorque);

methods = {'numerical','symbolic'}; 
method = methods{1};

%% --- choose between numerical or symbolic differentiation and run ilqr/ddp --- 

if strcmp(method,'numerical')
    DYNCST = numerical_derivatives_cost_dynamics(f_symbolic,FullDDP); 
    [X_f,U_f,logs]=ilqr_ddp_numerical(DYNCST,params);
else
    symbolic_2_link_cost;
    [X_f,U_f,logs]=ilqr_ddp_symbolic(f_symbolic,cost_symbolic,params);
end

%% --- movement visualization --- 


Graph(X_f,dt,model_params,method)

%% --- visualization function --- 
  
  
function Graph(X_r,dt,model_params,method)


L1=model_params(3); L2=model_params(4);
if strcmp(method,'numerical')
    X_r = X_r';
end

f1=figure;
ax=axes(f1);ax.XLim=[-2 2];ax.YLim=[-2 2];ax.XGrid='ON';ax.YGrid='ON';
Link1 = line(ax,'XData',[0, L1*cos(X_r(1,1))],'YData',[0, L1*sin(X_r(1,1))],'lineWidth',4,'color','blue'); 
Link2 = line(ax,'XData',[L1*cos(X_r(1,1)), L2*cos(X_r(1,1) + X_r(1,2)) + L1*cos(X_r(1,1))],'YData',[L1*sin(X_r(1,1)), L2*sin(X_r(1,1) + X_r(1,2)) + L1*sin(X_r(1,1))], 'lineWidth',5,'color','blue'); 
hs=[Link1;Link2];
pause(4*dt)

    for i=2:size(X_r,1)
        
    set(hs,{'XData';'Ydata'},{[0, L1*cos(X_r(i,1))], [0, L1*sin(X_r(i,1))];[L1*cos(X_r(i,1)), L2*cos(X_r(i,1) + X_r(i,2)) + L1*cos(X_r(i,1))],[L1*sin(X_r(i,1)), L2*sin(X_r(i,1) + X_r(i,2)) + L1*sin(X_r(i,1))]})   

    pause(dt)

    end
    delete(Link1)
    delete(Link2)
    delete(f1)
end


    