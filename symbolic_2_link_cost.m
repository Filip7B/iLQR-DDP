% --- DESCRIPTION---
% This script generates running and final cost in the symbolic form. 
% The running cost (the cost at each time step up) is only the function of
% input torques U. Final cost penalizes deviations from the final state:
% [pi/2 0 0 0]. This script is used only if the symbolic differentiation
% method is selected.
% ------------------

syms X [4 1] real
syms U [2 1] real



    
R=1e-1*eye(2,2);
Qf=1000*eye(4);
Qf(1,1)=1e3*Qf(1,1);
Qf(2,2)=Qf(1,1);

l_run=(1/2)*U'*R*U;
l_fin=(1/2)*(X-[pi/2;0;0;0])'*Qf*(X-[pi/2;0;0;0]);

cost_symbolic={l_run,l_fin};

clear X1 X2 X3 X4 X U1 U2 U Qf R l_run l_fin