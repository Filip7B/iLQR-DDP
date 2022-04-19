% 


syms X [4 1] real
syms U [2 1] real

n=4; %size of state variables
m=2; %size of control variable


m1=5.5; m2 = 4.5; % kg
L1 = 0.9; L2 = 1;   % meter
Lc1 = 0.5; Lc2 = 0.5; % COM distance (meter)
I1 = (1/3)*(m1)*(L1^2); % Inertia
I2 = (1/3)*(m2)*(L2^2); % Inertia
g = 9.18; % Gravity Acceleration
dt=0.005; % time step

model_params= [m1,m2,L1,L2,Lc1,Lc2,I1,I2,g];

f_symbolic=Symbolic_Discrete_Dynamics_Euler(X,U,dt,model_params); % symbolic eq for next x(i+1)


clear m1 m2 L1 L2 Lc1 Lc2 I1 I2 g
clear X1 X2 X3 X4 X U1 U2 U


function X_dot=Symbolic_Forward_Dynamics_Arm(X,U,model_params)
% x=[ th1 th2 th1_dot th2_dot ]
% unpack the parameters
m1=model_params(1);m2=model_params(2); L1=model_params(3); L2=model_params(4);
Lc1=model_params(5);Lc2=model_params(6); I1=model_params(7); I2=model_params(8); g=model_params(9);
%f1=model_params(10); f2=model_params(11);

%Maxx matrix
syms M [2 2]
M(1,1)=I1+I2+m2*L1^2+2*m2*L1*Lc2*cos(X(2));
M(1,2)=I2+m2*L1*Lc2*cos(X(2));
M(2,1)=I2+m2*L1*Lc2*cos(X(2));
M(2,2)=I2;
% Coriolis
syms C [2 2]
C(1,1)=-2*m2*L1*Lc2*sin(X(2))*X(4);
C(1,2)=-m2*L1*Lc2*sin(X(2))*X(4);
C(2,1)=m2*L1*Lc2*sin(X(2))*X(3);
C(2,2)=0;
% gravity part
tau_g=[-m1*g*Lc1*cos(X(1))-m2*g*(L1*cos(X(1))+Lc2*cos(X(1)+X(2))) ; -m2*g*Lc2*cos(X(1)+X(2))];
% Input torque
B=[1 0; 0 1];
% continous derivative
X_dot=simplify([X(3);X(4);M\(tau_g+B*U-C*X(3:4))]);
end


function X_next= Symbolic_Discrete_Dynamics_Euler(X,U,dt,model_params)
X_next=X+dt*Symbolic_Forward_Dynamics_Arm(X,U,model_params);
end