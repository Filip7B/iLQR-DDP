%--------------- LICENSE-------------------------------
% This code is based on:
% https://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization
% Copyright (c) 2015, Yuval Tassa
% All rights reserved.

% The theoretical background of the algorithm is in:
% BIBTeX:
% @INPROCEEDINGS{
% author={Tassa, Y. and Mansard, N. and Todorov, E.},
% booktitle={Robotics and Automation (ICRA), 2014 IEEE International Conference on},
% title={Control-Limited Differential Dynamic Programming},
% year={2014}, month={May}, doi={10.1109/ICRA.2014.6907001}}
%-------------------------------------------------------------


%%% ---DESCRIPTION----
% This function returns DYNCST which is a function that bundles together: cost, forward
% stepping(dynamics) and the derivatives depending on the number of outputs
% (nargout). Dynamics and cost function derivatives are calculated with finite
% differences. The cost function structure and its parameters are the same
% as in symbolic_2_link_cost.m
%--------------------


function DYNCST = numerical_derivatives_cost_dynamics(f_symbolic,FullDDP)

  f=matlabFunction(f_symbolic);
  dynamics=@(x,u)f_next(f,x,u);
  DYNCST=@(x,u) dyn_cst(x,u,dynamics,FullDDP);

end


function c=cost(x,u)
    
    %--- Xf_d=[pi/2;0;0;0];
    final = isnan(u(1,:)); % creating a sequence of booleans [0 0 0 0 0...1]
    u(:,final)  = 0; % last rows that were nans become 0

    lu=0.5*(0.1*u(1,:).^2+0.1*u(2,:).^2);

    if any(final)
        llf      = 0.5*[1e6 1e6 1e3 1e3]*(x(:,final)-[pi/2;0;0;0]).^2;
        lf       = double(final); 
        lf(final)= llf; % assigns the last element the llf
    else
        lf    = 0;
    end

    % total cost
    c     = lu  + lf;

end     

function x_next=f_next(f,x,u)
    x_next=f(u(1,:),u(2,:),x(1,:),x(2,:),x(3,:),x(4,:));
end



function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = dyn_cst(x,u,dynamics,full_DDP)
% combine dynamics and cost
% use helper function finite_difference() to compute derivatives

    if nargout == 2 % computes next state and cost (when x and u are matrices 
        f = dynamics(x,u);
   
        c = cost(x,u);
    else
    % state and control indices
        ix = 1:4;
        iu = 5:6;
    
    % dynamics first derivatives
        xu_dyn  = @(xu) dynamics(xu(ix,:),xu(iu,:));
        J       = finite_difference(xu_dyn, [x; u]);
        fx      = J(:,ix,:); % size is n x (n+m) x (N+1)  for boths fx and fu where  (N+1) matrix is nan
        fu      = J(:,iu,:);
    
        % dynamics second derivatives
        if full_DDP
            xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
            JJ      = finite_difference(xu_Jcst, [x; u]);
            JJ      = reshape(JJ, [4 6 size(J,2,3)]); % fixed, previously size(J)
            JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4])); %symmetrize.... size is [4 6 6 N+1] where N+1 tensor is nan
            fxx     = JJ(:,ix,ix,:); % [4 4 4 N+1]
            fxu     = JJ(:,ix,iu,:); % [4 4 2 N+1]
            fuu     = JJ(:,iu,iu,:); % [4 2 2 N+1]
        else
            [fxx,fxu,fuu] = deal([]);
        end    
    
        % cost first derivatives
        xu_cost = @(xu) cost(xu(ix,:),xu(iu,:));
        J       = squeeze(finite_difference(xu_cost, [x; u]));
        cx      = J(ix,:); % [4, N+1] note these gradients are tall
        cu      = J(iu,:);  %[2,N+1]  N+1 is 0 because u doesnt affect the last cost
    
        % cost second derivatives
        xu_Jcst = @(xu) squeeze(finite_difference(xu_cost, xu));
        JJ      = finite_difference(xu_Jcst, [x; u]);
        JJ      = 0.5*(JJ + permute(JJ,[2 1 3])); %symmetrize ... of size (m+n)x(m+n)*(N+1)
        cxx     = JJ(ix,ix,:); 
        cxu     = JJ(ix,iu,:);
        cuu     = JJ(iu,iu,:);
    
        [f,c] = deal([]);
    end

end


function J = finite_difference(fun, x, h)
    % simple finite-difference derivatives
    % assumes the function fun() is vectorized

    if nargin < 3
        h = 2^-17;

    
    end

    [n, K]  = size(x);
    H       = [zeros(n,1) h*eye(n)];
    H       = permute(H, [1 3 2]);
    X        =x+H;
    X       = reshape(X, n, K*(n+1));
    Y       = fun(X);

    m       = numel(Y)/(K*(n+1));
    Y       = reshape(Y, m, K, n+1);
    J       = (Y(:,:,2:end)-Y(:,:,1)) / h;
    J       = permute(J, [1 3 2]);

end
