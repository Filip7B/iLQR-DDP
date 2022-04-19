function [X,U,logs]=ilqr_ddp_symbolic(f_symbolic,cost_symbolic,par)

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


%-----------------DESCRIPTION-------------------
% This code is a reduced version of the code cited above since it excludes
% control limits,  graphing functions, serial line search
% The key difference is that this version uses symbolic differentiation to create derivative functions.
% The inputs to the function are symbolic dynamics, symbolic cost and
% optimization parameters.
% The outputs are optimized state and control trajectories X,U and the logs
% file for storing optimization updates
%------------------------------------------------

l_r_symbolic=cost_symbolic{1}; l_f_symbolic=cost_symbolic{2};


N=par.N; x0=par.x0; max_iter=par.max_iter; lambda_init=par.lambda_init; lambda_min=par.lambda_min; lambda_max=par.lambda_max; lambda_factor=par.lambda_factor; n=par.n; m=par.m; tolVal=par.tolVal; tolGrad=par.tolGrad; alpha=par.alpha; FullDDP=par.FullDDP;
pass_params={n,m,N};





%% ---initialize logs and print structure ----
logs = struct('iter',nan,'lambda',nan,'cost',nan,...
        'alpha',nan,'grad_norm',nan,'improvement',nan,...
        'time_derivs',nan,'time_forward',nan,'time_backward',nan,'time_symbolic_diff',nan,'diverge_backward',nan);
logs= repmat(logs,[max_iter 1]);

print_head  = 6; % print headings every print_head lines
last_head   = print_head;



%% ---- get derivatives and dynamics as functions from their symbolic expressions----

t_symbolic_diff=tic;

f = matlabFunction(f_symbolic);
l_r = matlabFunction(l_r_symbolic);
l_f = matlabFunction(l_f_symbolic);
pass_f_handles = {f,l_r,l_f};

[fx,fu,lx,lxx,lu,luu,lux,lfx,lfxx,fxx,fux,fuu]=derivatives_fast(f_symbolic,cost_symbolic{1},cost_symbolic{2});
 

%% --- initialize constant part of the derivatives --- 
% The reason why derivative functions are split into constant and variable
% parts is to utalize vectorization during optimization. 
% This part of the code is tailored specifically to the structure of the dynamics and
% the cost function. If different dynamical system and/or cost function is
% used, this part should be changed.


fx_values = zeros(16,N-1);
fx_values(fx.fx_const_ind,:) = repmat(fx.fx_const,1,N-1);


fu_values = zeros(8,N-1);
fu_values(fu.fu_const_ind,:) = repmat(fu.fu_const,1,N-1);

fxx_values = zeros(64,N-1);
fxx_values(fxx.fxx_const_ind,:)=repmat(fxx.fxx_const,1,N-1); % assigning constnt values to fxx_values


fux_values = zeros(32,N-1);
fux_values(fux.fux_const_ind,:)=repmat(fux.fux_const,1,N-1); % assigning constnt values to fux_values

fuu_values = zeros(16,N-1);
fuu_values(fuu.fuu_const_ind,:)=repmat(fuu.fuu_const,1,N-1); % assigning constnt values to fuu_values

% ---zero at all times---
lxs=lx(:,:,ones(1,N-1));
lxxs=lxx(:,:,ones(1,N-1));
luxs=lux(:,:,ones(N-1));
luus = luu(:,:,ones(1,N-1));
lfxxs = lfxx;


logs(1).time_symbolic_diff=toc(t_symbolic_diff);
% -----------------------------------------


%% --- initialize  values ---
% U=zeros(N-1,m);
U=zeros(N-1,m);
X=Rollout(f,x0,U);
cost=total_cost(U,X(end,:),l_r,l_f);
lambda=lambda_init;
if par.SimulateZeroTorque == true
    fprintf('\n No optimization! Simulating dynamics with zero torque.\n')
    return
end




%% ------begin ddp-----

fprintf('\n=========== begin iLQG ===========\n');

for j=1:max_iter
    
    

    

  %% --- calculate variable part of the derivative around current control and state trajectory ---
    
    logs(j).iter=j;
    time_deriv=tic;

    fx_values(fx.fx_var_ind,:) = fx.fx_f(U(1:N-1,1)',U(1:N-1,2)',X(1:N-1,1)',X(1:N-1,2)',X(1:N-1,3)',X(1:N-1,4)');
    fxs = reshape(fx_values,[4,4,N-1]);
  
    fu_values(fu.fu_var_ind,:) = fu.fu_f(X(1:end-1,2)');
    fus = reshape(fu_values,[4,2,N-1]);
  
    lus = lu(U(:,1)',U(:,2)');
    lus = reshape(lus,[2,1,N-1]);
   
    lfxs=lfx(X(end,1)',X(end,2)',X(end,3)',X(end,4)'); 
  
    if FullDDP
      
        fxx_values(fxx.fxx_var_ind,:) = fxx.fxx_v(U(1:N-1,1)',U(1:N-1,2)',X(1:N-1,1)',X(1:N-1,2)',X(1:N-1,3)',X(1:N-1,4)');
        fxxs = permute(reshape(fxx_values,[4,4,4,N-1]),[3 1 2 4]);
  
        fux_values(fux.fux_var_ind,:) = fux.fux_v(X(1:N-1,2)');
        fxus= permute(reshape(fux_values,[4,2,4,N-1]),[3 1 2 4]);
  

        fuus = permute(reshape(fuu_values,[2,2,4,N-1]),[3 1 2 4]); %fuu is constant and zero everywhere
  
    else
      
        fxxs = [];
        fxus = [];
        fuus = [];
  
    end
   
    logs(j).time_derivs=toc(time_deriv);
    
  %% --- backward  pass ---
        backPassDone = 0;
   
    while ~backPassDone 
        t_backward=tic;
        [diverge,K,k,dV] = backward(fxs,fus,lxs,lus,lxxs,luus,luxs,lfxs,lfxxs,fxxs,fuus,fxus,lambda,pass_params);
        logs(j).time_backward=toc(t_backward);
        if diverge
           logs(j).diverge_backward=diverge;
           fprintf('Cholesky failed at timestep %d during iteration %d!\n',diverge,j);
           lambda=max(lambda*lambda_factor,lambda_min);
           if lambda>lambda_max
               fprintf('Backward pass failed. Lambda reached maximum, the algorithm terminates!');
               return;
           end
           continue;
        end
        backPassDone = 1;                         
    end
    
    g_norm = mean(max(abs(k) ./ (abs(U')+1),[],1)); % =[ max(|d1|./|u1|) max(|du2|/|u2|) ......max(|dun|/|un|)] -> it checks the magnitude of the gradient norm when this is small we can be sure that every update will be very small
    logs(j).grad_norm = g_norm;
    if g_norm < tolGrad && lambda < 1e-5
        fprintf(' \n Success: gradient norm is smaller than tolGrad! \n');
        break
    end
    
    
    %% ---forward pass--- 
     t_forward=tic;
     forwardPassDone=0;
     [Xnew,Unew,Costnew]=forward(X,U,K,k,alpha,pass_params,pass_f_handles);
     Dcost=cost*ones(length(alpha),1)-Costnew;
     [dcost, w]=max(Dcost);
     alpha_best=alpha(w);
     expected= -alpha_best*(dV(1)+alpha_best*dV(2));
    
     logs(j).time_forward=toc(t_forward);
     
     if expected > 0
         z = dcost/expected;
     else
         z=sign(dcost);
         warning('Negative expected change of cost should not occur!');
     end
     
     if z > 0
        forwardPassDone=1;
        X=Xnew(:,:,w);
        U=Unew(:,:,w);
        cost=Costnew(w);
        logs(j).alpha = alpha(w);
     end 
     logs(j).time_forward=toc(t_forward);
     
     
     
  
     %% --- accept new trajectories (or not), change lambda and print ---
     
     logs(j).lambda=lambda;
     logs(j).cost=cost; %in case we didn't lower the cost write the old cost in the logs struct
     logs(j).grad_norm=g_norm;
     logs(j).improvement=dcost;
     
     if  last_head == print_head
        last_head = 0;
        fprintf('%-12s','iteration','cost','reduction','expected','gradient','lambda')
        fprintf('\n');
     end
     
     if forwardPassDone
          fprintf('%-12d%-12.6g%-12.3g%-12.3g%-12.3g%-12.1f\n',j, cost, dcost, expected, g_norm, lambda);
          last_head = last_head+1;
          
          lambda = lambda*(lambda > lambda_min)/lambda_factor;
          
          if dcost<tolVal
              fprintf('\nSuccess! Cost change < tolVal');
              break
          end
          
          
     else
              
         fprintf('%-12d%-12s%-12.3g%-12.3g%-12.3g%-12.1f\n',j,'NO STEP', dcost, expected, g_norm, lambda);           
         last_head = last_head+1;
         
         lambda=max(lambda*lambda_factor,lambda_min);
         
         
         if lambda > lambda_max
             fprintf('\nForward pass terminated without cost improvement and Lambda > LambdaMax!\n');
             break
         end
     end 
     
 end

    
    
    
    
    
    
    


%% --- functions ---

function [xnew,unew,cnew]=forward(X,U,K,k,alpha,pass_params,pass_f_handles)
    n=pass_params{1};m=pass_params{2};N=pass_params{3};
    
    
    alphal= length(alpha);
    dx=zeros(alphal,n);
    Kl=ones(alphal,1); % useful for expansion
    xnew        = zeros(alphal,n,N); % third dimension is time
    xnew(:,:,1) = X(Kl,:); % copy all the x0s len alpha times in xnew
    unew        = zeros(alphal,m,N-1);
    cnew        = zeros(1,alphal,N);
    
   
    for i=1:N-1
        unew(:,:,i)=U(i*Kl,:)+(k(:,i)*alpha)'+dx*K(:,:,i)';
       
        [cnew(:,:,i),xnew(:,:,i+1)]=Step_and_Cost(xnew(:,:,i),unew(:,:,i),pass_f_handles);
        dx=xnew(:,:,i+1)-X((i+1)*Kl,:);
    end

    cnew(:,:,N)=Step_and_Cost(xnew(:,:,N),nan(alphal,m,1),pass_f_handles);
    
    xnew = permute(xnew, [3 2 1]);
    unew = permute(unew, [3 2 1]);
    cnew = permute(cnew, [1 3 2]);
    cnew=sum(cnew,2);
    cnew=reshape(cnew,[size(cnew,3),1]);


    
    function [cost,xnext]=Step_and_Cost(x,u,pass_f_handles) % here watch out about the order of output variables when using nargout: if xnext is not returned in certain case then put it second!!
        f=pass_f_handles{1};l_r=pass_f_handles{2};l_f=pass_f_handles{3};
        if nargout==2
            
            xnext=(f(u(:,1)',u(:,2)',x(:,1)',x(:,2)',x(:,3)',x(:,4)'))';
            cost=l_r(u(:,1)',u(:,2)');
        else
            
            cost=l_f(x(:,1)',x(:,2)',x(:,3)',x(:,4)');
        end
            


function [diverge,K,k,dV]=backward(fxs,fus,lxs,lus,lxxs,luus,luxs,lfxs,lfxxs,fxx,fuu,fxu,regu,pass_params)
    n=pass_params{1};m=pass_params{2};N=pass_params{3};
    Vx=zeros(n,N);
    Vxx=zeros(n,n,N);
    dV=[0 0];
    k=zeros(m,N-1);
    K=zeros(m,n,N-1);
    diverge=0;

    %%-- Initialize V_N--
    Vx(:,N)=lfxs;
    Vxx(:,:,N)=lfxxs;

        for p=N-1:-1:1
 
            Qx=lxs(:,:,p)+fxs(:,:,p)'*Vx(:,p+1);
            Qu=lus(:,:,p)+fus(:,:,p)'*Vx(:,p+1);
        
        
            Qxx=lxxs(:,:,p)+fxs(:,:,p)'*Vxx(:,:,p+1)*fxs(:,:,p);
        
             if ~isempty(fxx)
                Qxx = Qxx + permute(sum(Vx(:,p+1).*fxx(:,:,:,p),1),[3 2 1]);
             end
        
            Qux=luxs(:,:,p)+fus(:,:,p)'*Vxx(:,:,p+1)*fxs(:,:,p);
        
            if ~isempty(fxu)
                Qux   = Qux + permute(sum(Vx(:,p+1).*fxu(:,:,:,p),1),[3 2 1]);
            end
        
        
            Quu=luus(:,:,p)+fus(:,:,p)'*Vxx(:,:,p+1)*fus(:,:,p);
        
            if ~isempty(fuu)
                Quu   = Quu + permute(sum(Vx(:,p+1).*fuu(:,:,:,p),1),[3 2 1]);
            end
        
        
            Quu_reg=Quu+eye(m)*regu;
        
            [R,d] = chol(Quu_reg);
            if d ~= 0
                diverge  = p; % at which backward step out of N did inversion fail
                return;
            end
        
        % calculate gains
            kK= -R\(R'\[Qu Qux]);
            k(:,p) = kK(:,1);
            K(:,:,p) = kK(:,2:end);
        
            k_p = k(:,p);
            K_p = K(:,:,p);
                
            Vx(:,p)     = Qx  + K_p'*Quu*k_p + K_p'*Qu  + Qux'*k_p;
            Vxx(:,:,p)  = Qxx + K_p'*Quu*K_p + K_p'*Qux + Qux'*K_p;
        
            Vxx(:,:,p)=0.5*(Vxx(:,:,p)'+Vxx(:,:,p)); % this is symmetrize operation
            dV= dV  + [k_p'*Qu, 0.5*k_p'*Quu*k_p]; 
        
        
        end    

       

function [Fx,Fu,Lx,Lxx,Lu,Luu,Lux,Lfx,Lfxx,fxx,fux,fuu]=derivatives_fast(f_symbolic,l_r_symbolic,l_f_symbolic)

    var_f=symvar(f_symbolic); var_l_r=symvar(l_r_symbolic); var_l_f=symvar(l_f_symbolic); % get symbolic variables


    %% fx ([4 4 N-1]) and fxx ([4 4 4 N-1])
    fx=jacobian(f_symbolic,var_f(3:6)); % calculate jacobian of dynamics function this is 4x4 symbolic matrix
    fx_v = reshape(fx,16,1); % fx vectorized
    fx_var_ind = find(hasSymType(fx_v,'variable'));
    fx_const_ind = find(~hasSymType(fx_v,'variable'));
    fx_const = double(fx_v(fx_const_ind));
    fx_f = matlabFunction(fx_v(fx_var_ind));
    Fx = struct('fx_f',fx_f,'fx_var_ind',fx_var_ind,'fx_const',fx_const,'fx_const_ind',fx_const_ind);


    fxx_1 = jacobian(fx(1,:),var_f(3:6));
    fxx_2 = jacobian(fx(2,:),var_f(3:6));
    fxx_3 = jacobian(fx(3,:),var_f(3:6));
    fxx_4 = jacobian(fx(4,:),var_f(3:6));

    fxx = [reshape(fxx_1,16,1);reshape(fxx_2,16,1);reshape(fxx_3,16,1);reshape(fxx_4,16,1)]; % vectorize fxx

    fxx_var_ind = find(hasSymType(fxx,'variable')); % find indices that correspond to parts of fxx that depend on x1..x4
    fxx_const_ind = find(~hasSymType(fxx,'variable'));% find indices that are constant

    fxx_const = double(fxx(fxx_const_ind)); % extract constant parts of fxx
    fxx_v = matlabFunction(fxx(fxx_var_ind)); % create function from variable parts of fxx only. 

    fxx = struct('fxx_v',fxx_v,'fxx_var_ind',fxx_var_ind,'fxx_const',fxx_const,'fxx_const_ind',fxx_const_ind); 



    %% fux (final form [4 4 2 N-1])


    fux_1 = jacobian(fx(1,:),var_f(1:2)); 
    fux_2 = jacobian(fx(2,:),var_f(1:2));
    fux_3 = jacobian(fx(3,:),var_f(1:2));
    fux_4 = jacobian(fx(4,:),var_f(1:2));

    fux = [reshape(fux_1,8,1);reshape(fux_2,8,1);reshape(fux_3,8,1);reshape(fux_4,8,1)];

    fux_var_ind = find(hasSymType(fux,'variable'));
    fux_const_ind = find(~hasSymType(fux,'variable'));

    fux_const = double(fux(fux_const_ind));
    fux_v = matlabFunction(fux(fux_var_ind));

    fux = struct('fux_v',fux_v,'fux_var_ind',fux_var_ind,'fux_const',fux_const,'fux_const_ind',fux_const_ind);




    %% fu (final form [4 2 N-1]) and fuu (final form [4 2 2 N-1])

    fu=jacobian(f_symbolic,var_f(1:2));


    fu_v = reshape(fu,8,1); 
    fu_var_ind = find(hasSymType(fu_v,'variable'));
    fu_const_ind = find(~hasSymType(fu_v,'variable'));
    fu_const = double(fu_v(fu_const_ind));
    fu_f = matlabFunction(fu_v(fu_var_ind));
    Fu = struct('fu_f',fu_f,'fu_var_ind',fu_var_ind,'fu_const',fu_const,'fu_const_ind',fu_const_ind);




    fuu_1 = jacobian(fu(1,:),var_f(1:2));
    fuu_2 = jacobian(fu(2,:),var_f(1:2));
    fuu_3 = jacobian(fu(3,:),var_f(1:2));
    fuu_4 = jacobian(fu(4,:),var_f(1:2));

    fuu = [reshape(fuu_1,4,1);reshape(fuu_2,4,1);reshape(fuu_3,4,1);reshape(fuu_4,4,1)];

    fuu_var_ind = find(hasSymType(fuu,'variable'));
    fuu_const_ind = find(~hasSymType(fuu,'variable'));

    fuu_const = double(fuu(fuu_const_ind));
    fuu_v = matlabFunction(fuu(fuu_var_ind));

    fuu = struct('fuu_v',fuu_v,'fuu_var_ind',fuu_var_ind,'fuu_const',fuu_const,'fuu_const_ind',fuu_const_ind);



    %% -- some predictable terms-- found after examining the functions

    lu=gradient(l_r_symbolic,var_l_r); % strictly function of u
    Lu = matlabFunction(lu);

    Luu = double(jacobian(lu,var_l_r)); % constant

    lfx=gradient(l_f_symbolic,var_l_f);
    Lfx=matlabFunction(lfx); % strictly function of x
    Lfxx=double(jacobian(lfx,var_l_f)); % constant

    % zero terms 
    Lx=zeros(4,1); 
    Lxx=zeros(4,4); 
    Lux=zeros(2,4);



function X_r=Rollout(f,x0,U_r)

    X_r=zeros(size(U_r,1)+1,size(x0,1));

    X_r(1,:)=x0';

        for i=1:size(U_r,1)
            X_r(i+1,:)=f(U_r(i,1),U_r(i,2),X_r(i,1),X_r(i,2),X_r(i,3),X_r(i,4));
        end


function C=total_cost(U,Xend,l_r,l_f)
    C=sum(l_r(U(:,1),U(:,2)))+l_f(Xend(1),Xend(2),Xend(3),Xend(4));

