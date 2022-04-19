function [X,U,logs]=ilqr_ddp_numerical(DYNCST,par)


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


%-------------------DESCRITPION-----------------------
% This code is a reduced version of the code cited above since it excludes
% control limits,  graphing functions, serial line search.
% The outputs are optimized state and control trajectories X,U and the logs
% file for storing optimization updates
%-----------------------------------------------------

N=par.N; x0=par.x0; max_iter=par.max_iter; lambda_init=par.lambda_init; lambda_min=par.lambda_min; lambda_max=par.lambda_max; lambda_factor=par.lambda_factor; n=par.n; m=par.m; tolVal=par.tolVal; tolGrad=par.tolGrad; alpha=par.alpha; FullDDP=par.FullDDP;



pass_params={n,m,N};




%% ---initialize logs and print structure ----
logs = struct('iter',nan,'lambda',nan,'cost',nan,...
        'alpha',nan,'grad_norm',nan,'improvement',nan,...
        'time_derivs',nan,'time_forward',nan,'time_backward',nan,'time_symbolic_diff',nan,'diverge_backward',nan);
logs= repmat(logs,[max_iter 1]);

print_head  = 6; % print headings every print_head lines
last_head   = print_head;



%% --- initialize values---
U=zeros(m,N-1);
[X,cost]=Rollout(DYNCST,x0,U);
lambda=lambda_init;
if par.SimulateZeroTorque == true
    fprintf('\n No optimization! Simulating dynamics with zero torque.\n')
    return
end


%% ------begin ddp-----

fprintf('\n=========== begin iLQG ===========\n');

for j=1:max_iter
    
    
    logs(j).iter=j;
    time_deriv=tic;    
    
    [~,~,fxs,fus,fxxs,fxus,fuus,lxs,lus,lxxs,lxus,luus]   = DYNCST(X, [U nan(m,1)]);
    
    
    logs(j).time_derivs=toc(time_deriv);
    
    %% --- backward ---
    
   
     t_backward=tic;

    
    backPassDone = 0;
   
    while ~backPassDone 
        t_backward=tic;
        
        lfxs = lxs(:,N);
        lfxxs = lxxs(:,:,N);
        
        [diverge,K,k,dV] = backward(fxs,fus,lxs,lus,lxxs,luus,lxus,lfxs,lfxxs,fxxs,fuus,fxus,lambda,pass_params);
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
    

    logs(j).time_backward=toc(t_backward);
    
    
    g_norm = mean(max(abs(k) ./ (abs(U)+1),[],1)); % =[ max(|d1|./|u1|) max(|du2|/|u2|) ......max(|dun|/|un|)] -> it checks the magnitude of the gradient norm when this is small -> every update will be  small
    logs(j).grad_norm = g_norm;
    if g_norm < tolGrad && lambda < 1e-5
        fprintf(' \n Success: gradient norm is smaller than tolGrad! \n');
        break
    end
        
   
    %% ---forward--- 
     t_forward=tic;
     forwardPassDone=0;
    
     [Xnew,Unew,Costnew]=forward(X,U,K,k,alpha,pass_params,DYNCST);
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

     % update logs
     
     logs(j).lambda=lambda;
     logs(j).cost=cost; % in case we didn't lower the cost  just write the old cost
     logs(j).grad_norm=g_norm;
     logs(j).improvement=dcost;
     
     
     if  last_head == print_head
        last_head = 0;
        fprintf('%-12s','iteration','cost','reduction','expected','gradient','lambda')
        fprintf('\n');
     end
    
     % print status
     if forwardPassDone
         
         fprintf('%-12d%-12.6g%-12.3g%-12.3g%-12.3g%-12.1f\n',j, cost, dcost, expected, g_norm, lambda);

         last_head = last_head+1;
         lambda = lambda*(lambda > lambda_min)/lambda_factor;
         if dcost<tolVal
              fprintf('\nSuccess! Cost change < tolVal');
              break
         end
     
     else
          fprintf('%-12d%-12.6g%-12.3g%-12.3g%-12.3g%-12.1f\n',j, cost, dcost, expected, g_norm, lambda);
          last_head = last_head+1;
          lambda=max(lambda*lambda_factor,lambda_min);

          if lambda > lambda_max
             fprintf('\nForward pass terminated without cost improvement and Lambda > LambdaMax!\n');
             break
         end
     end
     
 end


%% --- functions ---

function [xnew,unew,cnew]=forward(X,U,K,k,alpha,pass_params,DYNCST)
    n=pass_params{1};m=pass_params{2};N=pass_params{3};
    
    
    alphal= length(alpha);
    dx=zeros(n,alphal);
    Kl=ones(1,alphal); % useful for expansion
    xnew        = zeros(n,alphal,N); % third dimension is time
    xnew(:,:,1) = X(:,Kl); % copy all the x0s len alpha times in xnew
    unew        = zeros(m,alphal,N-1);
    cnew        = zeros(1,alphal,N);
    
   
    for i=1:N-1
        unew(:,:,i)=U(:,i*Kl)+k(:,i)*alpha+K(:,:,i)*dx;
       
        [xnew(:,:,i+1),cnew(:,:,i)]=DYNCST(xnew(:,:,i),unew(:,:,i));
        dx=xnew(:,:,i+1)-X(:,(i+1)*Kl);
    end

    [~,cnew(:,:,N)]=DYNCST(xnew(:,:,N),nan(m,alphal));
    
    xnew = permute(xnew, [1 3 2]);
    unew = permute(unew, [1 3 2]);
    cnew=squeeze(cnew);
    cnew=sum(cnew,2);




function [diverge,K,k,dV]=backward(fxs,fus,lxs,lus,lxxs,luus,lxus,lfxs,lfxxs,fxxs,fuus,fxus,lambda,pass_params)
    n=pass_params{1};m=pass_params{2};N=pass_params{3};
    Vx=zeros(n,N);
    Vxx=zeros(n,n,N);
    dV=[0 0];
    k=zeros(m,N-1);
    K=zeros(m,n,N-1);
    diverge=0;

    Vx(:,N) = lfxs;
    Vxx(:,:,N) = lfxxs;

    for p=N-1:-1:1
 
        Qx=lxs(:,p)+fxs(:,:,p)'*Vx(:,p+1);
        Qu=lus(:,p)+fus(:,:,p)'*Vx(:,p+1);
        
        Qxx=lxxs(:,:,p)+fxs(:,:,p)'*Vxx(:,:,p+1)*fxs(:,:,p);
        
        if ~isempty(fxxs)
            Qxx = Qxx + permute(sum(Vx(:,p+1).*fxxs(:,:,:,p),1),[3 2 1]);
        end
        
        Qux=lxus(:,:,p)'+fus(:,:,p)'*Vxx(:,:,p+1)*fxs(:,:,p);
        
        if ~isempty(fxus)
            Qux   = Qux + permute(sum(Vx(:,p+1).*fxus(:,:,:,p),1),[3 2 1]);
        end
        
        
        
        
        Quu=luus(:,:,p)+fus(:,:,p)'*Vxx(:,:,p+1)*fus(:,:,p);
        

        if ~isempty(fuus)   
            Quu   = Quu + permute(sum(Vx(:,p+1).*fuus(:,:,:,p),1),[3 2 1]);
        end
        
        
        Quu_reg=Quu+eye(m)*lambda;
        
        
        
  
        
        [R,d] = chol(Quu_reg);
        if d ~= 0
            diverge  = p; % at which backward step out of N does the inversion fail
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
        
        Vxx(:,:,p)=0.5*(Vxx(:,:,p)'+Vxx(:,:,p)); %  symmetrize operation
        dV= dV  + [k_p'*Qu, 0.5*k_p'*Quu*k_p]; 
          
    end    

      
    
    
function [X_r c_total]=Rollout(DYNCST,x0,U_r)
    %U is tall
    X_r=zeros(size(x0,1),size(U_r,2)+1);
    X_r(:,1)=x0;
    c_total=0;
    for i=1:size(U_r,2)
        [X_r(:,i+1),c]=DYNCST(X_r(:,i),U_r(:,i));
        c_total=c_total+c;
    end
    U_r(:,end)=nan;
    [~,c]=DYNCST(X_r(:,end),U_r(:,end)); %extracting the last cost
    c_total=c_total+c;
