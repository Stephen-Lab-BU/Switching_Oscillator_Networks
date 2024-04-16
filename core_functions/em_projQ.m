function [mle_Q,X_RTS,SW,Q_func] = em_projQ(y,tol,iter,A,H,Q,R,C,X_0)

% Correlated Noise Model
% =========================================================================
% Author:          Wan-Chi Hsin
% Date:            April 16, 2024
% MATLAB Version:  R2021b
% This script is for EM on Q (Sigma) matrices with switching components.
% It includes filter (skf.m) and smoother (smoother.m) steps, so
% make sure these two functions can be called.

% Model:
% S_t: Markov chain, P(St=j|St-1) = Zij, transition matrix
% X_t = A * X_t-1 + u_t  , u_t ~ N(0, Qj)  -- oscillatory latent state
% y_t = B * X_t + v_t    , v_t ~ N(0, R)   -- observation

% Input:
% y: observations
% tol: tolerance
% iter: number of iteration
% X_0: initial value of the latent state

% Output:
% mle_Q: estimated Q matrices
% X_RTS: E(Xt|y_1:T)
% SW: estimated prob of the switching states
% =========================================================================

M = size(C,1);
dim_state = size(A,1);
T = size(y,1);

pi0 = ones(1,M)/M;

Qj = zeros(dim_state,dim_state,M,iter+1);
Qj_new = zeros(dim_state,dim_state,M,iter+1);

for i=1:M
    Qj(:,:,i,1)=Q(:,:,i); 
end
k = dim_state/2;

val = zeros(T,M);
Q_func = zeros(iter,1);
cst=zeros(M,1);
ins = zeros(dim_state,dim_state,T,M);
SW = zeros(T,M,iter+1);

[W_j,X_j,V_j,KT] = skf(y,A,H,Qj(:,:,:,1),R,X_0,C,pi0);
[M_j,X_RTS,V_RTS,V_cov] = smoother(y,A,H,Qj(:,:,:,1),R,C,X_j,V_j,W_j,KT);

SW(:,:,1) = M_j;

for itr=1:iter
    % E-step
    for j=1:M
        A_j = A(:,:,j);
        Q_j = Qj(:,:,j,itr);
        R_j = R(:,:,j);
        B_j = H(:,:,j);
        
        for t=2:T
            val(t,j) = SW(t,j,itr) * (log(det(Q_j)) + log(det(R_j)) + ...
                trace(inv(R_j)*(y(t,:)'*y(t,:) + ...
                    B_j*(V_RTS(:,:,t) + X_RTS(:,t)*X_RTS(:,t)')*B_j' - ...
                    B_j*X_RTS(:,t)*y(t,:) - y(t,:)'*X_RTS(:,t)'*B_j')) + ...
                trace(inv(Q_j)*(V_RTS(:,:,t)+X_RTS(:,t)*X_RTS(:,t)' - ...
                    A_j*(V_cov(:,:,t)+X_RTS(:,t)*X_RTS(:,t-1)')' - ...
                    (V_cov(:,:,t)+X_RTS(:,t)*X_RTS(:,t-1)')*A_j' + ...
                    A_j*(V_RTS(:,:,t-1) + X_RTS(:,t-1)*X_RTS(:,t-1)')*A_j')));
        end
    end
    lik = -0.5*sum(val,'all'); 
    
    Q_func(itr) = lik;
    fprintf('iter %g Q function %g \n', itr, lik);
    
    % M-step
    for j=1:M
        A_j = A(:,:,j);
        cst(j) = sum(SW(:,j,itr));
        
        for t=2:T
            ins(:,:,t,j) = SW(t,j,itr)*((V_RTS(:,:,t)+X_RTS(:,t)*X_RTS(:,t)') - ...
                A_j*(V_cov(:,:,t)+X_RTS(:,t)*X_RTS(:,t-1)')' - ...
                (V_cov(:,:,t)+X_RTS(:,t)*X_RTS(:,t-1)')*A_j' + ...
                A_j*(V_RTS(:,:,t-1) + X_RTS(:,t-1)*X_RTS(:,t-1)')*A_j');
        end
        Qj(:,:,j,itr+1) = (1/cst(j)) * sum(ins(:,:,:,j),3);
    end

    % Projection
    % off-diagonal
    tmp = Qj(:,:,:,itr+1);
    scaling = zeros(k,k,M);
    for u=1:M
        for i=1:k
            for j=1:k
                if i<j
                    [U, S, V] = svd(tmp(2*i-1:2*i,2*j-1:2*j,u));
                    scaling(i,j,u) = sqrt(S(1,1)*S(2,2));
                    Qj_new(2*i-1:2*i,2*j-1:2*j,u,itr+1) = scaling(i,j,u)*U*V';
                elseif i==j
                    Qj_new(2*i-1:2*i,2*j-1:2*j,u,itr+1) = sqrt(Qj(2*i-1,2*j-1,u,itr+1)*Qj(2*i,2*j,u,itr+1))*eye(2);
                end
            end
        end
    end

    for u=1:M
        for i=1:k
            for j=1:k
                if i>j
                    Qj_new(2*i-1,2*j-1,u,itr+1) = Qj_new(2*j-1,2*i-1,u,itr+1);
                    Qj_new(2*i-1,2*j,u,itr+1) = Qj_new(2*j,2*i-1,u,itr+1); 
                    Qj_new(2*i,2*j-1,u,itr+1) = Qj_new(2*j-1,2*i,u,itr+1); 
                    Qj_new(2*i,2*j,u,itr+1) = Qj_new(2*j,2*i,u,itr+1); 
                end
            end
        end
    end
                    
    if(abs(Qj(:,:,1,itr+1)-Qj(:,:,1,itr)) < tol)
        break
    end
    [W_j,X_j,V_j,KT] = skf(y,A,H,Qj_new(:,:,:,itr+1),R,X_0,C,pi0);
    [M_j,X_RTS,V_RTS,V_cov] = smoother(y,A,H,Qj_new(:,:,:,itr+1),R,C,X_j,V_j,W_j,KT);
    SW(:,:,itr+1) = M_j;
end
mle_Q = Qj_new(:,:,:,2:itr+1);
end
