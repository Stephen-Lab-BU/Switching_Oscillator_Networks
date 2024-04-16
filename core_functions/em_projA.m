function [mle_A,X_RTS,SW,Q_func] = em_projA(y,tol,iter,A,B,Q,R,Z,X_0)

% Directed Influence Model
% =========================================================================
% Author:          Wan-Chi Hsin
% Date:            April 16, 2024
% MATLAB Version:  R2021b
% This script is for EM on A matrices with switching components.
% It includes filter (skf.m) and smoother (smoother.m) steps, so
% make sure these two functions can be called.

% Model:
% S_t: Markov chain, P(St=j|St-1) = Zij, transition matrix
% X_t = Aj * X_t-1 + u_t , u_t ~ N(0, Q)  -- oscillatory latent state
% y_t = B * X_t + v_t    , v_t ~ N(0, R)  -- observation

% Input:
% y: observations
% tol: tolerance
% iter: number of iteration
% X_0: initial value of the latent state

% Output:
% mle_A: estimated A matrices
% X_RTS: E(Xt|y_1:T)
% SW: estimated prob of the switching states
% =========================================================================

M = size(Z,1); 
dim_state = size(A,1);
T = size(y,1);

pi0 = ones(1,M)/M;   

Aj = zeros(dim_state,dim_state,M,iter+1);
Aj_new = zeros(dim_state,dim_state,M,iter+1);
% set up initial A
for i=1:M
    Aj(:,:,i,1) = A(:,:,i);
end
k = dim_state/2;

A1 = zeros(dim_state,dim_state,M,T);
A2 = zeros(dim_state,dim_state,M,T);
A1sum=zeros(dim_state,dim_state,T);
A2sum=zeros(dim_state,dim_state,T);

val = zeros(T,M);
Q_func = zeros(iter,1);
SW = zeros(T,M,iter+1);

[W_j,X_j,V_j,KT] = skf(y,Aj(:,:,:,1),B,Q,R,X_0,Z,pi0);
[M_j,X_RTS,V_RTS,V_cov] = smoother(y,Aj(:,:,:,1),B,Q,R,Z,X_j,V_j,W_j,KT);

SW(:,:,1) = M_j;

for itr=1:iter
    % E-step
    for j=1:M
        A_j = Aj(:,:,j,itr);
        Q_j = Q(:,:,j);
        R_j = R(:,:,j);
        B_j = B(:,:,j);
        
        for t=2:T
            val(t,j) = SW(t,j,itr) * (log(det(Q_j)) + log(det(R_j)) + ...
                trace(inv(R_j)*(y(t,:)'*y(t,:) + ...
                    B_j*(V_RTS(:,:,t) + X_RTS(:,t)*X_RTS(:,t)')*B_j' - ...
                    B_j*X_RTS(:,t)*y(t,:) - y(t,:)'*X_RTS(:,t)'*B_j')) + ...
                trace(inv(Q_j)*(V_RTS(:,:,t)+X_RTS(:,t)*X_RTS(:,t)' - ...
                    A_j*(V_cov(:,:,t)+X_RTS(:,t-1)*X_RTS(:,t)') - ...
                    (V_cov(:,:,t)+X_RTS(:,t)*X_RTS(:,t-1)')*A_j' + ...
                    A_j*(V_RTS(:,:,t-1) + X_RTS(:,t-1)*X_RTS(:,t-1)')*A_j')));           
        end
    end
    lik = -0.5*sum(val,'all');  
    Q_func(itr,1) = lik;

    fprintf('iter %g Q function %g \n', itr, lik);

    % M-step 
    for j=1:M
        for t=2:T
            A2(:,:,j,t) = M_j(t,j) * (V_cov(:,:,t) + X_RTS(:,t)*X_RTS(:,t-1)');
            A1(:,:,j,t) = M_j(t,j) * (V_RTS(:,:,t-1) + X_RTS(:,t-1)*X_RTS(:,t-1)');
        end
        A1sum(:,:,j) = sum(A1(:,:,j,:),4);
        A2sum(:,:,j) = sum(A2(:,:,j,:),4);
        
        Aj(:,:,j,itr+1) = A2sum(:,:,j)*inv(A1sum(:,:,j));
    end

    % Projection
    % off-diagonal
    tmp = Aj(:,:,:,itr+1);
    scaling = zeros(k,k,M);
    for u=1:M
        for i=1:k
            for j=1:k
                if i~=j
                    [U, S, V] = svd(tmp(2*i-1:2*i,2*j-1:2*j,u));
                    scaling(i,j,u) = sqrt(S(1,1)*S(2,2));
                    Aj_new(2*i-1:2*i,2*j-1:2*j,u,itr+1) = scaling(i,j,u)*U*V';
                end
            end
        end
    end

    % diagonal: plus scaling, projection, subtract scaling
    row_sum = sum(scaling,2);
    for u=1:M
        for i=1:k
              [U, S, V] = svd(tmp(2*i-1:2*i,2*i-1:2*i,u) + row_sum(i,u)*eye(2));
              Aj_new(2*i-1:2*i,2*i-1:2*i,u,itr+1) = sqrt(S(1,1)*S(2,2))*U*V' - row_sum(i,u)*eye(2);
        end
    end

    if(abs(Aj(:,:,1,itr+1)-Aj(:,:,1,itr)) < tol) 
        break
    end
    [W_j,X_j,V_j,KT] = skf(y,Aj_new(:,:,:,itr+1),B,Q,R,X_0,Z,pi0);
    [M_j,X_RTS,V_RTS,V_cov] = smoother(y,Aj_new(:,:,:,itr+1),B,Q,R,Z,X_j,V_j,W_j,KT);      
    SW(:,:,itr+1) = M_j;
end
mle_A = Aj_new(:,:,:,2:itr+1);
end
