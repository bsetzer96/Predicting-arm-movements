function [mu_update, Sigma_update] = kalman_filter(x, mu_onestep_old, Sigma_onestep_old, A, Q, R, C)


%% one step predicton
mu_onestep_new=A*mu_onestep_old;
Sigma_onestep_new=A*Sigma_onestep_old*A'+Q;

%% update
%plot determinent or trace of covariance
%try not updating if it becomes singular
K=Sigma_onestep_new*C'*inv(C*Sigma_onestep_new*C'+R);

%try fixing K as constant
%K=eye(size(K));

mu_update=mu_onestep_new+K*(x-C*mu_onestep_new);
Sigma_update=Sigma_onestep_new-K*C*Sigma_onestep_new;
