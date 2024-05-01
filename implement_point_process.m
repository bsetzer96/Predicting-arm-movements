%% point process

%A and variance for the movement are the same

%estimate C
%use GLM fit
%design matrix xyz (3xT)
%response matrix is neual data
% poisson

load('trial_cropped.mat')
trial=trial_cropped;


% train 
%can concatonate trails 
p=1400; %number of trials to train on
%for i=1:1410
    %trial(i).handPos=diff(trial(i).handPos')';
    %trial(i).spikes=trial(i).spikes(:,2:end);
%end
[F, D, ~, Q, pi, v]=train_kalman(trial, p);


%R=R+0.001*eye(size(R));

%% fit C
p=5; %training data
m=4; %number of movement directions + 1, dim 4 is the background rate (rate when the hand isn't moving at all)
n=97; %number of neurons
C=zeros(m,n);

%loop through neurons and fit for each neuron
%binomial or poisson 
%will return 4 parameters for each neuron

%concatonate data
spikes_concat=[];
handPos_concat=[];
neuron=struct;
for j=1:n
    for i=1:p
        %spikes=trial(i).spikes(j,2:end);
        spikes=trial(i).spikes(j,:);
        spikes_concat=[spikes_concat spikes];
        %handPos=diff(trial(i).handPos')';
        handPos=trial(i).handPos;
        handPos_concat=[handPos_concat handPos];
    end
    [C(:,j), dev, stats]= glmfit(handPos_concat', spikes_concat',  'poisson');
    neuron(j).stats=stats;
end



%%
p_value=zeros(n,3);
for j=1:n
    p_value(j,:)=neuron(j).stats.p(2:4);
end
signf=p_value<0.05;



%% estimate arm trajectories
% %e^c0
% %e^c1 - how much firing rate modulated for each direction
% % intercept term
% q=10;
% 
% for i=p:(p+q)
%     spikes=trial(i).spikes;
%     x_update=pi;
%     W_update=v;
%     trial_estimate(i).handPos(:,1)=x_update;
%     for j=2:T
%        %one-step
%        mu_onestep_new=F*x_update;
%        W_onestep_new=F'*W_update*F+D;
%        lambda_t=zeros(n,1);
%        for k=1:n
%             lambda_t(k) = exp(C(1,k)+C(2:4,k)'*mu_onestep_new);
%        end
%        %update
%        delta_N= spikes(:,j);
%        W_update=inv(inv(W_onestep_new)+C(2:4,:)*diag(lambda_t)*C(2:4,:)');
%        x_update=mu_onestep_new+W_update*C(2:4,:)*(delta_N-lambda_t);
%        trial_estimate(i).handPos(:,j)=x_update;
%     end
% end

%%

%lambda * t = e^[c0 c1 c2 c3]T[1 x y z]for each neuron along diagonal
q=10;
%F=eye(3);
for i=p:(p+q)
    spikes=trial(i).spikes;
    x_update=pi;
    W_update=v;
    trial_estimate(i).handPos(:,1)=x_update;
    T=length(spikes);
    for j=2:T
       %one-step
       mu_onestep_new=F*x_update;
       W_onestep_new=F'*W_update*F+D;
       lambda_t=zeros(n,1);
       W_sum=0;
       x_sum=0;
       for k=1:n
           lambda_t(k) = exp(C(1,k)+C(2:4,k)'*mu_onestep_new);
           delta_N= spikes(k,j);
           W_sum=W_sum+C(2:4,k)*lambda_t(k)*C(2:4,k)';
           x_sum=x_sum+C(2:4,k)*(delta_N-lambda_t(k));    
       end
       W_update=inv(inv(W_onestep_new)+W_sum);
       %W_update=eye(3)/10000;
       x_update=mu_onestep_new+W_update*x_sum;
       trial_estimate(i).handPos(:,j)=x_update;
    end
end

%% plot estiamte versus true

for i=p:(p+10)
    est=trial_estimate(i).handPos;
    true=trial(i).handPos;
    figure()
    plot_hand_pos(true); hold on
    plot_hand_pos(est); 
end
%%
for i=p:(p+10)
    est=trial_estimate(i).handPos;
    true=trial(i).handPos;
    figure()
    plot(est(1,:), 'k'); hold on
    plot(est(2,:), 'k')
    plot(est(3,:), 'k')

    plot(true(1,:), 'r'); hold on
    plot(true(2,:), 'r')
    plot(true(3,:), 'r')
end


