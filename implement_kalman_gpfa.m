%% load data
%arm trajectories
load('dec_G20040508.mat')

%neural trajectories
load(['seqTrain8.mat'])
trial_gpfa=reshape(seqTrain,182,8);
% traj=seqTrain';
% trial_gpfa=traj;
% 
% for i=2:8
%     load(['mat_results/mat_results/run0' num2str(i+20) '/gpfa_xDim08.mat'])
%     traj=seqTrain';
%     trial_gpfa(:,i)=traj;
% end


%% train 
%bin armPos to match gpfa
time_step=10;
trial_binned=struct;
for i=1:size(trial,1)
    for j=1:size(trial,2)
        handPos=trial(i,j).handPos;
        binned_handPos=bin_data(handPos, time_step);
        neural_traj=trial_gpfa(i,j).xorth;
        trial_binned(i,j).spikes=neural_traj;
        trial_binned(i,j).handPos=binned_handPos;
    end
end



%%
p=1400; %number of trials to train on
[A, Q, C, R, pi, v]=train_kalman(trial_binned, p);


R=R+0.001*eye(size(R));

%% test & estimate variables

n=size(trial,1)*size(trial,2);
trail_estimate=struct;

for i=p:n
    spikes=trial_binned(i).spikes;
    handPos_true=trial_binned(i).handPos;

    mu_onestep_old=pi;
    Sigma_onestep_old=v;
    T=length(handPos_true);
    %save determinant of sigma to see how it evolves over time
    %should come to some relatively fixed point
    % will give us an idea of undertainty over time
    sigma_det=zeros(1,T);
    trail_estimate(i).handPos(:,1)=pi;
    for j=2:T
        [mu_update, Sigma_update] = kalman_filter(spikes(:,j), mu_onestep_old, Sigma_onestep_old, A, Q, R, C);
        trail_estimate(i).handPos(:,j)=mu_update;
        mu_onestep_old=mu_update;
        Sigma_onestep_old=Sigma_update;
        sigma_det(j)=det(Sigma_update);
        sigma_norm(j)=norm(Sigma_update);
    end
    trail_estimate(i).sigma_det=sigma_det;
    trail_estimate(i).sigma_norm=sigma_norm;

end

% %% plot determinent
% figure()
% for i=p:p+10
%     sigma_det=trail_estimate(i).sigma_det;
%     subplot(11,1,i-p+1); plot(sigma_det)
% end
% %%
% figure() 
% 
% for i=p:p+3
%     sigma_norm=trail_estimate(i).sigma_norm;
%     subplot(4,1,i-p+1); plot(log(sigma_norm))
% end

%% plot estimates
for i=p+50:p+100
    est=trail_estimate(i).handPos;
    %smooth_est=smoothdata(est,100);
    true=trial_binned(i).handPos;
    figure()
    %plot_hand_pos(smooth_est(:,1:10:end));
    plot_hand_pos(est);
    plot_hand_pos(true);
end


%% variance explained


r_sq=zeros(n-p,1);
ind_vec=p:n;
for i = 1:(n-p)
    k=ind_vec(i);
    est=trail_estimate(k).handPos;
    true=trial_binned(k).handPos;
    r=true-est;
    r_sum=sum(sum(r.^2));
    mean_y=mean(true,2);
    y_sum=sum(sum((true-mean_y).^2));
    r_sq(i)=1-r_sum/y_sum;
end

%average variance explained
avg_r_sq=mean(r_sq);

fprintf('%f of variance explained\n', avg_r_sq*100)

%% use kalman filter toolbox and use their test data
% 
% Ts = -1;
% sys = ss(A,[B B],C,D,Ts,'InputName',{'dir 1' 'dir 2' 'dir3'},'OutputName','y');  % Plant dynamics and additive input noise w
% 
% kalmf.InputName={'z', 'x'};
% kalmf.OutputName='ze';
% 
% sys = ss(A, B, C, D, Ts,'InputName','y','OutputName','z');

