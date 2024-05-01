%% load data

load('dec_G20040508.mat')

% dont update kalman gain at all and see if youre being pushed in right
% direction by comparison from observation to expected observation

%% train on one trial
%can concatonate trails 
p=20; %number of trials to train on%timestep
time_step=30; %ms
[spikes_train, handPos_train]=concatonate_and_bin_training_data(p, time_step, trial);


%% train & estimate parameters

[A, Q, C, R]=train_kalman(spikes_train, handPos_train);
pi=mean(handPos_train,2);
v=cov(handPos_train');

%% test & estimate variables

n=length(trial);
T=length(handPos_train);

trail_estimate=struct;
n=p+5;

for i=p:n
    spikes=trial(i).spikes;
    handPos_true=trial(i).handPos;
    mu_onestep_old=pi;
    Sigma_onestep_old=v;
    T=length(handPos_true);

    %bin over timestep
    [binned_spikes] = bin_data(spikes, time_step);
    binned_handPos = bin_data(handPos_true, time_step);
    trail_estimate(i).true_handPos=binned_handPos;

    for j=1:size(binned_spikes,2)
        %estimate arm position
        [mu_update, Sigma_update] = kalman_filter(binned_spikes(:,j), mu_onestep_old, Sigma_onestep_old, A, Q, R, C);
        trail_estimate(i).handPos(:,j)=mu_update;

        %update
        mu_onestep_old=mu_update;
        Sigma_onestep_old=Sigma_update; 
    %save determinant of sigma to see how it evolves over time
    %should come to some relatively fixed point
    % will give us an idea of undertainty over time
        trail_estimate(i).sigma_det(j)=det(Sigma_update);
    end

end

%% plot

for i=p:n
    est=trail_estimate(i).handPos;
    %smooth_est=smoothdata(est,100);
    true=trial(i).handPos;
    figure()
    %plot_hand_pos(smooth_est(:,1:10:end));
    plot_hand_pos(est);
    plot_hand_pos(true);
end

%% variance explained


r_sq=zeros(n,1);
for i = p:n
    est=trail_estimate(i).handPos;
    true=trail_estimate(i).true_handPos;
    r=true-est;
    r_sum=sum(r.^2);
    mean_y=mean(true);
    y_sum=sum((true-mean_y).^2);
    r_sq(i)=1-r_sum/y_sum;
end

%average variance explained
avg_r_sq=mean(r_sq);

fprintf('%f of variance explained\n', avg_r_sq*100)

