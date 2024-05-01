%% load data

%load('dec_G20040508.mat')
load('trial_cropped.mat')
trial=trial_cropped;


%% train 
%can concatonate trails 
p=1400; %number of trials to train on


[A, Q, C, R, pi, v]=train_kalman(trial, p);


R=R+0.001*eye(size(R));

%% test & estimate variables

n=size(trial,1)*size(trial,2);
trail_estimate=struct;

for i=p:n
    spikes=trial(i).spikes;
    handPos_true=trial(i).handPos;
    mu_onestep_old=pi;
    Sigma_onestep_old=v;
    T=length(handPos_true);
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


%% plot estimates
for i=p+50:p+100
    est=trail_estimate(i).handPos;
    %smooth_est=smoothdata(est,100);
    true=trial(i).handPos;
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
    true=trial(k).handPos;
    r=true-est;
    r_sum=sum(sum(r.^2));
    mean_y=mean(true,2);
    y_sum=sum(sum((true-mean_y).^2));
    r_sq(i)=1-r_sum/y_sum;
end

%average variance explained
avg_r_sq=mean(r_sq);

fprintf('%f of variance explained\n', avg_r_sq*100)

