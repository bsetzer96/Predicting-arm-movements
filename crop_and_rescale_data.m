%% loop through dataset and crop


load('dec_G20040508.mat')

%%

n=length(trial);


trial_cropped=struct;
for i=1:n
    for j=1:8
        handPos=trial(i,j).handPos;
        true_diff=zeros(3, length(handPos)-1);
        for k=1:3
            true_diff(k,:)=diff(handPos(k,:));
        end
        total_movement=sum(true_diff);
        %figure(); plot(total_movement);
        start_ind= find(abs(total_movement)>0.05);
        spikes_new=trial(i,j).spikes(:,start_ind(1):(start_ind(1)+200));
        handPos_new=trial(i,j).handPos(:,start_ind(1):(start_ind(1)+200));

        %var_pos=var(handPos_new, [], 2);

        %handPos_new=handPos_new./var_pos;

        trial_cropped(i,j).spikes=spikes_new;
        trial_cropped(i,j).handPos=handPos_new;
    end
end

save('trial_cropped.mat', 'trial_cropped')