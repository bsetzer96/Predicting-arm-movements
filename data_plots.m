spikes=trial.spikes;
handPos=trial.handPos;
trialID=trial.trialId;

%% plot movement
vert_pos=handPos(1,:);
horz_pos=handPos(2,:);
perp_pos=handPos(3,:);

figure(); plot3(vert_pos, horz_pos, perp_pos, 'lineWidth',3); xlabel('horizontal'); ylabel('vertical'); zlabel('perpendicular')
ax=gca
set(ax, 'FontSize', 35)

%% plot spike trains
figure()
for i=1:96
    spike_times=find(spikes(i,:)==1);
    subplot(96,1,i); plot(spike_times, ones(length(spike_times)), '*')
    hAx = gca;
    set(hAx, 'Visible','off')
end
 xlabel('Time (mS)'); 

 figure()
 n=10;
for i=1:n
    subplot(n,1,i); plot(spikes(i,:))
    hAx = gca;
    set(hAx, 'Visible','off')
end
 xlabel('Time (mS)'); 

%% plotting smaller amount of time:
timeMax=100;
timeRange=1:timeMax;
 figure()
figure()
for i=1:10
    spike_times=find(spikes(i,timeRange)==1);
    subplot(10,1,i); plot(spike_times, ones(length(spike_times)), '*')
    hAx = gca;
    set(hAx, 'Visible','off')
end
 xlabel('Time (mS)'); 
 xlabel('Time (mS)'); 