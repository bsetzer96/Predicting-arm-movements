%% GPFA on trajectory data

addpath(genpath('gpfa_v0203'));

%%
%load data
load('dec_G20040508.mat')
%load('trial_cropped.mat')

for i=1:1410
    trial(i).spikes=trial_cropped(i).spikes;
    trial(i).handPos=trial_cropped(i).handPos;
end

trial_reshape=trial(:,1);
%change shape 
for i=2:8
    trial_traj=trial(:,i);
    rows=182*i;
    trial_reshape(rows-181:rows,1)=trial_traj;
end
%%


% n=length(trial);
% trial_cropped=struct;
% for i=1:n
%     for j=1:8
%         spikes=trial(i).spikes;
%         runIdx=trial(i).trialId;
%         sum_spikes=sum(spikes');
%         ind=sum_spikes>5;
%         spikes_more=spikes(ind,:);
%         trial_cropped(i,j).spikes=spikes_more;
%         trial_cropped(i,j).trialId=runIdx;
%     end
% end




%% fa
%lambda = factoran(spikes_more,3); %returns the factor loadings matrix lambda for the data matrix X with m common factors.


%% GPFA

method = 'gpfa';

% Select number of latent dimensions
xDim = 8;
% NOTE: The optimal dimensionality should be found using 
%       cross-validation (Section 2) below.

% If using a two-stage method ('fa', 'ppca', or 'pca'), select
% standard deviation (in msec) of Gaussian smoothing kernel.
kernSD = 30;
% NOTE: The optimal kernel width should be found using 
%       cross-validation (Section 2) below.


    runIdx=20+i;
    
    % Extract neural trajectories
    result = neuralTraj(runIdx, trial_reshape, 'method', method, 'xDim', xDim,... 
                        'kernSDList', kernSD, 'binWidth', 10);
    % NOTE: This function does most of the heavy lifting.
    
    % Orthonormalize neural trajectories
    [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
    % NOTE: The importance of orthnormalization is described on 
    %       pp.621-622 of Yu et al., J Neurophysiol, 2009.

    save('seqTrain8', 'seqTrain')


    %%
    % Plot neural trajectories in 3D space
    plot3D(seqTrain, 'xorth', 'dimsToPlot', 1:3, 'nPlotMax', 10);
