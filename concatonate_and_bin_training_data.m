function [spikes_train, handPos_train]=concatonate_and_bin_training_data(p, time_step, trial)

spikes_train=[];
handPos_train=[];
for i=1:p
    sp_t=trial(p).spikes;
    hp_t=trial(p).handPos;
    sp_binned=bin_data(sp_t, time_step);
    hp_binned=bin_data(hp_t, time_step);
    spikes_train=[spikes_train, sp_binned];
    handPos_train=[handPos_train hp_binned];
end