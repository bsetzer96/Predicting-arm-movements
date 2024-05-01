function [binned_data] = bin_data(data, time_step)

    T=size(data,2);
    bin_start=1:time_step:(T-time_step);
    num_bins=length(bin_start);
    binned_data=zeros(size(data,1), num_bins);
    %step through and average in time bins
    for j=1:num_bins
        time_start=bin_start(j);
        binned_data(:,j)=mean(data(:,time_start:(time_start+time_step)),2);
    end