function [adata,alabels] = averagetrials(data,labels,average)
    %% function [averageddata, averagedlabels] = AVERAGETRIALS(data, labels, number)
    %    Average data to make supertrials. For each unique class in labels,
    %    find data belonging to this class and average these trials into
    %    supertrials. The number of supertrials depends on the number of
    %    trials per class in the data. The data is truncated to fit the
    %    maximum number of supertrials.
    %    For example, if class 1 has 124 trials and class 2 has 117,
    %    averaging 10 trials would result in 12 supertrials of class 1, and
    %    11 supertrials of class 2.
    %
    % INPUT ARGUMENTS:
    %
    % data
    %     ntrials*ntimepoints*nfeatures matrix of trials
    % labels
    %     vector of length ntrials of class labels to be used for averaging
    % average
    %     number of trials to average
    %
    % OUTPUT:
    %
    % averageddata
    %     naveragedtrials*ntimepoints*nfeatures matrix of trials
    % averagedlabels
    %     vector of length naveragedtrials of class labels belonging to the
    %     averaged trials
    %     NOTE: output data will be sorted by ascending label
    %
    %
    % Tijl Grootswagers
    % 17-03-2015

    % simple case
    if average<=1
        adata = data;
        alabels = labels;
        return
    end
    
    %sort data based on label
    [slabels,idx] = sort(labels(:));
    sdata = data(idx,:,:);
    %compute the numbers and the bin locations
    counts = histc(slabels,unique(slabels));
    if any(counts<average)
        warning('Not all classes have enough trials! These classes will not be included.')
    end;
    %starts of the bins
    start = 1+cumsum(counts)-counts;
    %ends of the bins
    stop = start + average*floor(counts/average) -1; 
    %some trials won't fit
    retain = cell2mat(arrayfun(@(a,b) a:b,start,stop,'UniformOutput',0)');
    slabels = slabels(retain);
    sdata = sdata(retain,:,:);
    %do the averaging
    adata = squeeze(mean(reshape(sdata,average,size(sdata,1)/average,size(sdata,2),size(sdata,3))));
    alabels = slabels(1:average:end);
        
end

    
