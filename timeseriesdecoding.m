function result = timeseriesdecoding(data,labels,varargin)
    %%
    % function result = TIMESERIESDECODING(data, labels, ...)
    %     perform 1-dimensional cross decoding on the data by training and testing a classifier on
    %     all time points
    % 
    % INPUT ARGUMENTS:  
    %
    % data
    %     ntrials*ntimepoints*ncomponents megdata matrix of the trials that
    %     will be used to train the classifier
    % labels
    %     vector of length ntrials, containing class labels
    %
    % OUTPUT:
    %
    % result
    %     struct containing for each time point, the cross-validated 
    %     classifier predictions, confusion matrices, posterior (distance)
    %     and measures of accuracy
    %
    % OPTIONAL ARGUMENTS:
    %
    % testdata
    %     trial*time*components megdata matrix of the data
    %     that will be used to test the classifier, where components must be
    %     equal to the train data. For example if the classifier has to train
    %     on data from one condition and generalize to another
    % testlabels
    %     vector of length trials, containing 0 and 1 for class labels
    %     to test on, where 0 will be left out of the classification
    % exemplarlabels
    %     vector of length trials, containing exemplar labels of the data,
    %     to perform leave one exemplar out cross validation
    % testexemplarlabels
    %     vector of length trials, containing exemplar labels of the test
    %     data
    % cvfolds
    %     the number of folds that will be performed in cross-validation (default:10),
    %     use 0 to not perform cross-validation (e.g. if train and testset are 
    %     independent)
    % weights
    %     store the classifier weights for each timepoint (default: false)
    %     Setting this to true will also store corrected weights using the
    %     covariance in the data, to allow for better interpretation.
    %     Weights are computed using all available data for each timepoint.
    % classifier
    %     which classifier type to use (default: diagLinear) 
    %     see HELP FITCDISCR for the options
    % windowsize
    %     the size of the sliding window, how many time points to take
    %     before and including t. Default: 1
    % standardize
    %     standardize data (default: true) 
    % anovaselectfeatures
    %     if 1, select features using anovas
    % pcavariance
    %     if >0, select features using pca, retaining 'pcavariance' percentage of
    %     the explained variance (value between 0-100)
    % parallel
    %     use parallel computation (default: true)
    % verbose
    %     more verbose output (default: false)
    %
    % Tijl Grootswagers
    
    %% parse optional arguments
    param = inputParser();
    addOptional(param,'timevect',[]);
    addOptional(param,'windowsize',1);
    addOptional(param,'cvfolds',10);
    addOptional(param,'parallel',true);
    addOptional(param,'verbose',false);
    addOptional(param,'weights',false);
    addOptional(param,'standardize',1);
    addOptional(param,'anovaselectfeatures',0);
    addOptional(param,'pcavariance',0);
    addOptional(param,'classifier','diaglinear');
    addOptional(param,'exemplarlabels',[]);
    addOptional(param,'testdata',[]);
    addOptional(param,'testlabels',[]);
    addOptional(param,'testexemplarlabels',[]);
    parse(param,varargin{:});
        
    %get params
    timevect = param.Results.timevect;
    windowsize = param.Results.windowsize; if windowsize<1;windowsize=1;end
    cvfolds = param.Results.cvfolds;
    verbose = param.Results.verbose;
    classifier = param.Results.classifier;
    parallel = param.Results.parallel;
    standardize = param.Results.standardize;
    anovaselectfeatures = param.Results.anovaselectfeatures;
    pcavariance = param.Results.pcavariance;
    storeweights = param.Results.weights && length(unique(labels))==2;
    
    %init traindata
    traindata = data;
    trainlabels = double(labels);
    trainexemplarlabels = double(param.Results.exemplarlabels);
    
    %init testdata
    testdata = param.Results.testdata;
    if isempty(testdata);
        testdata = traindata;
    end
    testlabels = double(param.Results.testlabels);
    if isempty(testlabels);
        testlabels = trainlabels;
    end
    testexemplarlabels = double(param.Results.testexemplarlabels);
    if isempty(testexemplarlabels);
        testexemplarlabels = trainexemplarlabels;
    end
    
    rng('shuffle');
    
    %if we are parallel, count the workers
    if parallel
        try
            p = gcp();
            parallel = p.NumWorkers;
        catch e
            fprintf('Cannot start parallel pool. Message: %s\n',e.message)
            parallel = 0;
        end
    else
        parallel = 0;
    end
    
    result = {};
    result.windowsize = windowsize;
    result.classifier = classifier;
    result.pcavariance = pcavariance;
    if ~isempty(timevect)
        result.timevect = timevect;
    end
    
    %% shuffle data
    % if train & test are the same size, shuffle them the same way,
    % otherwise, shuffle them independently
    
    if verbose
        fprintf('Shuffle data..\n');
    end
    idx1 = randperm(length(trainlabels));
    if length(trainlabels)==length(testlabels)
        idx2 = idx1;
    else
        idx2 = randperm(length(testlabels));
    end
    traindata = traindata(idx1,:,:);
    trainlabels = trainlabels(idx1);
    testdata = testdata(idx2,:,:);
    testlabels = testlabels(idx2);
    result.trainlabels = trainlabels;
    result.testlabels = testlabels;
    result.classes = unique(trainlabels);
    nclasses = length(result.classes);
        
    %% set up crossval folds
    if ~isempty(trainexemplarlabels) %leave one exemplar out
        result.cvmethod = 'leaveoneexemplarout';
        trainexemplarlabels = trainexemplarlabels(idx1);
        testexemplarlabels = testexemplarlabels(idx2);
        result.trainexemplarlabels = trainexemplarlabels;
        result.testexemplarlabels = testexemplarlabels;
        result.exemplars = unique(trainexemplarlabels);
        trainidx = zeros(length(result.exemplars),length(trainexemplarlabels));
        testidx = zeros(length(result.exemplars),length(testexemplarlabels));
        for c=1:length(result.exemplars)
            trainidx(c,:) = trainexemplarlabels~=result.exemplars(c);
            testidx(c,:) = testexemplarlabels==result.exemplars(c);
        end
    elseif cvfolds>1
        result.cvmethod = sprintf('%i-fold',cvfolds);
        cv = cvpartition(trainlabels,'kfold',cvfolds);
        trainidx = zeros(cv.NumTestSets,length(trainlabels));
        testidx = zeros(size(trainidx));
        for c=1:cv.NumTestSets
            trainidx(c,:) = cv.training(c);
            testidx(c,:) = cv.test(c);
        end
    else %no crossval
        result.cvmethod = 'none';
        trainidx = ones(1,length(trainlabels));
        testidx = ones(1,length(testlabels));
    end
    if verbose
        fprintf('Using cvmethod: %s\n',result.cvmethod);
    end
    trainidx = logical(trainidx);
    testidx = logical(testidx);
    result.cvtrainidx = trainidx;
    result.cvtestidx = testidx;
      
    %% CLASSIFY
    
    if verbose
        fprintf('Classifying..\n');
    end
    %init empty predictions
    predictions = zeros(size(testlabels,1),size(traindata,2))-1;
    weights = zeros(size(traindata,2),size(traindata,3));
    posterior = zeros(size(testlabels,1),size(traindata,2),nclasses);
    correctedweights = zeros(size(traindata,2),size(traindata,3));
    nfeatures = zeros(size(traindata,2),size(trainidx,1));
    %loop over time points
    if verbose==2
        parfor_progress(size(traindata,2));
    end
    parfor (timepoint=1:size(traindata,2),parallel)
    %for (timepoint=1:size(traindata,2)) %for debugging
        %compute the time points in the current sliding window
        slidingwindow = timepoint-windowsize+1:timepoint;
        slidingwindow(slidingwindow<1)=[];
        %get the data for this time point
        train = squeeze(traindata(:,slidingwindow,:)); %#ok<PFBNS>
        test = squeeze(testdata(:,slidingwindow,:)); %#ok<PFBNS>
        %init empty predictions
        pred = zeros(size(testlabels))-1;
        post = zeros(size(testlabels,1),nclasses);
        c=[]; %#ok<NASGU>
        if storeweights
            X = train(:,:);
            w = zeros(1,size(X,2));
            cw = zeros(1,size(X,2));
            C=ones(size(X,2));
            % standardize data
            if standardize;
                [X,mu,sig] = zscore(X);
            end
            % feature selection by pca
            if pcavariance>0;
                [C, pcaX, ~, ~, e, mu] = pca(X);
                r = cumsum(e)<=pcavariance;
                X = pcaX(:,r);
            end
            c = fitcdiscr(X,trainlabels,'DiscrimType',classifier,'Prior','uniform');
            w(1:size(X,2)) = c.Coeffs(1,2).Linear;
            if strcmpi(classifier,'diaglinear');
                cw(1:size(X,2)) = w(1:size(X,2));
            else
                cw(1:size(X,2)) = cov(X) * w(1:size(X,2))' * inv(cov(w(1:size(X,2)) * X')); %#ok<MINV>
            end
            %transform back from pca space
            weights(timepoint,:) = w * C';
            correctedweights(timepoint,:) = cw * C';
        end
        foldfeatures = zeros(1,size(trainidx,1));
        for fold = 1:size(trainidx,1)
            %grab training data for fold
            X = train(trainidx(fold,:),:);
            %grab labels for fold
            Y = trainlabels(trainidx(fold,:));
            %grab testdata for fold
            Z = test(testidx(fold,:),:); %#ok<PFBNS>
            % standardize data
            if standardize;
                [X,mu,sig] = zscore(X); %compute mu and sig on the training data
                Z = (Z-repmat(mu,size(Z,1),1))./repmat(sig,size(Z,1),1);
            end
            % feature selection by anova
            if anovaselectfeatures && anovaselectfeatures < size(X,2)
                explainedvar = zeros(1,size(X,2));
                unexplainedvar = zeros(1,size(X,2));
                for i=result.classes'
                    dat = X(Y==i,:);
                    mu = mean(dat);
                    ni = size(dat,1);
                    explainedvar = explainedvar + ni.*((mu-mean(X)).^2)./nclasses-1;
                    unexplainedvar = unexplainedvar + sum(((dat-repmat(mean(dat),ni,1)).^2)./(length(Y)-nclasses));
                end
                F = explainedvar./unexplainedvar;
                
                [~,f] = sort(F,'descend');
                features = f(1:anovaselectfeatures);
                X = X(:,features);
                Z = Z(:,features);
            end
            % feature selection by pca
            if pcavariance>0;
                [C, pcaX, ~, ~, e, mu] = pca(X);
                r = cumsum(e)<=pcavariance;
                pcaZ = (Z-repmat(mu,size(Z,1),1))*C;
                X = pcaX(:,r);
                Z = pcaZ(:,r);
            end
            foldfeatures(fold) = size(X,2);
            %train a classifier
            switch classifier
                case 'svm';
                    c = fitcsvm(X,Y,'KernelFunction','linear','Prior','uniform','Standardize',1);
                case 'svmpolynomial';
                    c = fitcsvm(X,Y,'KernelFunction','polynomial','Prior','uniform','Standardize',1);
                case {'correlation','spearman'};
                    muX=zeros(nclasses,size(X,2));
                    muY=unique(Y);
                    for cc=1:nclasses
                        muX(cc,:) = mean(X(Y==muY(cc),:),1);
                    end
                    c = fitcknn(muX,muY,'Distance',classifier,'Prior','uniform');
                otherwise; %e.g. default diaglinear
                    c = fitcdiscr(X,Y,'DiscrimType',classifier,'Prior','uniform');
            end
            post(testidx(fold,:),1) = c.Coeffs(1,2).Const + Z*c.Coeffs(1,2).Linear;
            post(testidx(fold,:),2) = c.Coeffs(2,1).Const + Z*c.Coeffs(2,1).Linear;
            pred(testidx(fold,:)) = c.predict(Z);
        end
        %store number of features
        nfeatures(timepoint,:) = foldfeatures;
        %store predictions
        predictions(:,timepoint) = pred;
        %store distances (posterior)
        posterior(:,timepoint,:) = post;
        %report progress
        if verbose==2
            parfor_progress();
        end
    end
    
    if verbose==2
        parfor_progress(0);
    end
    if verbose
        fprintf('Writing results..');
    end
    %shuffle back
    [~,trainidx] = sort(idx1);
    [~,testidx] = sort(idx2);
    result.trainlabels = result.trainlabels(trainidx);
    result.testlabels = result.testlabels(testidx);
    result.cvtrainidx = result.cvtrainidx(:,trainidx);
    result.cvtestidx = result.cvtestidx(:,testidx);
    result.nfeatures = nfeatures;
    
    %store predictions and results
    result.predictions = predictions(testidx,:);
    result.correct = result.predictions==repmat(result.testlabels,1,size(result.predictions,2));
    result.pcorr = mean(result.correct);
    result.classpcorr = zeros(length(result.classes),length(result.pcorr));
    for i=1:length(result.classes);
        result.classpcorr(i,:) = mean(result.correct(result.testlabels==result.classes(i),:));
    end
    result.balancedpcorr = mean(result.classpcorr);
    if storeweights
        result.weights = weights;
        result.correctedweights = correctedweights;
    end
    %posterior (distance/confidence)
    result.posterior = posterior(testidx,:,:);
    for i=1:size(result.posterior,1)
        result.targetposterior(i,:) = result.posterior(i,:,result.classes==result.testlabels(i));
    end
    
    %confusion matrices and posterior
    for t=1:size(result.predictions,2)
        result.CM(:,:,t) = confusionmat(result.testlabels,result.predictions(:,t));
    end
    
    %per cvfold performance, for leave one exemplar out
    for c=1:size(result.cvtestidx,1)
        result.cvpcorr(c,:) = mean(result.correct(result.cvtestidx(c,:),:));
        if isempty(setdiff(result.classes,[0 1]))
            %treat as signal detection task, so compute dprime    
            % per fold dpr
            pcv = result.predictions(result.cvtestidx(c,:),:,:);
            tlab = repmat(result.testlabels(result.cvtestidx(c,:)),1,size(pcv,2));
            hr = sum(pcv==1 & tlab==1)./sum(tlab==1);
            hr(hr==1) = 1-1/length(result.testlabels);
            hr(hr==0) = 1/length(result.testlabels);
            fa = sum(pcv==1 & tlab==0)./sum(tlab==0);
            fa(fa==1) = 1-1/length(result.testlabels);
            fa(fa==0) = 1/length(result.testlabels);
            result.cvdpr(c,:) = norminv(hr)-norminv(fa);
        end
    end
        
    if isempty(setdiff(result.classes,[0 1]))
        %treat as signal detection task, so compute dprime    
        tlab = repmat(result.testlabels,1,size(result.predictions,2));
        hr = sum(result.predictions==1 & tlab==1)./sum(tlab==1);
        hr(hr==1) = 1-1/length(result.testlabels);
        hr(hr==0) = 1/length(result.testlabels);
        fa = sum(result.predictions==1 & tlab==0)./sum(tlab==0);
        fa(fa==1) = 1-1/length(result.testlabels);
        fa(fa==0) = 1/length(result.testlabels);
        result.dpr = squeeze(norminv(hr)-norminv(fa));
    end
    
    if verbose
        fprintf('done\n');
    end
    
    
    
    
    
    
    
    