%% This script was used to compare several options along the analysis pipeline for Grootswagers, Wardle, & Carlson
% It creates a massive RESULTS.mat file, that was used in makeplots.m to
% create the figures used in the paper.
%
% These scripts are supplied for reference only.
% Tijl Grootswagers, June 2015

%% find all subject ID's
subdir = '~/DATA/MEGMethodsPaper/';
files = dir([subdir '*_1000Hz_none.mat']);
nsubjects = length(files);subjects=[];
for i=1:nsubjects
    subjects{i} = files(i).name(19:20);
end

defaults.method='none';
defaults.samplerate=200;
defaults.averaging=4;
defaults.pcavariance=99;
TOTALTIME = [];

RESULTS = {};
RESULTS.dimredmethod = cell(nsubjects,1);
RESULTS.windowsize = cell(nsubjects,1);
RESULTS.samplerate = cell(nsubjects,1);
RESULTS.averaging = cell(nsubjects,1);
RESULTS.classifier = cell(nsubjects,1);
RESULTS.classifiernopca = cell(nsubjects,1);
RESULTS.cvmethod = cell(nsubjects,1);
RESULTS.crossdecoding = cell(nsubjects,1);
RESULTS.dsms = cell(nsubjects,1);

%%
for s=1:length(subjects);
    starttime = tic;
    
    %% 1.1 subsampling
    values = {50,100,200,1000};
    for i=1:length(values)
        fprintf('%s  ',datestr(now))
        fprintf('subject %i/%i - subsampling %i/%i ',s,nsubjects,i,length(values));tic;
        samplerate = values{i};
        [data,B] = loaddata(subjects{s},samplerate,defaults.method);
        [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar,defaults.averaging);
        res = timeseriesdecoding(avdata,avlabels>24,'exemplarlabels',avlabels,'pcavariance',defaults.pcavariance,'timevect',data.timevect);
        res.samplerate=samplerate;
        res.subject = subjects{s};
        RESULTS.samplerate{s,i} = res;
        fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
    end
    
    %% 1.2 sliding windows (we last loaded the 1000Hz data)
    values = {1,5,10,20};
    [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar,defaults.averaging);
    for i=1:length(values)
        fprintf('%s  ',datestr(now))
        fprintf('subject %i/%i - windowsizes %i/%i ',s,nsubjects,i,length(values));tic;
        windowsize = values{i};
        res = timeseriesdecoding(avdata,avlabels>24,'exemplarlabels',avlabels,'pcavariance',defaults.pcavariance,'windowsize',windowsize,'timevect',data.timevect);
        res.subject = subjects{s};
        RESULTS.windowsize{s,i} = res;
        fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
    end
    
    %% 1.3 averaging
    %from now on use the default (200Hz) data
    [data,B] = loaddata(subjects{s},defaults.samplerate,defaults.method);
    values = {1,4,8,16,32};
    for i=1:length(values)
        fprintf('%s  ',datestr(now))
        fprintf('subject %i/%i - averaging %i/%i ',s,nsubjects,i,length(values));tic;
        averaging = values{i};
        [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar,averaging);
        res = timeseriesdecoding(avdata,avlabels>24,'exemplarlabels',avlabels,'pcavariance',defaults.pcavariance,'timevect',data.timevect);
        res.averaging=averaging;
        res.subject = subjects{s};
        RESULTS.averaging{s,i} = res;
        fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
    end
    
    %% 2.1 none, pca, anova
    values = {'none','pca','anova'};
    %from now on use the default averaged data
    [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar,defaults.averaging);
    for i=1:length(values)
        fprintf('%s  ',datestr(now))
        fprintf('subject %i/%i - dimred %i/%i ',s,nsubjects,i,length(values));tic;
        dimredmethod = values{i};
        switch i
            case 1
                res = timeseriesdecoding(avdata,avlabels>24,'exemplarlabels',avlabels,'timevect',data.timevect);
            case 2
                res = timeseriesdecoding(avdata,avlabels>24,'pcavariance',defaults.pcavariance,'exemplarlabels',avlabels,'timevect',data.timevect);
            case 3
                res = timeseriesdecoding(avdata,avlabels>24,'anovaselectfeatures',25,'exemplarlabels',avlabels,'timevect',data.timevect);
        end
        res.dimredmethod=dimredmethod;
        res.subject = subjects{s};
        RESULTS.dimredmethod{s,i} = res;
        fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
    end
        
    %% 2.2 classifiers
    [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar,defaults.averaging);
    values = {'Linear','diagLinear','svm','Quadratic','svmpolynomial','correlation','spearman'};
    for i=1:length(values)
        fprintf('%s  ',datestr(now))
        fprintf('subject %i/%i - classifiers %i/%i ',s,nsubjects,i,length(values));tic;
        classifier = values{i};
        res = timeseriesdecoding(avdata,avlabels>24,'exemplarlabels',avlabels,'pcavariance',defaults.pcavariance,'classifier',classifier,'timevect',data.timevect);
        res.classifier=classifier;
        res.subject = subjects{s};
        RESULTS.classifier{s,i} = res;
        fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
    end
    
    %% 2.2.1 classifiers without pca
    [avdata,avlabels] = averagetrials(data.class_dat,B.exemplar,defaults.averaging);
    values = {'Linear','diagLinear','svm','correlation','spearman'};
    for i=1:length(values)
        fprintf('%s  ',datestr(now))
        fprintf('subject %i/%i - classifiers nopca %i/%i ',s,nsubjects,i,length(values));tic;
        classifier = values{i};
        res = timeseriesdecoding(avdata,avlabels>24,'exemplarlabels',avlabels,'pcavariance',0,'classifier',classifier,'timevect',data.timevect);
        res.classifier=classifier;
        res.subject = subjects{s};
        RESULTS.classifiernopca{s,i} = res;
        fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
    end
    
    %% 2.3 cross-validation
    values = {0,2,10,-1,length(avlabels)};
    for i=1:length(values)
        fprintf('%s  ',datestr(now))
        fprintf('subject %i/%i - cross-validation %i/%i ',s,nsubjects,i,length(values));tic;
        cvmethod = values{i};
        if cvmethod==-1 %% leave one out
            res = timeseriesdecoding(avdata,avlabels>24,'exemplarlabels',avlabels,'pcavariance',defaults.pcavariance,'timevect',data.timevect);
        else
            res = timeseriesdecoding(avdata,avlabels>24,'cvfolds',cvmethod,'pcavariance',defaults.pcavariance,'timevect',data.timevect);
        end
        res.cvmethod=cvmethod;
        res.subject = subjects{s};
        RESULTS.cvmethod{s,i} = res;
        fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
    end
    
    %% 3.1 cross-decoding
    fprintf('%s  ',datestr(now))
    fprintf('subject %i/%i - cross-decoding',s,nsubjects);tic;
    res = timeseriescrossdecoding(avdata,avlabels>24,'exemplarlabels',avlabels,'pcavariance',defaults.pcavariance,'timevect',data.timevect);
    res.subject = subjects{s};
    RESULTS.crossdecoding{s} = res;
    fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
        
    %% 4.1 make dsms (decoding metric)
    fprintf('%s  ',datestr(now))
    fprintf('subject %i/%i - make dsm',s,nsubjects);tic;
    DSM = zeros(141,48,48);
    for e1=1:48
        fprintf('.')
        for e2=e1+1:48
            sub = ismember(avlabels,[e1 e2]);
            r = timeseriesdecoding(avdata(sub,:,:),avlabels(sub)==e1,'cvfolds',2,'pcavariance',defaults.pcavariance,'classifier','diagLinear');
            DSM(:,e1,e2) = r.balancedpcorr;
            DSM(:,e2,e1) = r.balancedpcorr;
        end
    end
    res = {};
    res.DSM=zeros(141,1128);
    for t=1:141
        res.DSM(t,:) = squareform(squeeze(DSM(t,:,:)));
    end
    res.subject = subjects{s};
    res.timevect = data.timevect;
    RESULTS.dsms{s} = res;
    fprintf('- %s\n',datestr(toc*1/24/3600,'DD-HH:MM:SS'))
    
    %% write out
    fprintf('%s  ',datestr(now))
    fprintf('subject %i/%i - writing results\n',s,nsubjects);
    save('RESULTS.mat','RESULTS','s','defaults','subjects','-v7.3')
    
    TOTALTIME(s) = toc(starttime); %#ok<SAGROW>
    prt = [sprintf('%s  ',datestr(now))...
    sprintf('subject %i/%i ',s,nsubjects)...
    sprintf('- TIME: %s ', datestr(TOTALTIME(s)*1/24/3600,'DD-HH:MM:SS'))...
    sprintf('- TOTALTIME: %s ', datestr(sum(TOTALTIME)*1/24/3600,'DD-HH:MM:SS'))...
    sprintf('- ETA: %s\n',datestr(mean(TOTALTIME(1:s))*(nsubjects-s)*1/24/3600,'DD-HH:MM:SS'))];
    fprintf(prt);
    
    %% send mail
    if ~isempty(which('sendemail'))
        sendemail([],'Matlab status report',sprintf('Subject %i finished!\n\n%s',s,prt))
    end
end
    