
function []=preprocess(subdir,datadir,eventlength,projectordelay)
    %%
    % One whole preprocessing pipeline
    % @subdir: the location of the subjects' files
    % @datadir: where to write the output
    % @subsamplefactor: the factor by which the data should be subsampled,
    %    e.g. factor 5 would be re-sampling to 200Hz (=1000Hz/5)
    %%
    
    gcp;
    %% 1 read meg files
    fprintf('\nWorking in %s\n\n',subdir)
    
    behav = dir([subdir '/*.mat']);
    fprintf('1/9 Reading %s (%.2f mb)...\n',behav.name,behav.bytes/1024/1024)
    load([subdir filesep behav.name]);
    
    targetfile = [datadir filesep 'prep_' behav.name];
    dataset = struct2table(data);
    behavdata = data;
    save(targetfile);
        
    megfiles = dir([subdir '/*_B*.con']);
    MEGdata = []; %cell(length(megfiles),1);
    MEGheaders = cell(length(megfiles),1);
    for j = 1:length(megfiles)
        fname = [subdir filesep megfiles(j).name];
        fprintf('Reading %s (%.2f mb)...\n',fname,megfiles(j).bytes/1024/1024)
        MEGdata = [MEGdata getYkgwData(fname)]; %#ok<AGROW> cause loop is small
        MEGheaders{j} = getYkgwHdrAcqCond(fname);
    end
    
    p.usedtriggers = [p.triggerstimon p.triggerstimoff p.triggeranimate ...
        p.triggerinanimate p.triggerclear p.triggerblurry p.triggerexamplar];
    
    p.triggerthresh = .5;
        
    trigmatrix = zeros([size(MEGdata,2)-1 length(p.usedtriggers)]);% comes from MEG

    %% 2 make triggermatrix
    fprintf('\n2/9 Making triggermatrix')
    for triggernum = 1:length(p.usedtriggers)
        % find the triggers
        temp = (diff(MEGdata(p.usedtriggers(triggernum)+1,:)))>p.triggerthresh;
        trigmatrix(:,triggernum) = temp*p.usedtriggers(triggernum);
        fprintf('.')
    end
    fprintf('Done\n\n')
    
    exemplarvector = sum(trigmatrix(:,1:length(p.usedtriggers)),2)';
    triggertimes = find(exemplarvector);
    triggerlabels = exemplarvector(triggertimes);
    
    %clear trigmatrix;
    
    expectedtriggers = zeros(1,p.totaltrials*5); %5 triggers per stim
    for j=1:p.totaltrials
        k = 1+(j-1)*5;
        expectedtriggers(k) = p.triggerstimon;
        if data.animate(j);t=p.triggeranimate;else t=p.triggerinanimate;end
        expectedtriggers(k+1) = t;
        t = stimuli.triggernumbers(find(strcmp(data.stimulus(j), stimuli.names),1)); %#ok<FNDSB>
        expectedtriggers(k+2) = t; 
        if data.blurred(j);t=p.triggerblurry;else t=p.triggerclear;end
        expectedtriggers(k+3) = t;
        expectedtriggers(k+4) = p.triggerstimoff;
    end
    
    %% 3 remove duplicate triggers (with diff of <5 in between)   
    %dup = find(diff(triggertimes)<=5 & diff(triggerlabels)==0);   
    dup = find(diff(triggerlabels)==0);
    fprintf('3/9 Removing these %d duplicate triggers:\n',length(dup));
    disp([[triggertimes(dup);triggerlabels(dup)]' [triggertimes(dup+1);triggerlabels(dup+1)]'])
    
    triggertimes(dup+1) = [];
    triggerlabels(dup+1) = [];
    
    %% 4 check missing triggers
    fprintf('4/9 Checking for missing triggers...');
    if all(expectedtriggers(1:5)==triggerlabels(1:5))
        fprintf('First triggers match, lets roll!\n')
    else
        fprintf('First triggers do not match\n')
    	return;
    end
    
    triallist = zeros(sum(expectedtriggers==p.triggerstimon),2);
    triggerlist = zeros(length(expectedtriggers),2);
    for j = 1:length(expectedtriggers)
        if expectedtriggers(j) ~= triggerlabels(j)
            fprintf('#%d: missing trigger %d...',j,expectedtriggers(j));
            %times should be x 1000 as MEG is in ms (data is in seconds)
            if expectedtriggers(j)==p.triggerstimon
                %stimon trigger, see if the next trigger is there
                if expectedtriggers(j+1) == triggerlabels(j)
                    t = triggertimes(j)-p.timebetweentriggers*1000;
                    fprintf(' all good, the next trigger was there :)\n')
                elseif expectedtriggers(j+2) == triggerlabels(j)
                    t = triggertimes(j)-2*p.timebetweentriggers*1000;
                    fprintf(' all good, the 2nd next trigger was there :)\n')
                elseif expectedtriggers(j+3) == triggerlabels(j)
                    t = triggertimes(j)-3*p.timebetweentriggers*1000;
                    fprintf(' all good, the 3nd next trigger was there :)\n')
                else
                    fprintf(' uh oh, we seem to be missing 4 triggers in a row here. Aborting..\n')
                    return;
                end
           elseif expectedtriggers(j)==p.triggerstimoff
                dur = 1000*(dataset(j/5,:).stimoff-dataset(j/5,:).stimon);
                fprintf(' Stim lasted for %.5f, stimon at %.5f; setting stimoff at %.5f\n',...
                    dur,triggertimes(j-4),triggertimes(j-4)+dur)
                t = triggertimes(j-4)+dur;
            else
                %take previous trigger and put the missing one after it
                t = triggertimes(j-1)+p.timebetweentriggers*1000; %or NaN?
                fprintf(' all good :)\n')
            end
            triggertimes = [triggertimes(1:j-1),t,triggertimes(j:end)];
            triggerlabels = [triggerlabels(1:j-1),expectedtriggers(j),triggerlabels(j:end)];
        end 
        triggerlist(j,:) = [expectedtriggers(j) triggertimes(j)];
    end
    triallist = triggerlist(triggerlist(:,1)==p.triggerstimon,:);
    
    %% 5 trigger check
    fprintf('\n5/9 Final trigger check..');
    if all(triggerlabels==expectedtriggers)
        fprintf(' passed!\n\n');
    else
        fprintf(' failed!\n\n');
    end
    
        

    p.eventlength = eventlength;
    p.projectordelay = projectordelay; %delay of the projector, in ms


    fulldata = zeros([size(triallist,1) diff(p.eventlength)+1 160]);
    fulllist = zeros([size(triallist,1) 2]);
    %disp(p)

    %% 6 slice
    fprintf('6/9 Processing raw data')
    for trialnum = 1:length(triallist)-1
        ebegin  = triallist(trialnum,end)+p.projectordelay;
        ebounds = ebegin + p.eventlength;

        fulldata(trialnum,:,:) = MEGdata(1:160,ebounds(1):ebounds(2))';
        fulllist(trialnum,:) = triallist(trialnum,1:2);

        if ~mod(trialnum,round(length(triallist)/50))
            fprintf('.')
        end
    end
    fprintf('Done\n\n')
    
    %% subset to clear only
    
    subset = dataset.blurred==0;
    fulldata = fulldata(subset,:,:);
    fulllist = fulllist(subset,:);
    dataset = dataset(subset,:);

    %% make a proper triallist
    %add stimnum to dataset
    stimnumber = zeros(length(dataset.stimulus),1);
    exemplar = zeros(length(dataset.stimulus),1);
    for i=1:length(dataset.stimulus)
        stim = dataset.stimulus{i};
        stimnumber(i) = stimuli.numbers(find(strcmp(stimuli.names,stim),1));
        exemplar(i) = stimnumber(i)+(24*dataset.animate(i));
    end
    dataset = [dataset table(stimnumber,exemplar)]; %#ok<AGROW>
    %add triggeron and triggeroff to dataset
%     triggeron = triggertimes(triggerlabels==p.triggerstimon)';
%     dataset = [dataset table(triggeron)]; %#ok<AGROW>
%     triggeroff = triggertimes(triggerlabels==p.triggerstimoff)';
%     dataset = [dataset table(triggeroff)]; %#ok<AGROW>
    %add subject number
    [~,snum] = fileparts(subdir);
    dataset = [table(ones(height(dataset),1).*str2double(snum),'VariableNames',{'subjectnr'}) dataset]; %#ok<AGROW>
    %remove string vars, they are not possible in triallists
    d2 = dataset;
    d2.subject = [];
    d2.stimulus = [];
    for i=1:length(d2.Properties.VariableNames)
        VarMap{i}.column = d2.Properties.VariableNames{i}; %#ok<SAGROW>
    end
    TrialList = table2array(d2); 
    
    
    
    %% 9 resample
    for s=5%[1,5,10,20]
    
        p.subsamp = s; % 200Hz (1000Hz/s)
        fprintf('9/9 resampling to %dHz\n',MEGheaders{1}.sample_rate/p.subsamp)
        if s==1
            cd = fulldata;
        else
            cd = zeros([size(fulldata,1) ceil(size(fulldata,2)/p.subsamp) size(fulldata,3)]);
            n = size(fulldata,3);
            t = size(fulldata,1);
            parfor i = 1: size(fulldata,1) % loop over the trials
                for j = 1:n
                    cd(i,:,j) = decimate(squeeze(fulldata(i,:,j)),p.subsamp);
                end
                if ~mod(i,500)
                    fprintf('%i\n',i)
                end
            end
        end
        
        fprintf('Done\n\n')

        %clear MEGdata;

        p.PCAVarianceThreshold=.99;

        %% 8 dimred
        methods = {'pca'};
        if s==5
            %methods = {'none','pca','ica'};
        end
        
        for m = 1:length(methods)
            p.dimredmethod = methods{m};
             fprintf('8/9 Applying %s...',p.dimredmethod) 
            dims = size(cd);
            X = permute(cd,[2 1 3]);
            X = reshape(X,[dims(1)*dims(2) dims(3)]);
            fprintf('.')


            data = {};
            switch lower(p.dimredmethod)
                case 'pca'
                    data.method = 'PCA';
                    [COEFF, SCORE, LATENT, TSQUARED] = pca(X);fprintf('.')
                    clear X;
                    retain_ncomp = find(diff(cumsum(LATENT)/sum((LATENT))<p.PCAVarianceThreshold));fprintf('.')
                    pca_dat = SCORE(:,1:retain_ncomp);fprintf('.')
                    clear SCORE;
                    class_dat = reshape(pca_dat,[dims(2) dims(1) retain_ncomp]);fprintf('.')
                    clear pca_dat;
                    class_dat = permute(class_dat,[2 1 3]);fprintf('.')
                    data.PCAVarianceThreshold = p.PCAVarianceThreshold;
                    data.COEFF = COEFF;
                    data.LATENT = LATENT;
                case 'ica'
                    data.method = 'ICA';
                    X = X./max(abs(X(:))); %scale 0-1 or fastica will complain
                    %note: retain the same number of components as pca did
                    [icasig, A, W]  = fastica(X','lastEig',retain_ncomp); %the rows of icasig contain the estimated independent components.
                    clear X;
                    class_dat = reshape(icasig',[dims(2) dims(1) size(icasig,1)]);
                    clear icasig;
                    class_dat = permute(class_dat,[2 1 3]);
                    data.PCAVarianceThreshold = p.PCAVarianceThreshold;
                    data.A = A;
                    data.W = W;
                case 'none'
                    data.method = 'None';
                    class_dat = cd;
            end
            fprintf('Done\n\n')

            data.class_dat = class_dat;
            data.sampleFreq = MEGheaders{1}.sample_rate/p.subsamp;
            data.timevect = p.eventlength(1):1000/data.sampleFreq:p.eventlength(2);
            data.TrialList = TrialList;
            data.FullLabels = TrialList;
            data.triggertimes = triggertimes;
            data.triggerlabels = triggerlabels;
            data.EventLength = p.eventlength;
            data.VarMap = VarMap;
            data.MEGheaders = MEGheaders;

            targetfileh = sprintf('%s_%dHz_%s.mat',targetfile(1:end-4),data.sampleFreq,p.dimredmethod);
            fprintf('File completed..\n\nSaving to %s\n\n',targetfileh);
            save(targetfileh,'p','data','dataset');
        end
    end
      
end %function
