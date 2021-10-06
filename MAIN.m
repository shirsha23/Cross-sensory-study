% Main code for cross-sensory study
% ----------------------------------
% This code contains several steps:
% 1) Analyzing raw data and generating psychometric curves
%    1.1) Fix stimulus types (so that vestibular is always '1' and visual 
%         is always '2'
% 2) Priors Model - extract beta coefficients
% 3) Expanded priors model - extract betas when the 3 priors are separated
% 4) Plot figures 

%% STEP 1
%  Analyzing raw data and generating psychometric curves

clear;clc;close all;
currDir = pwd;

condLabel={'ves-ves' 'vis-vis' 'vis-ves' 'ves-vis'};

global is_new_pfit %shir - for using the new pfit
is_new_pfit = 1;
HO = 1; %heading offset (m=1 is the headings for INstring. INnum, however, starts at 1)

for condInd=1:4 
    DATA_DIR = strcat(currDir,'\',condLabel{condInd});    
    %read in subject information
    [INnum, INstring] = xlsread(strcat(DATA_DIR,'\All input_',condLabel{condInd},'.xlsx'),'subjects');
    [~,IX] = sort(INnum(:,1)); %in order of ascending subj_NUM
    subj_NUM  = INnum(IX,1);
    subj_ID = INnum(IX,2);
    subj_AGE = INnum(IX,3);
    subj_SEX = INstring(IX+HO,4);
    subj_EXCL = INstring(IX+HO,5);
    subj_COMM = INstring(IX+HO,6);
    
    % Read in session information
    [INnum, INstring] = xlsread(strcat(DATA_DIR,'\\All input_',condLabel{condInd},'.xlsx'),'sessions');
    sess_n = size(INnum,1);
    i=0;
    % Get session info
    for sess_i=1:sess_n
        i = i+1;
        disp(sprintf('Processing session %u',sess_i));
        filename = INstring{sess_i+HO,4};
        % Reading data file
        disp(sprintf('Reading file - %s',filename));
        filename_full_in=sprintf('%s\\Raw\\%s.mat',DATA_DIR,filename);
        filename_full_out=sprintf('%s\\Psychometric curves\\%s.mat',DATA_DIR,filename);
        data = importdata(filename_full_in);
        if isfield(data, 'SavedInfo') %sometimes the data is saved in a field called "SavedInfo" - I don't know why it sometimes is/not?
            data=data.SavedInfo;
        end
        exclude_sess = INstring{sess_i+HO,5};
        exclude_subj = subj_EXCL{subj_NUM==INnum(sess_i,1)};
        exclude = ~isempty(exclude_sess) | ~isempty(exclude_subj);
        % general properties
        sess_subjNUM(i) = INnum(sess_i,1);
        session_date = INstring{sess_i+HO,2};
        sess_date(i) = str2num(strrep((session_date(1:10)),'_',''));
        sess_AGE(i)=subj_AGE(subj_NUM==sess_subjNUM(i));
        sess_SEX(i)=subj_SEX(subj_NUM==sess_subjNUM(i));
        sess_coherence(i)=100;
        sess_exclude(i) = exclude;
        if (~sess_exclude(i))
            % analyze data and create psychometric curves
            AnaData{i} = analyze_single(filename_full_in,filename_full_out,data);
            % Arrange all data in one continuous record - shir 05/2018
            % it will include: trialCount, stimType, dir, response, coherence
            dataRec{i}=oneDataRec(data);
        end
    end
    % Save all data structures
    save(strcat(DATA_DIR,'\\All\Results'),'sess_*','subj_*','dataRec*','AnaData*');
end

%% STEP 1.1
% Fix stim types
%
% This code changes the stimType in the priors experiments:
% 1. In "vis-vis" and "vis-ves" conditions the vis is accidently marked as
%    stim type "1" (ves) so it is changed here to "2" (vis)
% 2. In all conditions stimTypes 6-10 are changed to 1 and 2

clear;clc;close all;
currDir = pwd;

condLabel={'vis-vis' 'ves-ves' 'ves-vis' 'vis-ves'};

for cond=1:4
    DATA_DIR=strcat(currDir,'\',condLabel{cond});
    ALL_DIR = strcat(DATA_DIR,'\All');
    load(strcat(ALL_DIR,'\Results'));
    
    dataRecWithOriginalStimType = dataRec;
    
    if cond==1|cond==4
        for i=1:length(dataRec)
            if ~isempty(dataRec{i})
                % visual trials are marked by mistake as "1" -> change to "2"
                dataRec{i}.stimType(dataRec{i}.stimType==1)=2;
                % change other stimTypes (6-10) to stadard (1-2)
                dataRec{i}.stimType(dataRec{i}.stimType==6)=1;
                dataRec{i}.stimType(dataRec{i}.stimType==9)=1;
                dataRec{i}.stimType(dataRec{i}.stimType==7)=2;
                dataRec{i}.stimType(dataRec{i}.stimType==10)=2;
            end
        end
    elseif cond==2|cond==3
        for i=1:length(dataRec)
            if ~isempty(dataRec{i})
                % change other stimTypes (6-10) to stadard (1-2)
                dataRec{i}.stimType(dataRec{i}.stimType==6)=1;
                dataRec{i}.stimType(dataRec{i}.stimType==9)=1;
                dataRec{i}.stimType(dataRec{i}.stimType==7)=2;
                dataRec{i}.stimType(dataRec{i}.stimType==10)=2;
            end
        end
    end
    
    save(strcat(DATA_DIR,'\\All\Results'),'sess_*','subj_*','dataRec','AnaData');
end

%% STEP 2

% Priors Model - extract beta coefficients
%---------------------------------------------------------------------
clear;clc;close all;
currDir = pwd;
%---------------------------------------------------------------------
% my variables
%---------------------------------------------------------------------
condLabel={'vis-vis' 'ves-ves' 'ves-vis' 'vis-ves'};
stimuliLabels={'Vestibular' 'Visual' 'Comb'};
coh=[100];
avgChoices = 0; % 1 - prevC is the average of 3 choices, 0 - mode
%---------------------------------------------------------------------
% calculate RMS for normalizing HD - load data from all conditions
allHD=[];
for condInd=1:4
    load(strcat(currDir,'\',condLabel{condInd},'\All\Results'));
    for i=1:length(dataRec)
        if ~isempty(dataRec{i})
            allHD=[allHD dataRec{i}.dir];
        end
    end
end
% calc RMS
rmsHD=rms(allHD);
%----------------------------------------------------------------------------------------------------
for condInd=1:4
    load(strcat(currDir,'\',condLabel{condInd},'\All\Results'));
    switch condInd
        case 1 %vis-vis
            stimuli=[2]; % 1-vestibular, 2-visual, 3-combined
        case 2 %ves-ves
            stimuli=[1]; % 1-vestibular, 2-visual, 3-combined
        case 3 %ves-vis
            stimuli=[2]; % 1-vestibular, 2-visual, 3-combined
        case 4 %vis-ves
            stimuli=[1]; % 1-vestibular, 2-visual, 3-combined
        otherwise
            warning('CONDITION NUMBER NOT DEFINED');
    end
    %----------------------------------------------------------
    % getting the subjects in the model
    display('Subjects in the model:');
    m=0; % counter for subjects in the model
    for sub=1:length(subj_NUM)
        currSum=sum(sess_subjNUM==subj_NUM(sub)&~sess_exclude);
        if currSum>0 % the subjects had some sessions
            disp(subj_NUM(sub));
            m=m+1;
            model_subjNUM(m)=subj_NUM(sub);
            model_subjSEX(m)=subj_SEX(sub);
        end
    end
    
    %% *** Priors model ***
    
    % getting beta coefficients for each subject in the model,
    % on each coherence, on each stimulus type (in this experiment
    % there's one coherence and one stimulus type)
    for s=1:length(stimuli)
        for c=1:length(coh)
            for m=1:length(model_subjNUM)
                subInd=find(sess_subjNUM==model_subjNUM(m)&sess_coherence==coh(c)&~sess_exclude);
                if ~isempty(subInd)
                    for currSubInd=1:length(subInd) % when there are repeating coherences sessions, subInd is >1
                        
                        % Using all responses as data points, we need these vectors:
                        dirVec=dataRec{subInd(currSubInd)}.dir; % vector of HD
                        dirVecNormRMS=dirVec./rmsHD; % normalize HD vector
                        respVec=dataRec{subInd(currSubInd)}.response; % vector of participant's responses
                        stimVec=dataRec{subInd(currSubInd)}.stimType; % vector of stimulus types
                        isTestTrial=dataRec{subInd(currSubInd)}.isTestTrial; % indicates which trials are test (1) or prior (0)
                        
                        % because the first 20 trials are "test" with no priors (for convergence)
                        % we need to mark the interesting test trials:
                        intTestTrials=zeros(length(isTestTrial),1);
                        % every test trial that is preceded by 3 priors is an
                        % interesting test trial
                        for t=4:length(intTestTrials)
                            % check if preceded by 3 priors
                            if isTestTrial(t)==1&(isTestTrial(t-1)+isTestTrial(t-2)+isTestTrial(t-3)==0)
                                intTestTrials(t)=1;
                            end
                        end
                        
                        % shir - check if the priors are from the correct stimType
                        % (there were some mistakes when running the exp)
                        correctPriorType=0;
                        priorsInd=find(isTestTrial==0);
                        if condInd==1|condInd==4
                            % priors supposed to be visual
                            if stimVec(priorsInd(1))==2
                                correctPriorType=1;
                            end
                        elseif condInd==2|condInd==3
                            % priors supposed to be vestibular
                            if stimVec(priorsInd(1))==1
                                correctPriorType=1;
                            end
                        end
                        
                        if sum(intTestTrials)>0&correctPriorType;
                            respVec01=respVec;
                            respVec01(respVec01==1)=-1; % left will be -1
                            respVec01(respVec01==2)=1; % right will be 1
                            
                            data_for_modelling.subjNum = [];
                            data_for_modelling.stimType = [];
                            data_for_modelling.trialNum = [];
                            data_for_modelling.currChoice = [];
                            data_for_modelling.currStim = [];
                            data_for_modelling.prevChoice = [];
                            data_for_modelling.prevChoice = [];
                            data_for_modelling.prevStim = [];
                            
                            % model test trials only
                            for t=1:length(intTestTrials)
                                if intTestTrials(t)==1
                                    data_for_modelling.subjNum = [data_for_modelling.subjNum model_subjNUM(m)];
                                    data_for_modelling.stimType = [data_for_modelling.stimType stimVec(t)];
                                    data_for_modelling.trialNum = [data_for_modelling.trialNum t];
                                    data_for_modelling.currChoice = [data_for_modelling.currChoice respVec01(t)];
                                    data_for_modelling.currStim = [data_for_modelling.currStim dirVecNormRMS(t)];
                                    if avgChoices
                                        avgPrevC = (respVec01(t-1)+respVec01(t-2)+respVec01(t-3))/3;
                                        data_for_modelling.prevChoice=[data_for_modelling.prevChoice avgPrevC];
                                    else
                                        if respVec01(t-1)+respVec01(t-2)+respVec01(t-3) < 0 % previous 3 trials were left
                                            data_for_modelling.prevChoice=[data_for_modelling.prevChoice -1];
                                        else % previous 3 trials were right
                                            data_for_modelling.prevChoice=[data_for_modelling.prevChoice 1];
                                        end
                                    end
                                    avg3prevStim=(dirVecNormRMS(t-1)+dirVecNormRMS(t-2)+dirVecNormRMS(t-3))/3;
                                    data_for_modelling.prevStim = [data_for_modelling.prevStim avg3prevStim];
                                end
                            end
                            data_for_modelling.currChoice(data_for_modelling.currChoice==-1)=0; %dependent variable is 0/1
                            % LOGISTIC REGRESSION MODEL
                            [betas,dev,stats] = glmfit([data_for_modelling.currStim' data_for_modelling.prevChoice'...
                                data_for_modelling.prevStim'],data_for_modelling.currChoice','binomial');
                            % use the new matlab function - fitglm - that also gives the AIC and BIC
                            % for my model (m3) : beta0 + betaStim + betaPrevChoice + betaPrevStim
                            mdl3mat=[data_for_modelling.currStim' data_for_modelling.prevChoice' data_for_modelling.prevStim'];
                            mdl3res = fitglm(mdl3mat,data_for_modelling.currChoice','Distribution','binomial','link','logit');
                            % and for the basic model (m0) : beta0 + betaStim
                            mdl0Mat=[data_for_modelling.currStim'];
                            mdl0res = fitglm(mdl0Mat,data_for_modelling.currChoice','Distribution','binomial','link','logit');
                            % and for another model (m1) : beta0 + betaStim + betaPrevStim
                            mdl1Mat=[data_for_modelling.currStim' data_for_modelling.prevStim'];
                            mdl1res = fitglm(mdl1Mat,data_for_modelling.currChoice','Distribution','binomial','link','logit');
                            % and for another model (m2) : beta0 + betaStim + betaPrevChoice
                            mdl2Mat=[data_for_modelling.currStim' data_for_modelling.prevChoice'];
                            mdl2res = fitglm(mdl2Mat,data_for_modelling.currChoice','Distribution','binomial','link','logit');
                            
                            model_m3BIC(m,c+currSubInd-1,s)=mdl3res.ModelCriterion.BIC;
                            model_m0BIC(m,c+currSubInd-1,s)=mdl0res.ModelCriterion.BIC;
                            model_m1BIC(m,c+currSubInd-1,s)=mdl1res.ModelCriterion.BIC;
                            model_m2BIC(m,c+currSubInd-1,s)=mdl2res.ModelCriterion.BIC;
                            model_m3AIC(m,c+currSubInd-1,s)=mdl3res.ModelCriterion.AIC;
                            model_m0AIC(m,c+currSubInd-1,s)=mdl0res.ModelCriterion.AIC;
                            model_m1AIC(m,c+currSubInd-1,s)=mdl1res.ModelCriterion.AIC;
                            model_m2AIC(m,c+currSubInd-1,s)=mdl2res.ModelCriterion.AIC;
                            model_stats(m,c+currSubInd-1,s)=stats;
                            model_beta0(m,c+currSubInd-1,s)=betas(1);
                            model_beta0PV(m,c+currSubInd-1,s)=stats.p(1);
                            model_betaStim(m,c+currSubInd-1,s)=betas(2);
                            model_betaStimPV(m,c+currSubInd-1,s)=stats.p(2);
                            model_betaPrevC(m,c+currSubInd-1,s)=betas(3);
                            model_betaPrevCPV(m,c+currSubInd-1,s)=stats.p(3);
                            if condInd==1|condInd==4 % only prev vis is relevant
                                model_betaPrevStimVes(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVesPV(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVis(m,c+currSubInd-1,s)=betas(4);
                                model_betaPrevStimVisPV(m,c+currSubInd-1,s)=stats.p(4);
                            elseif condInd==2|condInd==3 % only prev ves is relevant
                                model_betaPrevStimVes(m,c+currSubInd-1,s)=betas(4);
                                model_betaPrevStimVesPV(m,c+currSubInd-1,s)=stats.p(4);
                                model_betaPrevStimVis(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVisPV(m,c+currSubInd-1,s)=NaN;
                            end
                        else
                            model_m3BIC(m,c+currSubInd-1,s)=NaN;
                            model_m0BIC(m,c+currSubInd-1,s)=NaN;
                            model_m1BIC(m,c+currSubInd-1,s)=NaN;
                            model_m2BIC(m,c+currSubInd-1,s)=NaN;
                            model_beta0(m,c+currSubInd-1,s)=NaN;
                            model_m3AIC(m,c+currSubInd-1,s)=NaN;
                            model_m0AIC(m,c+currSubInd-1,s)=NaN;
                            model_m1AIC(m,c+currSubInd-1,s)=NaN;
                            model_m2AIC(m,c+currSubInd-1,s)=NaN;
                            model_beta0(m,c+currSubInd-1,s)=NaN;
                            model_betaStim(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevC(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevStimVes(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevStimVis(m,c+currSubInd-1,s)=NaN;
                        end
                    end
                else
                    model_m3BIC(m,c,s)=NaN;
                    model_m0BIC(m,c,s)=NaN;
                    model_m1BIC(m,c,s)=NaN;
                    model_m2BIC(m,c,s)=NaN;
                    model_m3AIC(m,c,s)=NaN;
                    model_m0AIC(m,c,s)=NaN;
                    model_m1AIC(m,c,s)=NaN;
                    model_m2AIC(m,c,s)=NaN;
                    model_beta0(m,c,s)=NaN;
                    model_betaStim(m,c,s)=NaN;
                    model_betaPrevC(m,c,s)=NaN;
                    model_betaPrevStimVes(m,c,s)=NaN;
                    model_betaPrevStimVis(m,c,s)=NaN;
                end
            end
        end
    end
    
    % cancel zeros that matlab created when expanding the matrix
    model_beta0(model_beta0==0)=NaN;
    model_betaStim(model_betaStim==0)=NaN;
    model_PrevC(model_betaPrevC==0)=NaN;
    model_PrevStimVes(model_betaPrevStimVes==0)=NaN;
    model_PrevStimVis(model_betaPrevStimVis==0)=NaN;
    
    if avgChoices
        save(strcat(currDir,'\',condLabel{condInd},'\Priors model\avg3choices\ModelResults_avg3priors_avgChoices+RMS'),'model_*');
    else
        save(strcat(currDir,'\',condLabel{condInd},'\Priors model\ModelResults_avg3priors+RMS'),'model_*');
    end
end

%% STEP 3
% Priors Model - extract beta coefficients
%                but here we look at 3 separate priors!
%---------------------------------------------------------------------
clear;clc;close all;
currDir = pwd;
%---------------------------------------------------------------------
% my variables
%---------------------------------------------------------------------
avgChoices = 0; % 1 - prevC is the average of 3 choices, 0 - mode
condLabel={'vis-vis' 'ves-ves' 'ves-vis' 'vis-ves'};
%---------------------------------------------------------------------
% calculate RMS for normalizing HD - load data from all conditions
allHD=[];
for condInd=1:4
    load(strcat(currDir,'\',condLabel{condInd},'\All\Results'));
    for i=1:length(dataRec)
        if ~isempty(dataRec{i})
            allHD=[allHD dataRec{i}.dir];
        end
    end
end
% calc RMS
rmsHD=rms(allHD);
%----------------------------------------------------------------------------------------------------

for condInd=1:4
    load(strcat(currDir,'\',condLabel{condInd},'\All\Results'));
    
    switch condInd
        case 1 %vis-vis
            stimuli=[2]; % 1-vestibular, 2-visual, 3-combined
        case 2 %ves-ves
            stimuli=[1]; % 1-vestibular, 2-visual, 3-combined
        case 3 %ves-vis
            stimuli=[2]; % 1-vestibular, 2-visual, 3-combined
        case 4 %vis-ves
            stimuli=[1]; % 1-vestibular, 2-visual, 3-combined
        otherwise
            warning('CONDITION NUMBER NOT DEFINED');
    end
    coh=[100];
    %----------------------------------------------------------
    % getting the subjects in the model
    display('Subjects in the model:');
    m=0; % counter for subjects in the model
    for sub=1:length(subj_NUM)
        currSum=sum(sess_subjNUM==subj_NUM(sub)&~sess_exclude);
        if currSum>0 % the subjects had some sessions
            disp(subj_NUM(sub));
            m=m+1;
            model_subjNUM(m)=subj_NUM(sub);
            model_subjSEX(m)=subj_SEX(sub);
        end
    end
    
    %% *** Priors model ***
    
    % getting the paramteres for each subject in the model, on each coherence,
    % on each stimulus type
    for s=1:length(stimuli)
        for c=1:length(coh)
            for m=1:length(model_subjNUM)
                subInd=find(sess_subjNUM==model_subjNUM(m)&sess_coherence==coh(c)&~sess_exclude);
                if ~isempty(subInd)
                    for currSubInd=1:length(subInd) % when there are repeating coherences sessions, subInd is >1
                        
                        % Using all responses as data points, we need these vectors:
                        dirVec=dataRec{subInd(currSubInd)}.dir; % vector of HD
                        dirVecNormRMS=dirVec./rmsHD; % normalize HD vector
                        respVec=dataRec{subInd(currSubInd)}.response; % vector of participant's responses
                        stimVec=dataRec{subInd(currSubInd)}.stimType; % vector of stimulus types
                        isTestTrial=dataRec{subInd(currSubInd)}.isTestTrial; % indicates which trials are test (1) or prior (0)
                        
                        % because the first 20 trials are "test" with no priors (for convergence)
                        % we need to mark the interesting test trials (those with priors):
                        intTestTrials=zeros(length(isTestTrial),1);
                        % every test trial that is preceded by 3 priors is an
                        % interesting test trial
                        for t=4:length(intTestTrials)
                            % check if preceded by 3 priors
                            if isTestTrial(t)==1&(isTestTrial(t-1)+isTestTrial(t-2)+isTestTrial(t-3)==0)
                                intTestTrials(t)=1;
                            end
                        end
                        % shir - check if the priors are from the correct stimType
                        % (there were some mistakes when running the exp)
                        correctPriorType=0;
                        priorsInd=find(isTestTrial==0);
                        if condInd==1|condInd==4
                            % priors supposed to be visual
                            if stimVec(priorsInd(1))==2
                                correctPriorType=1;
                            end
                        elseif condInd==2|condInd==3
                            % priors supposed to be vestibular
                            if stimVec(priorsInd(1))==1
                                correctPriorType=1;
                            end
                        end
                        
                        if sum(intTestTrials)>0&correctPriorType;
                            respVec01=respVec;
                            respVec01(respVec01==1)=-1; % left will be -1
                            respVec01(respVec01==2)=1; % right will be 1
                            
                            data_for_modelling.subjNum = [];
                            data_for_modelling.stimType = [];
                            data_for_modelling.trialNum = [];
                            data_for_modelling.currChoice = [];
                            data_for_modelling.currStim = [];
                            data_for_modelling.prevChoice = [];
                            data_for_modelling.prevChoice = [];
                            data_for_modelling.prevStim1 = [];
                            data_for_modelling.prevStim2 = [];
                            data_for_modelling.prevStim3 = [];
                            
                            % model test trials only
                            for t=1:length(intTestTrials)
                                if intTestTrials(t)==1
                                    %data_for_modelling.subjNum = [data_for_modelling.subjNum model_subjNUM(m)];
                                    data_for_modelling.stimType = [data_for_modelling.stimType stimVec(t)];
                                    data_for_modelling.trialNum = [data_for_modelling.trialNum t];
                                    data_for_modelling.currChoice = [data_for_modelling.currChoice respVec01(t)];
                                    %%% shir - just for simulation
                                    %data_for_modelling.currChoice = [data_for_modelling.currChoice simulatedCurrChoices(t)];
                                    %%%
                                    data_for_modelling.currStim = [data_for_modelling.currStim dirVecNormRMS(t)];
                                    if avgChoices
                                        avgPrevC = (respVec01(t-1)+respVec01(t-2)+respVec01(t-3))/3;
                                        data_for_modelling.prevChoice=[data_for_modelling.prevChoice avgPrevC];
                                    else
                                        % previous 3 choices become 1 choice
                                        if respVec01(t-1)+respVec01(t-2)+respVec01(t-3) < 0 % previous 3 trials were left
                                            data_for_modelling.prevChoice=[data_for_modelling.prevChoice -1];
                                        else % previous 3 trials were right
                                            data_for_modelling.prevChoice=[data_for_modelling.prevChoice 1];
                                        end
                                    end
                                    % prior stim - taken from last 3 priors
                                    data_for_modelling.prevStim1 = [data_for_modelling.prevStim1 dirVecNormRMS(t-1)];
                                    data_for_modelling.prevStim2 = [data_for_modelling.prevStim2 dirVecNormRMS(t-2)];
                                    data_for_modelling.prevStim3 = [data_for_modelling.prevStim3 dirVecNormRMS(t-3)];
                                end
                            end
                            data_for_modelling.currChoice(data_for_modelling.currChoice==-1)=0; %dependent variable is 0/1
                            % LOGISTIC REGRESSION MODEL
                            [betas,dev,stats] = glmfit([data_for_modelling.currStim' data_for_modelling.prevChoice'...
                                data_for_modelling.prevStim1' data_for_modelling.prevStim2' data_for_modelling.prevStim3'],...
                                data_for_modelling.currChoice','binomial');
                            % use the new matlab function - fitglm - that also gives the BIC
                            % for my model (m1) : beta0 + betaStim + betaPrevChoice + betaPrevStim
                            mdl1mat=[data_for_modelling.currStim' data_for_modelling.prevChoice'...
                                data_for_modelling.prevStim1' data_for_modelling.prevStim2' data_for_modelling.prevStim3'];
                            mdl1res = fitglm(mdl1mat,data_for_modelling.currChoice','Distribution','binomial','link','logit');
                            % and for the basic model (m0) : beta0 + betaStim
                            mdl0Mat=[data_for_modelling.currStim'];
                            mdl0res = fitglm(mdl0Mat,data_for_modelling.currChoice','Distribution','binomial','link','logit');
                            % and for another model (m2) : beta0 + betaStim + betaPrevStim
                            mdl2Mat=[data_for_modelling.currStim' data_for_modelling.prevStim1' data_for_modelling.prevStim2'...
                                data_for_modelling.prevStim3'];
                            mdl2res = fitglm(mdl2Mat,data_for_modelling.currChoice','Distribution','binomial','link','logit');
                            % and for another model (m3) : beta0 + betaStim + betaPrevChoice
                            mdl3Mat=[data_for_modelling.currStim' data_for_modelling.prevChoice'];
                            mdl3res = fitglm(mdl3Mat,data_for_modelling.currChoice','Distribution','binomial','link','logit');
                            
                            model_m1BIC(m,c+currSubInd-1,s)=mdl1res.ModelCriterion.BIC;
                            model_m0BIC(m,c+currSubInd-1,s)=mdl0res.ModelCriterion.BIC;
                            model_m2BIC(m,c+currSubInd-1,s)=mdl2res.ModelCriterion.BIC;
                            model_m3BIC(m,c+currSubInd-1,s)=mdl3res.ModelCriterion.BIC;
                            model_stats(m,c+currSubInd-1,s)=stats;
                            model_beta0(m,c+currSubInd-1,s)=betas(1);
                            model_beta0PV(m,c+currSubInd-1,s)=stats.p(1);
                            model_betaStim(m,c+currSubInd-1,s)=betas(2);
                            model_betaStimPV(m,c+currSubInd-1,s)=stats.p(2);
                            model_betaPrevC(m,c+currSubInd-1,s)=betas(3);
                            model_betaPrevCPV(m,c+currSubInd-1,s)=stats.p(3);
                            if condInd==1|condInd==4 % only prev vis is relevant
                                model_betaPrevStimVes1(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVes2(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVes3(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVes1PV(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVes2PV(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVes3PV(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVis1(m,c+currSubInd-1,s)=betas(4);
                                model_betaPrevStimVis2(m,c+currSubInd-1,s)=betas(5);
                                model_betaPrevStimVis3(m,c+currSubInd-1,s)=betas(6);
                                model_betaPrevStimVis1PV(m,c+currSubInd-1,s)=stats.p(4);
                                model_betaPrevStimVis2PV(m,c+currSubInd-1,s)=stats.p(5);
                                model_betaPrevStimVis3PV(m,c+currSubInd-1,s)=stats.p(6);
                            elseif condInd==2|condInd==3 % only prev ves is relevant
                                model_betaPrevStimVes1(m,c+currSubInd-1,s)=betas(4);
                                model_betaPrevStimVes2(m,c+currSubInd-1,s)=betas(5);
                                model_betaPrevStimVes3(m,c+currSubInd-1,s)=betas(6);
                                model_betaPrevStimVes1PV(m,c+currSubInd-1,s)=stats.p(4);
                                model_betaPrevStimVes2PV(m,c+currSubInd-1,s)=stats.p(5);
                                model_betaPrevStimVes3PV(m,c+currSubInd-1,s)=stats.p(6);
                                model_betaPrevStimVis1(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVis2(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVis3(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVis1PV(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVis2PV(m,c+currSubInd-1,s)=NaN;
                                model_betaPrevStimVis3PV(m,c+currSubInd-1,s)=NaN;
                            end
                        else
                            model_m1BIC(m,c+currSubInd-1,s)=NaN;
                            model_m0BIC(m,c+currSubInd-1,s)=NaN;
                            model_m2BIC(m,c+currSubInd-1,s)=NaN;
                            model_m3BIC(m,c+currSubInd-1,s)=NaN;
                            model_beta0(m,c+currSubInd-1,s)=NaN;
                            model_betaStim(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevC(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevStimVes1(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevStimVes2(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevStimVes3(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevStimVis1(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevStimVis2(m,c+currSubInd-1,s)=NaN;
                            model_betaPrevStimVis3(m,c+currSubInd-1,s)=NaN;
                        end
                    end
                else
                    model_m1BIC(m,c,s)=NaN;
                    model_m0BIC(m,c,s)=NaN;
                    model_m2BIC(m,c,s)=NaN;
                    model_m3BIC(m,c,s)=NaN;
                    model_beta0(m,c,s)=NaN;
                    model_betaStim(m,c,s)=NaN;
                    model_betaPrevC(m,c,s)=NaN;
                    model_betaPrevStimVes1(m,c,s)=NaN;
                    model_betaPrevStimVes2(m,c,s)=NaN;
                    model_betaPrevStimVes3(m,c,s)=NaN;
                    model_betaPrevStimVis1(m,c,s)=NaN;
                    model_betaPrevStimVis2(m,c,s)=NaN;
                    model_betaPrevStimVis3(m,c,s)=NaN;
                end
            end
        end
    end
    
    % cancel zeros that matlab created when expanding the matrix
    model_beta0(model_beta0==0)=NaN;
    model_betaStim(model_betaStim==0)=NaN;
    model_PrevC(model_betaPrevC==0)=NaN;
    model_PrevStimVes1(model_betaPrevStimVes1==0)=NaN;
    model_PrevStimVes2(model_betaPrevStimVes2==0)=NaN;
    model_PrevStimVes3(model_betaPrevStimVes3==0)=NaN;
    model_PrevStimVis1(model_betaPrevStimVis1==0)=NaN;
    model_PrevStimVis2(model_betaPrevStimVis2==0)=NaN;
    model_PrevStimVis3(model_betaPrevStimVis3==0)=NaN;
    
    if avgChoices
        save(strcat(currDir,'\',condLabel{condInd},'\Priors model\ModelResults_3priors_avgChoices+RMS'),'model_*');
    else
        save(strcat(currDir,'\',condLabel{condInd},'\Priors model\ModelResults_3priors+RMS'),'model_*');
    end
end

%% STEP 4
% Plot figures for cross-sensory study
%---------------------------------------------------------------------
clear;clc;close all;
currDir = pwd;

AVG_CHOICES_DIR = strcat(currDir,'\Priors model\avg3choices');

figs=[1 2 3 4 5];
% figs options are:
%                   1) PSE differences - summary
%                   2) M3 model - Four beta coefficients - boxplots
%                   3) Model comparisons - AIC
%                   4) M3E - Six beta coefficients - boxplot (3 separate priors)
%                   5) model comparisons - BIC
% models:
% m0 - "no history" - beta0 + current stim
% m1 - beta0 + current stim + prev stim
% m2 - beta0 + current stim + prev choice
% m3 - beta0 + current stim + prev choice + prev stim
%---------------------------------------------------------------------
% model_betaX is sized (m,c,s) m=subjects, c=coherence, s=stimuli
%---------------------------------------------------------------------
% my variables - CHECK IF NEED TO CHANGE
avgChoices = 0; % 1 - prevC is the average of 3 choices, 0 - mode
highlimit=35;
lowlimit=-35;
%---------------------------------------------------------------------
% general variables
condLabel={'vis-vis' 'vis-ves' 'ves-vis' 'ves-ves'};
stimuliLabel={'Vestibular' 'Visual' 'Comb'};
priorCondVec=[1 2 3 4]; %number of different priors conditions (vis-vis etc.)
coh=[100];
priorCoh=[100];
stimuli=[1]; %all results are in one layer
group={''}; %no groups
cohIdx=1:length(coh);
stimuliIdx=1:length(stimuli);
fontsize=12;

%%
if ismember(1,figs) % PSE differences - summary
    biasDiff=nan(20,4);
    meanBiasDiff=[];
    pVec=[];
    for priorCond=1:4
        biasL=[];
        biasR=[];
        load(strcat(currDir,'\',condLabel{priorCond},...
            '\All\for PSEs\ResultsPriors_',condLabel{priorCond}));
        if priorCond==1|priorCond==3
            biasL=sess_biasVisualpriorL;
            biasR=sess_biasVisualpriorR;
        elseif priorCond==2|priorCond==4
            biasL=sess_biasVespriorL;
            biasR=sess_biasVespriorR;
        end
        biasL=biasL(~isnan(biasL));
        biasR=biasR(~isnan(biasR));
        % some structures are shorter than 20 subjects so..
        biasDiff(1:length(biasR),priorCond)=biasR-biasL;
        meanBiasDiff(priorCond)=mean(biasDiff(:,priorCond));
        semBiasDiff(priorCond)=std(biasDiff(:,priorCond))/sqrt(length(biasDiff(:,priorCond)));
        display(sprintf('ttest for PSE difference against 0 - %s:',condLabel{priorCond}));
        [h,pv,ci,stats] = ttest(biasDiff(:,priorCond))
        % effect size - cohen's d = tstat/sqrt(df+1)
        cohensD = stats.tstat/sqrt(stats.df+1)
        display(sprintf('t%u = %.2f, p = %.3f, Cohen''s d = %.2f, 95%% CI = [%.2f %.2f], t-test',...
            stats.df,stats.tstat,pv,cohensD,ci(1),ci(2)));
        pVec(priorCond)=pv;
        
    end
    
    linewidth=4;
    
    % violin plots
    f=figure; hold on;
    title('Priors Effect','fontsize',fontsize);
    set(gcf,'Units','normalized','color','w','PaperPositionMode','auto','NumberTitle','off',...
        'Name','PSE_violin');
    
    line([0 4.5],[0 0],'linestyle',':','color',[.5 .5 .5]);
    distributionPlot(biasDiff,'showMM',4);
    
    set(gca,'xtick',[1:1:4],'xticklabel',({condLabel{1},condLabel{2},condLabel{3},condLabel{4}}),'fontsize',fontsize);
    set(gca,'ytick',[-10:5:10]);
    ylabel('delta PSE (R-L) [deg]','fontsize',12);
    xlabel('Condition','fontsize',12);
    ylim([-10 10]);
    xlim([0.5 4.5]);
    set(f,'Units','normalized','color','w','PaperPositionMode','auto');
    set(f,'name','PSE violin','NumberTitle','off');
    
    filename_fig=strcat(currDir,'\Plots\PSEs_summary_violin');
    print(gcf,'-dtiff','-r300', filename_fig);
    saveas(gcf,strcat(filename_fig,'.eps'));
    saveas(gcf,strcat(filename_fig,'.tif'));
    saveas(gcf,strcat(filename_fig,'.fig'));
end

%%
if ismember(2,figs) % M3 - four beta values
    condLabel={'ves-ves' 'vis-vis' 'vis-ves' 'ves-vis'};
    cohIdx=1;
    for g=1:length(group)
        for c=cohIdx
            for priorCond=1:length(priorCondVec)
                MODEL_DIR=strcat(currDir,'\',condLabel{priorCond},...
                    '\Priors model');
                if avgChoices
                    load(strcat(MODEL_DIR,'\ModelResults_avg3priors_avgChoices+RMS'));
                else
                    load(strcat(MODEL_DIR,'\ModelResults_avg3priors+RMS'));
                end
                
                % indicate group members (gm)
                % all subjects count
                gm=find(model_subjNUM>0);
                % get parameters
                beta0(gm,priorCond)=model_beta0(gm,c,1);
                betaStim(gm,priorCond)=model_betaStim(gm,c,1);
                betaPrevC(gm,priorCond)=model_betaPrevC(gm,c,1);
                betaPrevStimVis=model_betaPrevStimVis(gm,c,1);
                betaPrevStimVes=model_betaPrevStimVes(gm,c,1);
                % previous stim depends on the cond
                if priorCond==2|priorCond==3 % previous trial was vis
                    betaPrevStim(gm,priorCond)=betaPrevStimVis;
                elseif priorCond==1|priorCond==4 % previous trial was ves
                    betaPrevStim(gm,priorCond)=betaPrevStimVes;
                end
            end
            
            % box plots for each beta
            for currBeta=1:4
                f=figure; hold on;
                switch currBeta
                    case 1
                        beta=beta0;
                        betaLabel='beta0';
                        set(gca,'ytick',[-5:2.5:5],'yticklabel',[-5:2.5:5],'fontsize',fontsize);
                        ylim([-5 5]);
                    case 2
                        beta=betaStim;
                        betaLabel='betaStim';
                        set(gca,'ytick',[-5 0 10 20 30],'yticklabel',[-5 0 10 20 30],'fontsize',fontsize);
                        ylim([-5 30]);
                    case 3
                        beta=betaPrevC;
                        betaLabel='betaPrevC';
                        set(gca,'ytick',[-5:2.5:5],'yticklabel',[-5:2.5:5],'fontsize',fontsize);
                        ylim([-5 5]);
                    case 4
                        beta=betaPrevStim;
                        betaLabel='betaPrevStim';
                        set(gca,'ytick',[-5:2.5:5],'yticklabel',[-5:2.5:5],'fontsize',fontsize);
                        ylim([-5 5]);
                end
                M3title = ['Model M3: ' betaLabel];
                title(M3title,'fontsize',fontsize);
                boxplot(beta,'symbol','o');
                line([0 4.5],[0 0],'linestyle',':','color',[.5 .5 .5]);
                set(gca,'xtick',[1:1:4],'xticklabel',(condLabel),'fontsize',fontsize);
                ylabel('beta Coefficients [A.U.]','fontsize',12);
                xlabel('Condition','fontsize',12);
                xlim([0.5 4.5]);
                if currBeta==2
                    ylim([-5 30]);
                else
                    ylim([-5 5]);
                end
                % plot significance asterisks for each condition
                for b=1:4 % b - cond
                    disp(sprintf('condition: %s',condLabel{b}));
                    disp(sprintf('one sample T-test for %s #%u against 0',betaLabel,b));
                    [h,pv,ci,stats] = ttest(beta(:,b))
                    % effect size - cohen's d = tstat/sqrt(df+1)
                    cohensD = stats.tstat/sqrt(stats.df+1)
                    display(sprintf('t%u = %.2f, p = %.3f, Cohen''s d = %.2f, 95%% CI = [%.2f %.2f], t-test',...
                        stats.df,stats.tstat,pv,cohensD,ci(1),ci(2)));
                    if h
                        if currBeta==2 %betaStim's values are very high
                            astPos=27;
                        elseif nanmean(beta(:,b))>0
                            astPos=nanmean(beta(:,b))+1;
                        else
                            astPos=2;
                        end
                        if pv<0.001
                            text(b-0.05,astPos,'***','fontsize',12);
                        elseif pv<0.01
                            text(b-0.05,astPos,'**','fontsize',12);
                        elseif pv<0.05
                            text(b,astPos,'*','fontsize',12);
                        end
                    end
                end
                
                set(f,'Units','normalized','color','w','PaperPositionMode','auto');
                set(f,'name',sprintf('%s_boxplot',betaLabel),'NumberTitle','off');
                if avgChoices==1
                    filename_fig=strcat(AVG_CHOICES_DIR,'\',betaLabel,'_boxplot');
                elseif avgChoices==0
                    filename_fig=strcat(currDir,'\Plots\M3_',betaLabel,'_boxplot');
                end
                
                print(gcf,'-dtiff','-r300', filename_fig);
                saveas(gcf,strcat(filename_fig,'.eps'));
                saveas(gcf,strcat(filename_fig,'.tif'));
                saveas(gcf,strcat(filename_fig,'.fig'));
            end
        end
    end
end

%%
if ismember(3,figs) % model comparisons - AIC
    %order of plotting from bottom to top
    modelNumber=[3 2 0];
    modelLabel={'prior choices + prior stimuli' 'prior choices' 'no history'};
    %order of plotting from bottom to top
    condLabelforModels={'ves-vis' 'vis-ves' 'vis-vis' 'ves-ves'};
    
    bfMat=nan(20,3,4); % 20 subjects, 3 models, 4 conditions
    meanMat=nan(3,4);
    semMat=nan(3,4);
    for m=1:length(modelNumber)
        for priorCond=1:length(priorCondVec)
            MODEL_DIR=strcat(currDir,'\',condLabelforModels{priorCond},...
                '\Priors model');
            load(strcat(MODEL_DIR,'\ModelResults_avg3priors+RMS'));
            
            switch modelNumber(m)
                case 0
                    currModelAIC = model_m0AIC;
                case 2
                    currModelAIC = model_m2AIC;
                case 3
                    currModelAIC = model_m3AIC;
            end
            
            deltaAIC=model_m1AIC-currModelAIC;
            BFc2=exp(deltaAIC/2); %BFc2 - "c" for current model, 2 for m2
            BF2c=1./BFc2;
            logBF2c=log10(BF2c);
            bfMat(:,m,priorCond)=logBF2c;
            meanMat(m,priorCond)=nanmean(logBF2c);
            semMat(m,priorCond)=nanstd(logBF2c)./sqrt(length(logBF2c));
            % for paper - values before log
            bfMat_nolog(:,m,priorCond)=BF2c;
            meanMat_nolog(m,priorCond)=nanmean(BF2c);
        end
    end
    f=figure;hold on;
    title('Model Comparisons - AIC','fontsize',12);
    posVec=[1 2 3 4;7 8 9 10;13 14 15 16]; %bars position on y axis
    %bars=barh(posVec,meanMat); % short and works but with weird colors..
    % so i will do it manually:
    % (plotting order from bottom to top)
    bvesvis=barh([meanMat(1,1),0,0,0,0,0,...
        meanMat(2,1),0,0,0,0,0,...
        meanMat(3,1),0,0,0],'hist');
    bvisves=barh([0,meanMat(1,2),0,0,0,0,...
        0,meanMat(2,2),0,0,0,0,...
        0,meanMat(3,2),0,0],'hist');
    bvisvis=barh([0,0,meanMat(1,3),0,0,0,...
        0,0,meanMat(2,3),0,0,0,...
        0,0,meanMat(3,3),0],'hist');
    bvesves=barh([0,0,0,meanMat(1,4),0,0,...
        0,0,0,meanMat(2,4),0,0,...
        0,0,0,meanMat(3,4)],'hist');
    
    set(bvesvis,'FaceColor',[.9 .9 .9]);
    set(bvisves,'FaceColor',[.7 .7 .7]);
    set(bvisvis,'FaceColor',[.4 .4 .4]);
    set(bvesves,'FaceColor',[.2 .2 .2]);
    
    for m=1:length(modelNumber)
        for priorCond=1:length(priorCondVec)
            % errorbars
            line([meanMat(m,priorCond)-semMat(m,priorCond) ...
                meanMat(m,priorCond)+semMat(m,priorCond)],...
                [posVec(m,priorCond) posVec(m,priorCond)],'color','k');
            % scatter subjects
            currpos=posVec(m,priorCond);
            posSubj=currpos*ones(length(bfMat(:,m,priorCond)),1);
            %scatter(bfMat(:,m,priorCond),posSubj);
        end
    end
    xlabel('BF','fontsize',12);
    ylabel('Alternative Model','fontsize',12);
    set(gca,'ytick',[2.5,8.5,14.5],'yticklabel',{modelLabel{1},modelLabel{2},modelLabel{3}});
    set(gca,'xtick',[-1:1:3],'xticklabel',[{'10^{-1}'},{'10^{0}'},{'10^{1}'},{'10^{2}'},{'10^{3}'}]);
    legend({'ves-vis' 'vis-ves' 'vis-vis' 'ves-ves'},'fontsize',12,'location','southeast');
    xlim([-1 3]);
    ylim([0 17]);
    set(f,'Units','normalized','color','w','PaperPositionMode','auto');
    set(f,'name','Models Comparisons - AIC','NumberTitle','off');
    
    filename_fig=strcat(currDir,'\Plots\Model Comparisons - AIC');
    print(f,'-dtiff','-r300', filename_fig);
    saveas(f,strcat(filename_fig,'.eps'));
    saveas(f,strcat(filename_fig,'.tif'));
    saveas(f,strcat(filename_fig,'.fig'));
end

%%
if ismember(4,figs) % M3E - six beta values - 3 separate priors
    condLabel={'ves-ves' 'vis-vis' 'vis-ves' 'ves-vis'};
    cohIdx=1;
    for g=1:length(group)
        for c=cohIdx
            for priorCond=1:length(priorCondVec)
                MODEL_DIR=strcat(currDir,'\',condLabel{priorCond},...
                    '\Priors model');
                if avgChoices
                    load(strcat(MODEL_DIR,'\ModelResults_3priors_avgChoices+RMS'));
                else
                    load(strcat(MODEL_DIR,'\ModelResults_3priors+RMS'));
                end
                
                % indicate group members (gm)
                gm=find(model_subjNUM>0);
                % get parameters
                beta0(gm,priorCond)=model_beta0(gm,c,1);
                betaStim(gm,priorCond)=model_betaStim(gm,c,1);
                betaPrevC(gm,priorCond)=model_betaPrevC(gm,c,1);
                betaPrevStimVis1=model_betaPrevStimVis1(gm,c,1);
                betaPrevStimVis2=model_betaPrevStimVis2(gm,c,1);
                betaPrevStimVis3=model_betaPrevStimVis3(gm,c,1);
                betaPrevStimVes1=model_betaPrevStimVes1(gm,c,1);
                betaPrevStimVes2=model_betaPrevStimVes2(gm,c,1);
                betaPrevStimVes3=model_betaPrevStimVes3(gm,c,1);
                % previous stim depends on the cond
                if priorCond==2|priorCond==3 % previous trial was vis
                    betaPrevStim1(gm,priorCond)=betaPrevStimVis1;
                    betaPrevStim2(gm,priorCond)=betaPrevStimVis2;
                    betaPrevStim3(gm,priorCond)=betaPrevStimVis3;
                elseif priorCond==1|priorCond==4 % previous trial was ves
                    betaPrevStim1(gm,priorCond)=betaPrevStimVes1;
                    betaPrevStim2(gm,priorCond)=betaPrevStimVes2;
                    betaPrevStim3(gm,priorCond)=betaPrevStimVes3;
                end
            end
            
            % box plots for each beta
            for currBeta=1:6
                f=figure; hold on;
                switch currBeta
                    case 1
                        beta=beta0;
                        betaLabel='beta0';
                        set(gca,'ytick',[-5:2.5:5],'yticklabel',[-5:2.5:5],'fontsize',fontsize);
                        ylim([-5 5]);
                    case 2
                        beta=betaStim;
                        betaLabel='betaStim';
                        set(gca,'ytick',[-5 0 10 20 30],'yticklabel',[-5 0 10 20 30],'fontsize',fontsize);
                        ylim([-5 30]);
                    case 3
                        beta=betaPrevC;
                        betaLabel='betaPrevC';
                        set(gca,'ytick',[-5:2.5:5],'yticklabel',[-5:2.5:5],'fontsize',fontsize);
                        ylim([-5 5]);
                    case 4
                        beta=betaPrevStim1;
                        betaLabel='betaPrevStim1';
                        set(gca,'ytick',[-5:2.5:5],'yticklabel',[-5:2.5:5],'fontsize',fontsize);
                        ylim([-5 5]);
                    case 5
                        beta=betaPrevStim2;
                        betaLabel='betaPrevStim2';
                        set(gca,'ytick',[-5:2.5:5],'yticklabel',[-5:2.5:5],'fontsize',fontsize);
                        ylim([-5 5]);
                    case 6
                        beta=betaPrevStim3;
                        betaLabel='betaPrevStim3';
                        set(gca,'ytick',[-5:2.5:5],'yticklabel',[-5:2.5:5],'fontsize',fontsize);
                        ylim([-5 5]);
                end
                M3Etitle = ['Model M3E: ' betaLabel];
                title(M3Etitle,'fontsize',fontsize);
                boxplot(beta,'symbol','o');
                line([0 4.5],[0 0],'linestyle',':','color',[.5 .5 .5]);
                set(gca,'xtick',[1:1:4],'xticklabel',(condLabel),'fontsize',fontsize);
                ylabel('beta Coefficients [A.U.]','fontsize',12);
                xlabel('Condition','fontsize',12);
                xlim([0.5 4.5]);
                if currBeta==2
                    ylim([-5 30]);
                else
                    ylim([-5 5]);
                end
                % plot significance asterisks for each condition
                for b=1:4 % b - cond
                    [h,pv,ci,stats] = ttest(beta(:,b))
                    % effect size - cohen's d = tstat/sqrt(df+1)
                    cohensD = stats.tstat/sqrt(stats.df+1)
                    disp(sprintf('condition: %s',condLabel{b}));
                    disp(sprintf('one sample T-test for %s #%u against 0',betaLabel,b));
                    display(sprintf('t%u = %.2f, p = %.3f, Cohen''s d = %.2f, 95%% CI = [%.2f %.2f], t-test',...
                        stats.df,stats.tstat,pv,cohensD,ci(1),ci(2)));
                    
                    if h
                        if currBeta==2 %betaStim's values are very high
                            astPos=27;
                        elseif nanmean(beta(:,b))>0
                            astPos=nanmean(beta(:,b))+1;
                        else
                            astPos=2;
                        end
                        if pv<0.001
                            text(b-0.05,astPos,'***','fontsize',12);
                        elseif pv<0.01
                            text(b-0.05,astPos,'**','fontsize',12);
                        elseif pv<0.05
                            text(b,astPos,'*','fontsize',12);
                        end
                    end
                end
                
                set(f,'Units','normalized','color','w','PaperPositionMode','auto');
                set(f,'name',sprintf('%s_boxplot',betaLabel),'NumberTitle','off');
                filename_fig=strcat(currDir,'\Plots\M3E plots\M3E_',betaLabel,'_boxplot');
                
                print(gcf,'-dtiff','-r300', filename_fig);
                saveas(gcf,strcat(filename_fig,'.eps'));
                saveas(gcf,strcat(filename_fig,'.tif'));
                saveas(gcf,strcat(filename_fig,'.fig'));
            end
        end
    end
end

%%
if ismember(5,figs) % model comparisons - BIC
    %order of plotting from bottom to top
    modelNumber=[3 2 0];
    modelLabel={'prior choices + prior stimuli' 'prior choices' 'no history'};
    %order of plotting from bottom to top
    condLabelforModels={'ves-vis' 'vis-ves' 'vis-vis' 'ves-ves'};
    
    bfMat=nan(20,3,4); % 20 subjects, 3 models, 4 conditions
    meanMat=nan(3,4);
    semMat=nan(3,4);
    for m=1:length(modelNumber)
        for priorCond=1:length(priorCondVec)
            MODEL_DIR=strcat(currDir,'\',condLabelforModels{priorCond},...
                '\Priors model');
            load(strcat(MODEL_DIR,'\ModelResults_avg3priors+RMS'));
            
            switch modelNumber(m)
                case 0
                    currModelBIC = model_m0BIC;
                case 2
                    currModelBIC = model_m2BIC;
                case 3
                    currModelBIC = model_m3BIC;
            end
            
            deltaBIC=model_m1BIC-currModelBIC;
            BFc2=exp(deltaBIC/2); %BFc2 - "c" for current model, 2 for m2
            BF2c=1./BFc2;
            logBF2c=log10(BF2c);
            bfMat(:,m,priorCond)=logBF2c;
            meanMat(m,priorCond)=nanmean(logBF2c);
            semMat(m,priorCond)=nanstd(logBF2c)./sqrt(length(logBF2c));
            % for paper - values before log
            bfMat_nolog(:,m,priorCond)=BF2c;
            meanMat_nolog(m,priorCond)=nanmean(BF2c);
        end
    end
    f=figure;hold on;
    title('Model Comparisons - BIC','fontsize',12);
    posVec=[1 2 3 4;7 8 9 10;13 14 15 16]; %bars position on y axis
    %bars=barh(posVec,meanMat); % short and works but with weird colors..
    % so i will do it manually:
    % (plotting order from bottom to top)
    bvesvis=barh([meanMat(1,1),0,0,0,0,0,...
        meanMat(2,1),0,0,0,0,0,...
        meanMat(3,1),0,0,0],'hist');
    bvisves=barh([0,meanMat(1,2),0,0,0,0,...
        0,meanMat(2,2),0,0,0,0,...
        0,meanMat(3,2),0,0],'hist');
    bvisvis=barh([0,0,meanMat(1,3),0,0,0,...
        0,0,meanMat(2,3),0,0,0,...
        0,0,meanMat(3,3),0],'hist');
    bvesves=barh([0,0,0,meanMat(1,4),0,0,...
        0,0,0,meanMat(2,4),0,0,...
        0,0,0,meanMat(3,4)],'hist');
    
    set(bvesvis,'FaceColor',[.9 .9 .9]);
    set(bvisves,'FaceColor',[.7 .7 .7]);
    set(bvisvis,'FaceColor',[.4 .4 .4]);
    set(bvesves,'FaceColor',[.2 .2 .2]);
    
    for m=1:length(modelNumber)
        for priorCond=1:length(priorCondVec)
            % errorbars
            line([meanMat(m,priorCond)-semMat(m,priorCond) ...
                meanMat(m,priorCond)+semMat(m,priorCond)],...
                [posVec(m,priorCond) posVec(m,priorCond)],'color','k');
            % scatter subjects
            currpos=posVec(m,priorCond);
            posSubj=currpos*ones(length(bfMat(:,m,priorCond)),1);
            %scatter(bfMat(:,m,priorCond),posSubj);
        end
    end
    xlabel('BF','fontsize',12);
    ylabel('Alternative Model','fontsize',12);
    set(gca,'ytick',[2.5,8.5,14.5],'yticklabel',{modelLabel{1},modelLabel{2},modelLabel{3}});
    set(gca,'xtick',[-1:1:3],'xticklabel',[{'10^{-1}'},{'10^{0}'},{'10^{1}'},{'10^{2}'},{'10^{3}'}]);
    legend({'ves-vis' 'vis-ves' 'vis-vis' 'ves-ves'},'fontsize',12,'location','southeast');
    xlim([-1 3]);
    ylim([0 17]);
    set(f,'Units','normalized','color','w','PaperPositionMode','auto');
    set(f,'name','Models Comparisons - BIC','NumberTitle','off');
    
    filename_fig=strcat(currDir,'\Plots\Model Comparisons - BIC');
    print(f,'-dtiff','-r300', filename_fig);
    saveas(f,strcat(filename_fig,'.eps'));
    saveas(f,strcat(filename_fig,'.tif'));
    saveas(f,strcat(filename_fig,'.fig'));
end