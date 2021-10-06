% Simulation for priors model with 3 priors
% -----------------------------------------
% This code contains several steps:
% 1) Simulating participants' choices (how? detailed below)
% 2) Fitting the model to simulated data
% 3) Plot scatters to see original betas vs. recovered betas

%% STEP 1 
%  Simulation
%
% What the simulation does:
% 1. assign betas (real subject's betas)
% 2. simulate subject's HD and stimulus type vectors
% 3. compute the subject's current choices
% 4. repeat sections 2-3 many times for each subject
% ==> put the choices back in the model and check if betas were recovered
%-------------------------------------------------------------------------

% number of trials in this experiment:
% - 35 batches of 4 trials (3 priors + 1 test) => 140
% - this is happening for R and for L => 280
% - plus 10 convergence trials at the beginning for R and for L => 300
% - this is one block. 2 blocks => 600
%-------------------------------------------------------------------------
clear;clc;close all;
currDir = pwd;
%----------------------------------------
% my variables

% reverse order of prevStim - to double check simulation
revPrevStimVec = [0 1];

% by default simulation is running for vis-vis
% you can change it here if you want:
cond=1; %1)vis-vis 2)ves-ves 3)ves-vis 4)vis-ves

if cond==1
    priorStim=2; %vis
    testStim=2;  %vis
elseif cond==2
    priorStim=1; %ves
    testStim=1;  %ves
elseif cond==3
    priorStim=1; %ves
    testStim=2;  %vis
elseif cond==4
    priorStim=2; %vis
    testStim=1;  %ves
else
    warning('*** You need to assign correct stim types for this condition ***')
end
%----------------------------------------

% constants
hd=[16 8 4 2 1 0.5 0.25];
condLabel={'vis-vis' 'ves-ves' 'ves-vis' 'vis-ves'};

MODEL_DIR = strcat(currDir,'\',condLabel{cond},'\Priors model');
SIM_DIR = strcat(MODEL_DIR,'\Simulation');

load(strcat(MODEL_DIR,'\ModelResults_3priors+RMS'));
modelRMS=5.8396; % taken from model - to unnormalize the betas here

%% Simulate many times for each participant
%  parameters are like in the real experiment: 300 trials X 2 blocks etc.

iterNum=100; % many iterations for each subject
trialsPerBlock=300;         % like in the exp.
blocksNum=2;                % like in the exp.
trialsForConv=20;           % like in the exp.
priorsMu=5;priorsSigma=2.5; % like in the exp.
subjectsNum=20;             % like in the exp.
subjectsNum=subjectsNum+1;  % to add a fictitious subject
subjects=cell(subjectsNum,iterNum);

% simulate once in regular order and once in reverse order (to double
% check)
for revInd=1:length(revPrevStimVec)
    revPrevStim = revPrevStimVec(revInd);
    for subInd=1:subjectsNum
        % take the real subject's betas as input
        if subInd<subjectsNum % first 20 subjects are real, get their betas
            in_beta0(subInd)         = model_beta0(subInd);
            in_betaStim(subInd)      = model_betaStim(subInd)/modelRMS; % unnormalize
            in_betaPrevC(subInd)     = model_betaPrevC(subInd);
            if revPrevStim==0 % regular order of prevStim
                in_betaPrevStim1(subInd) = model_betaPrevStimVis1(subInd)/modelRMS; % unnormalize
                in_betaPrevStim2(subInd) = model_betaPrevStimVis2(subInd)/modelRMS; % unnormalize
                in_betaPrevStim3(subInd) = model_betaPrevStimVis3(subInd)/modelRMS; % unnormalize
            elseif revPrevStim==1 % reverse order of prevStim
                in_betaPrevStim1(subInd) = model_betaPrevStimVis3(subInd)/modelRMS; % unnormalize
                in_betaPrevStim2(subInd) = model_betaPrevStimVis2(subInd)/modelRMS; % unnormalize
                in_betaPrevStim3(subInd) = model_betaPrevStimVis1(subInd)/modelRMS; % unnormalize
            end
        elseif subInd==subjectsNum % subject 21 is fictitious
            in_beta0(subInd)         = 1.5;
            in_betaStim(subInd)      = 0.9;
            in_betaPrevC(subInd)     = 0.3;
            if revPrevStim==0 % regular order of prevStim
                in_betaPrevStim1(subInd) = -0.2;
                in_betaPrevStim2(subInd) = -0.1;
                in_betaPrevStim3(subInd) = 0;
            elseif revPrevStim==1 % reverse order of prevStim
                in_betaPrevStim1(subInd) = 0;
                in_betaPrevStim2(subInd) = -0.1;
                in_betaPrevStim3(subInd) = -0.2;
            end
        end
        
        for iterInd=1:iterNum
            for blockInd=1:blocksNum
                testsAfterR=45; % num of test trials after R prior. this is
                testsAfterL=45; % for each block - like in the exp.
                for i=1:trialsForConv % first 20 test trials for convergence (10 each R/L)
                    p = binornd(1,0.5); % draw "afterR" or "afterL"
                    if (p==1 && testsAfterR>35) | testsAfterL==35 % afterR and not finished 10 yet
                        subjects{subInd,iterInd,blockInd}.priorsType(i)='R';
                        if testsAfterR==45 % first trial of this kind
                            currHD = hd(1);
                        else
                            indR=find(subjects{subInd,iterInd,blockInd}.priorsType=='R');
                            lastRind=indR(length(indR)-1);
                            prevHDind=find(hd==abs(subjects{subInd,iterInd,blockInd}.dirVec(lastRind)));
                            if subjects{subInd,iterInd,blockInd}.indCorrect(lastRind)==1 % if previous test was correct
                                p = binornd(1,0.3); % 30% chance to go to a harder HD
                                if p==1 % go to a harder one
                                    if prevHDind==length(hd) % if it's the hardest one - stay
                                        currHD = hd(prevHDind);
                                    else % if there's a harder one - go there
                                        currHD = hd(prevHDind+1);
                                    end
                                elseif p==0 % 70% chance - stay in this HD
                                    currHD = hd(prevHDind);
                                end
                            elseif subjects{subInd,iterInd,blockInd}.indCorrect(lastRind)==0 % if previous test was wrong
                                p = binornd(1,0.8); % 80% chance to go to an easier HD
                                if p==1 % go to an easier one
                                    if prevHDind==1 % if it's the easiest one - stay
                                        currHD = hd(prevHDind);
                                    else % if there's an easier one - go there
                                        currHD = hd(prevHDind-1);
                                    end
                                elseif p==0 % 20% chance - stay in this HD
                                    currHD = hd(prevHDind);
                                end
                            end
                        end
                        testsAfterR=testsAfterR-1; % to count not more than 45
                    elseif (p==0 && testsAfterL>35) | testsAfterR==35 % afterL and not finished 10 yet
                        subjects{subInd,iterInd,blockInd}.priorsType(i)='L';
                        if testsAfterL==45 % first trial of this kind
                            currHD = hd(1);
                        else
                            indL=find(subjects{subInd,iterInd,blockInd}.priorsType=='L');
                            lastLind=indL(length(indL)-1);
                            prevHDind=find(hd==abs(subjects{subInd,iterInd,blockInd}.dirVec(lastLind)));
                            if subjects{subInd,iterInd,blockInd}.indCorrect(lastLind)==1 % if previous test was correct
                                p = binornd(1,0.3); % 30% chance to go to a harder HD
                                if p==1 % go to a harder one
                                    if prevHDind==length(hd) % if it's the hardest one - stay
                                        currHD = hd(prevHDind);
                                    else % if there's a harder one - go there
                                        currHD = hd(prevHDind+1);
                                    end
                                elseif p==0 % 70% chance - stay in this HD
                                    currHD = hd(prevHDind);
                                end
                            elseif subjects{subInd,iterInd,blockInd}.indCorrect(lastLind)==0 % if previous test was wrong
                                p = binornd(1,0.8); % 80% chance to go to an easier HD
                                if p==1 % go to an easier one
                                    if prevHDind==1 % if it's the easiest one - stay
                                        currHD = hd(prevHDind);
                                    else % if there's an easier one - go there
                                        currHD = hd(prevHDind-1);
                                    end
                                elseif p==0 % 20% chance - stay in this HD
                                    currHD = hd(prevHDind);
                                end
                            end
                        end
                        testsAfterL=testsAfterL-1; % to count not more than 45
                    end
                    % 50% chance to make the HD left/right
                    p = binornd(1,0.5);
                    if p==1 % make HD to the right
                        subjects{subInd,iterInd,blockInd}.dirVec(i) = currHD;
                    elseif p==0 % make HD to the left
                        subjects{subInd,iterInd,blockInd}.dirVec(i) = -currHD;
                    end
                    
                    % --- response ---
                    currStim = subjects{subInd,iterInd,blockInd}.dirVec(i);
                    if i<=3
                        % compute response for first three trials based on
                        % betaStim only
                        prevC = 0;
                        currStim = subjects{subInd,iterInd,blockInd}.dirVec(i);
                        comb=exp(in_beta0(subInd)+in_betaStim(subInd)*currStim+...
                            in_betaPrevC(subInd)*prevC);
                        prob1=comb/(comb+1);
                        subjects{subInd,iterInd,blockInd}.simulatedProb1(i)=prob1;
                        % draw the choice from a binomial dist.
                        r = binornd(1,prob1);
                        if r==1 % right was chosen
                            subjects{subInd,iterInd,blockInd}.respVec(i)=1;
                        else
                            subjects{subInd,iterInd,blockInd}.respVec(i)=-1;
                        end
                    else
                        % compute response for all trials from trial 4
                        prevC = mode(subjects{subInd,iterInd,blockInd}.respVec(i-3:i-1));
                        comb=exp(in_beta0(subInd)+in_betaStim(subInd)*currStim+...
                            in_betaPrevC(subInd)*prevC);
                        prob1=comb/(comb+1);
                        subjects{subInd,iterInd,blockInd}.simulatedProb1(i)=prob1;
                        % draw the choice from a binomial dist.
                        r = binornd(1,prob1);
                        if r==1 % right was chosen
                            subjects{subInd,iterInd,blockInd}.respVec(i)=1;
                        else
                            subjects{subInd,iterInd,blockInd}.respVec(i)=-1;
                        end
                    end
                    % check if response for test trial was correct (for the staircase)
                    if subjects{subInd,iterInd,blockInd}.respVec(i)==sign(subjects{subInd,iterInd,blockInd}.dirVec(i))
                        subjects{subInd,iterInd,blockInd}.indCorrect(i)=1;
                    else
                        subjects{subInd,iterInd,blockInd}.indCorrect(i)=0;
                    end
                    
                    % --- stim type ---
                    subjects{subInd,iterInd,blockInd}.stimVec(i) = testStim;
                    
                    % --- isTestTrial ---
                    % isTestTrial - indication for test trials
                    subjects{subInd,iterInd,blockInd}.isTestTrial(i)=1;
                    % intTestTrials - the interesting test trials (preceeded by 3 priors)
                    subjects{subInd,iterInd,blockInd}.intTestTrials(i)=0;
                end
                
                % now simulate all trials from 21 to 300
                % take into account separate staircase for R and L
                startPriors=trialsForConv+1;
                for i=startPriors:4:trialsPerBlock
                    testInd=i+3; % i,i+1,i+2 are priors trials, i+3 is a test trial
                    p = binornd(1,0.5); % draw "afterR" or "afterL"
                    if (p==1 && testsAfterR>0) | testsAfterL==0 % afterR and not finished yet
                        subjects{subInd,iterInd,blockInd}.dirVec(i:i+2)=...
                            normrnd(priorsMu,priorsSigma,1,3); % draw 3 priors from right priors dist.
                        subjects{subInd,iterInd,blockInd}.priorsType(testInd)='R';
                        indR=find(subjects{subInd,iterInd,blockInd}.priorsType=='R');
                        lastRind=indR(length(indR)-1);
                        prevHDind=find(hd==abs(subjects{subInd,iterInd,blockInd}.dirVec(lastRind)));
                        if subjects{subInd,iterInd,blockInd}.indCorrect(lastRind)==1 % if previous test was correct
                            p = binornd(1,0.3); % 30% chance to go to a harder HD
                            if p==1 % go to a harder one
                                if prevHDind==length(hd) % if it's the hardest one - stay
                                    currHD = hd(prevHDind);
                                else % if there's a harder one - go there
                                    currHD = hd(prevHDind+1);
                                end
                            elseif p==0 % 70% chance - stay in this HD
                                currHD = hd(prevHDind);
                            end
                        elseif subjects{subInd,iterInd,blockInd}.indCorrect(lastRind)==0 % if previous test was wrong
                            p = binornd(1,0.8); % 80% chance to go to an easier HD
                            if p==1 % go to an easier one
                                if prevHDind==1 % if it's the easiest one - stay
                                    currHD = hd(prevHDind);
                                else % if there's an easier one - go there
                                    currHD = hd(prevHDind-1);
                                end
                            elseif p==0 % 20% chance - stay in this HD
                                currHD = hd(prevHDind);
                            end
                        end
                        testsAfterR=testsAfterR-1; % to count not more than 45
                    elseif (p==0 && testsAfterL>0) | testsAfterR==0 % afterL and not finished yet
                        subjects{subInd,iterInd,blockInd}.dirVec(i:i+2)=...
                            normrnd(-priorsMu,priorsSigma,1,3); % draw 3 priors from left priors dist.
                        subjects{subInd,iterInd,blockInd}.priorsType(testInd)='L';
                        indL=find(subjects{subInd,iterInd,blockInd}.priorsType=='L');
                        lastLind=indL(length(indL)-1);
                        prevHDind=find(hd==abs(subjects{subInd,iterInd,blockInd}.dirVec(lastLind)));
                        if subjects{subInd,iterInd,blockInd}.indCorrect(lastLind)==1 % if previous test was correct
                            p = binornd(1,0.3); % 30% chance to go to a harder HD
                            if p==1 % go to a harder one
                                if prevHDind==length(hd) % if it's the hardest one - stay
                                    currHD = hd(prevHDind);
                                else % if there's a harder one - go there
                                    currHD = hd(prevHDind+1);
                                end
                            elseif p==0 % 70% chance - stay in this HD
                                currHD = hd(prevHDind);
                            end
                        elseif subjects{subInd,iterInd,blockInd}.indCorrect(lastLind)==0 % if previous test was wrong
                            p = binornd(1,0.8); % 80% chance to go to an easier HD
                            if p==1 % go to an easier one
                                if prevHDind==1 % if it's the easiest one - stay
                                    currHD = hd(prevHDind);
                                else % if there's an easier one - go there
                                    currHD = hd(prevHDind-1);
                                end
                            elseif p==0 % 20% chance - stay in this HD
                                currHD = hd(prevHDind);
                            end
                        end
                        testsAfterL=testsAfterL-1; % to count not more than 45
                    end
                    % 50% chance to make the test HD left/right
                    p = binornd(1,0.5);
                    if p==1 % make HD to the right
                        subjects{subInd,iterInd,blockInd}.dirVec(testInd) = currHD;
                    elseif p==0 % make HD to the left
                        subjects{subInd,iterInd,blockInd}.dirVec(testInd) = -currHD;
                    end
                    
                    % --- response ---
                    for ii=i:testInd % simulate response for priors and test
                        currStim = subjects{subInd,iterInd,blockInd}.dirVec(ii);
                        prevStim1 = subjects{subInd,iterInd,blockInd}.dirVec(ii-1);
                        prevStim2 = subjects{subInd,iterInd,blockInd}.dirVec(ii-2);
                        prevStim3 = subjects{subInd,iterInd,blockInd}.dirVec(ii-3);
                        prevC     = mode(subjects{subInd,iterInd,blockInd}.respVec(ii-3:ii-1));
                        comb=exp(in_beta0(subInd)+in_betaStim(subInd)*currStim+...
                            in_betaPrevC(subInd)*prevC+...
                            in_betaPrevStim1(subInd)*prevStim1+...
                            in_betaPrevStim2(subInd)*prevStim2+...
                            in_betaPrevStim3(subInd)*prevStim3);
                        prob1=comb/(comb+1);
                        subjects{subInd,iterInd,blockInd}.simulatedProb1(ii)=prob1;
                        % draw the choice from a binomial dist.
                        r = binornd(1,prob1);
                        if r==1 % right was chosen
                            subjects{subInd,iterInd,blockInd}.respVec(ii)=1;
                        else
                            subjects{subInd,iterInd,blockInd}.respVec(ii)=-1;
                        end
                    end
                    % check if response for test trial was correct (for the staircase)
                    if subjects{subInd,iterInd,blockInd}.respVec(testInd)==sign(subjects{subInd,iterInd,blockInd}.dirVec(testInd))
                        subjects{subInd,iterInd,blockInd}.indCorrect(testInd)=1;
                    else
                        subjects{subInd,iterInd,blockInd}.indCorrect(testInd)=0;
                    end
                    
                    % --- stim type ---
                    subjects{subInd,iterInd,blockInd}.stimVec(i:i+2) = priorStim;
                    subjects{subInd,iterInd,blockInd}.stimVec(testInd) = testStim;
                    
                    % --- isTestTrial ---
                    % isTestTrial - indication for test trials
                    subjects{subInd,iterInd,blockInd}.isTestTrial(i:i+2)=0; % priors
                    subjects{subInd,iterInd,blockInd}.isTestTrial(testInd)=1; % test
                    % intTestTrials - the interesting test trials (preceeded by 3 priors)
                    subjects{subInd,iterInd,blockInd}.isTestTrial(i:i+2)=0; % priors
                    subjects{subInd,iterInd,blockInd}.intTestTrials(testInd)=1; % test
                end
            end
            % merge the two blocks to one struct
            subjectsMerged{subInd,iterInd}.priorsType = ...
                [subjects{subInd,iterInd,1}.priorsType subjects{subInd,iterInd,2}.priorsType];
            subjectsMerged{subInd,iterInd}.dirVec = ...
                [subjects{subInd,iterInd,1}.dirVec subjects{subInd,iterInd,2}.dirVec];
            subjectsMerged{subInd,iterInd}.simulatedProb1 = ...
                [subjects{subInd,iterInd,1}.simulatedProb1 subjects{subInd,iterInd,2}.simulatedProb1];
            subjectsMerged{subInd,iterInd}.respVec = ...
                [subjects{subInd,iterInd,1}.respVec subjects{subInd,iterInd,2}.respVec];
            subjectsMerged{subInd,iterInd}.indCorrect = ...
                [subjects{subInd,iterInd,1}.indCorrect subjects{subInd,iterInd,2}.indCorrect];
            subjectsMerged{subInd,iterInd}.stimVec = ...
                [subjects{subInd,iterInd,1}.stimVec subjects{subInd,iterInd,2}.stimVec];
            subjectsMerged{subInd,iterInd}.isTestTrial = ...
                [subjects{subInd,iterInd,1}.isTestTrial subjects{subInd,iterInd,2}.isTestTrial];
            subjectsMerged{subInd,iterInd}.intTestTrials = ...
                [subjects{subInd,iterInd,1}.intTestTrials subjects{subInd,iterInd,2}.intTestTrials];
        end
    end
    
    if revPrevStim==0
        save(strcat(SIM_DIR,'\\SimulatedData'),'subjects','subjectsMerged','in_*','hd');
    elseif revPrevStim==1
        save(strcat(SIM_DIR,'\\SimulatedData_rev'),'subjects','subjectsMerged','in_*','hd');
    end
end

%% STEP 2+3
%  Priors Model on simulated data - extract beta coefficients
%                                   (we look at 3 separate priors)
%  And plot scatters

% run the model to extract betas for each iteration for each subject
% then check whether the extracted betas similar to the input ones
%---------------------------------------------------------------------
clear;clc;close all;
currDir = pwd;
%---------------------------------------------------------------------
% my variables
%---------------------------------------------------------------------
% by default simulation is running for vis-vis
% if you changed it then you need to change also this:
cond=1; %1)vis-vis 2)ves-ves 3)ves-vis 4)vis-ves

revPrevStimVec = [0 1 2];
% 0) regular order
% 1) reverse order of prevStim 
% 2) original betas not reversed vs. recovered reversed 
sections=[1 2];
% 1 - run model to get betas
% 2 - plot scatters of simulated betas vs. Original betas
%---------------------------------------------------------------------

condLabel={'vis-vis' 'ves-ves' 'ves-vis' 'vis-ves'};

MODEL_DIR = strcat(currDir,'\',condLabel{cond},'\Priors model');
SIM_DIR = strcat(MODEL_DIR,'\Simulation');

fontsize=12;
fontaxis=16;
linewidth=2;
circsize=300;

for revInd=1:length(revPrevStimVec)
    revPrevStim = revPrevStimVec(revInd);
    
    if revPrevStim==0 % regular order
        load(strcat(SIM_DIR,'\simulatedData'));
    elseif revPrevStim==1 % reverse order
        load(strcat(SIM_DIR,'\simulatedData_rev'));
    elseif revPrevStim==2 % regular betas vs. recovered reverse betas
        load(strcat(SIM_DIR,'\simulatedData'));
    end
    
    subjectsNum=size(subjectsMerged,1);
    iterNum=size(subjectsMerged,2);
    
    %% *** Priors model ***
    if ismember(1,sections)
        % getting the paramteres for each subject in the model
        for s=1:subjectsNum
            for i=1:iterNum
                if sum(subjectsMerged{s,i}.intTestTrials)>0
                    data_for_modelling.stimType = [];
                    data_for_modelling.trialNum = [];
                    data_for_modelling.currChoice = [];
                    data_for_modelling.currStim = [];
                    data_for_modelling.prevChoice = [];
                    data_for_modelling.prevChoice = [];
                    data_for_modelling.prevStim1 = [];
                    data_for_modelling.prevStim2 = [];
                    data_for_modelling.prevStim3 = [];
                    
                    % model fit for test trials only
                    for t=1:length(subjectsMerged{s,i}.intTestTrials)
                        if subjectsMerged{s,i}.intTestTrials(t)==1
                            data_for_modelling.stimType = [data_for_modelling.stimType subjectsMerged{s,i}.stimVec(t)];
                            data_for_modelling.trialNum = [data_for_modelling.trialNum t];
                            data_for_modelling.currChoice = [data_for_modelling.currChoice subjectsMerged{s,i}.respVec(t)];
                            data_for_modelling.currStim = [data_for_modelling.currStim subjectsMerged{s,i}.dirVec(t)];
                            % previous 3 choices become 1 choice
                            if subjectsMerged{s,i}.respVec(t-1)+subjectsMerged{s,i}.respVec(t-2)+subjectsMerged{s,i}.respVec(t-3) < 0 % previous 3 trials were left
                                data_for_modelling.prevChoice=[data_for_modelling.prevChoice -1];
                            else % previous 3 trials were right
                                data_for_modelling.prevChoice=[data_for_modelling.prevChoice 1];
                            end
                            % previous stim - taken from last 3 priors
                            data_for_modelling.prevStim1 = [data_for_modelling.prevStim1 subjectsMerged{s,i}.dirVec(t-1)];
                            data_for_modelling.prevStim2 = [data_for_modelling.prevStim2 subjectsMerged{s,i}.dirVec(t-2)];
                            data_for_modelling.prevStim3 = [data_for_modelling.prevStim3 subjectsMerged{s,i}.dirVec(t-3)];
                        end
                    end
                    data_for_modelling.currChoice(data_for_modelling.currChoice==-1)=0; %dependent variable is 0/1
                    % LOGISTIC REGRESSION MODEL
                    [betas,dev,stats] = glmfit([data_for_modelling.currStim' data_for_modelling.prevChoice'...
                        data_for_modelling.prevStim1' data_for_modelling.prevStim2' data_for_modelling.prevStim3'],...
                        data_for_modelling.currChoice','binomial');
                    
                    sim_stats(s,i)=stats;
                    sim_beta0(s,i)=betas(1);
                    sim_beta0PV(s,i)=stats.p(1);
                    sim_betaStim(s,i)=betas(2);
                    sim_betaStimPV(s,i)=stats.p(2);
                    sim_betaPrevC(s,i)=betas(3);
                    sim_betaPrevCPV(s,i)=stats.p(3);
                    if cond==1|cond==4 % only prev vis is relevant
                        sim_betaPrevStimVes1(s,i)=NaN;
                        sim_betaPrevStimVes2(s,i)=NaN;
                        sim_betaPrevStimVes3(s,i)=NaN;
                        sim_betaPrevStimVes1PV(s,i)=NaN;
                        sim_betaPrevStimVes2PV(s,i)=NaN;
                        sim_betaPrevStimVes3PV(s,i)=NaN;
                        sim_betaPrevStimVis1(s,i)=betas(4);
                        sim_betaPrevStimVis2(s,i)=betas(5);
                        sim_betaPrevStimVis3(s,i)=betas(6);
                        sim_betaPrevStimVis1PV(s,i)=stats.p(4);
                        sim_betaPrevStimVis2PV(s,i)=stats.p(5);
                        sim_betaPrevStimVis3PV(s,i)=stats.p(6);
                    elseif cond==2|cond==3 % only prev ves is relevant
                        sim_betaPrevStimVes1(s,i)=betas(4);
                        sim_betaPrevStimVes2(s,i)=betas(5);
                        sim_betaPrevStimVes3(s,i)=betas(6);
                        sim_betaPrevStimVes1PV(s,i)=stats.p(4);
                        sim_betaPrevStimVes2PV(s,i)=stats.p(5);
                        sim_betaPrevStimVes3PV(s,i)=stats.p(6);
                        sim_betaPrevStimVis1(s,i)=NaN;
                        sim_betaPrevStimVis2(s,i)=NaN;
                        sim_betaPrevStimVis3(s,i)=NaN;
                        sim_betaPrevStimVis1PV(s,i)=NaN;
                        sim_betaPrevStimVis2PV(s,i)=NaN;
                        sim_betaPrevStimVis3PV(s,i)=NaN;
                    end
                else
                    sim_m1BIC(s,i)=NaN;
                    sim_m0BIC(s,i)=NaN;
                    sim_m2BIC(s,i)=NaN;
                    sim_m3BIC(s,i)=NaN;
                    sim_beta0(s,i)=NaN;
                    sim_betaStim(s,i)=NaN;
                    sim_betaPrevC(s,i)=NaN;
                    sim_betaPrevStimVes1(s,i)=NaN;
                    sim_betaPrevStimVes2(s,i)=NaN;
                    sim_betaPrevStimVes3(s,i)=NaN;
                    sim_betaPrevStimVis1(s,i)=NaN;
                    sim_betaPrevStimVis2(s,i)=NaN;
                    sim_betaPrevStimVis3(s,i)=NaN;
                end
            end
        end
        
        % cancel zeros that matlab created when expanding the matrix
        sim_beta0(sim_beta0==0)=NaN;
        sim_betaStim(sim_betaStim==0)=NaN;
        sim_PrevC(sim_betaPrevC==0)=NaN;
        sim_PrevStimVes1(sim_betaPrevStimVes1==0)=NaN;
        sim_PrevStimVes2(sim_betaPrevStimVes2==0)=NaN;
        sim_PrevStimVes3(sim_betaPrevStimVes3==0)=NaN;
        sim_PrevStimVis1(sim_betaPrevStimVis1==0)=NaN;
        sim_PrevStimVis2(sim_betaPrevStimVis2==0)=NaN;
        sim_PrevStimVis3(sim_betaPrevStimVis3==0)=NaN;
        
        if revPrevStim==0
            save(strcat(SIM_DIR,'\\Betas for simulated data'),'sim_*');
        elseif revPrevStim==1
            save(strcat(SIM_DIR,'\\Betas for simulated data_rev'),'sim_*');
        end
    end
    %% scatter to check if betas from simulation = model betas
    if ismember(2,sections)
        
        if revPrevStim==0
            load(strcat(SIM_DIR,'\Betas for simulated data'));
            revString='';
        elseif revPrevStim==1
            load(strcat(SIM_DIR,'\Betas for simulated data_rev'));
            revString='_rev';
        elseif revPrevStim==2
            % here we want to plot the original regular order betas
            % vs the recovered reversed order betas
            
            % load reverse order to keep the recovered betas
            load(strcat(SIM_DIR,'\Betas for simulated data_rev'));
            revString='_origVSrev';
            % keep the recovered reversed betas in different name
            % (so it won't be replaced when loading the regular order recovered betas)
            sim_betaPrevStimVis1REV = sim_betaPrevStimVis1;
            sim_betaPrevStimVis2REV = sim_betaPrevStimVis2;
            sim_betaPrevStimVis3REV = sim_betaPrevStimVis3;
            sim_beta0REV = sim_beta0;
            sim_betaPrevCREV = sim_betaPrevC;
            sim_betaStimREV = sim_betaStim;
            % now load the regular order for the original betas
            load(strcat(SIM_DIR,'\Betas for simulated data'));
            % delete the regular order recovered betas, they are not needed
            clear sim_betaPrevStimVis1 sim_betaPrevStimVis2...
                sim_betaPrevStimVis3 sim_beta0 sim_betaPrevC sim_betaStim;
            % rename again the recovered betas (for the code that follows..)
            sim_betaPrevStimVis1 = sim_betaPrevStimVis1REV;
            sim_betaPrevStimVis2 = sim_betaPrevStimVis2REV;
            sim_betaPrevStimVis3 = sim_betaPrevStimVis3REV;
            sim_beta0 = sim_beta0REV;
            sim_betaPrevC = sim_betaPrevCREV;
            sim_betaStim = sim_betaStimREV;
        end
        
        % betaPrevStim1
        sim_betaPrevStimVis1mean=nanmean(sim_betaPrevStimVis1,2);
        sim_betaPrevStimVis1med=nanmedian(sim_betaPrevStimVis1,2);
        betaPrevStim1MAD=mad(sim_betaPrevStimVis1,1,2); % 1 for median absolute deviation, 2-dimension
        f=figure;hold on;
        errorbar(in_betaPrevStim1,sim_betaPrevStimVis1med,betaPrevStim1MAD,betaPrevStim1MAD,'o','color',[.5 .5 .5],'linewidth',linewidth);
        scatter(in_betaPrevStim1,sim_betaPrevStimVis1med,circsize,'markerfacecolor',[0 .8 .8],'markeredgecolor','k');
        scatter(in_betaPrevStim1(subjectsNum),sim_betaPrevStimVis1med(subjectsNum),circsize,'filled','red','markeredgecolor','k'); % fictitious subject
        title('\bf\it \beta \rm\bf previous stimulus (t-1)','fontsize',fontaxis);
        if revPrevStim==0 % regular order 
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==1 % both in reverse order
            xlabel('Original (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==2 % only recovered in reverse order
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        end
        xlim([-.6 .6]);ylim([-.6 .6]);
        set(gca,'xtick',[-.6:.3:.6],'xticklabel',[-.6:.3:.6],'fontsize',fontaxis);
        set(gca,'ytick',[-.6:.3:.6],'yticklabel',[-.6:.3:.6]);
        line([-.6 .6],[-.6 .6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([-.6 .6],[0 0],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([0 0],[-.6 .6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        [rho,pvr]=corr(sim_betaPrevStimVis1med,in_betaPrevStim1')
        text(-.5,.5,sprintf('R^{2} = %.2f',rho^2));%, pv = %.3f',rho,pvr));
        set(f,'Units','normalized','color','w','PaperPositionMode','auto');
        set(f,'name','betaPrevStim1','NumberTitle','off');
        filename_fig=strcat(SIM_DIR,'\scatter - betaPrevStim1',revString);
        print(f,'-dtiff','-r300', filename_fig);
        saveas(f,strcat(filename_fig,'.eps'));
        saveas(f,strcat(filename_fig,'.tif'));
        saveas(f,strcat(filename_fig,'.fig'));
        
        % betaPrevStim2
        sim_betaPrevStimVis2mean=nanmean(sim_betaPrevStimVis2,2);
        sim_betaPrevStimVis2med=nanmedian(sim_betaPrevStimVis2,2);
        betaPrevStim2MAD=mad(sim_betaPrevStimVis2,1,2); % 1 for median absolute deviation, 2-dimension
        f=figure;hold on;
        errorbar(in_betaPrevStim2,sim_betaPrevStimVis2med,betaPrevStim2MAD,betaPrevStim2MAD,'o','color',[.5 .5 .5],'linewidth',linewidth);
        scatter(in_betaPrevStim2,sim_betaPrevStimVis2med,circsize,'markerfacecolor',[0 .8 .8],'markeredgecolor','k');
        scatter(in_betaPrevStim2(subjectsNum),sim_betaPrevStimVis2med(subjectsNum),circsize,'filled','red','markeredgecolor','k'); % fictitious subject
        title('\bf\it \beta \rm\bf previous stimulus (t-2)','fontsize',fontaxis);
        if revPrevStim==0 % regular order
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==1 % both in reverse order
            xlabel('Original (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==2 % only recovered in reverse order
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        end
        xlim([-.6 .6]);ylim([-.6 .6]);
        set(gca,'xtick',[-.6:.3:.6],'xticklabel',[-.6:.3:.6],'fontsize',fontaxis);
        set(gca,'ytick',[-.6:.3:.6],'yticklabel',[-.6:.3:.6]);
        line([-.6 .6],[-.6 .6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([-.6 .6],[0 0],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([0 0],[-.6 .6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        [rho,pvr]=corr(sim_betaPrevStimVis2med,in_betaPrevStim2')
        text(-.5,.5,sprintf('R^{2} = %.2f',rho^2));%, pv = %.3f',rho,pvr));
        set(f,'Units','normalized','color','w','PaperPositionMode','auto');
        set(f,'name','betaPrevStim2','NumberTitle','off');
        filename_fig=strcat(SIM_DIR,'\scatter - betaPrevStim2',revString);
        print(f,'-dtiff','-r300', filename_fig);
        saveas(f,strcat(filename_fig,'.eps'));
        saveas(f,strcat(filename_fig,'.tif'));
        saveas(f,strcat(filename_fig,'.fig'));
        
        % betaPrevStim3
        sim_betaPrevStimVis3mean=nanmean(sim_betaPrevStimVis3,2);
        sim_betaPrevStimVis3med=nanmedian(sim_betaPrevStimVis3,2);
        betaPrevStim3MAD=mad(sim_betaPrevStimVis3,1,2); % 1 for median absolute deviation, 2-dimension
        f=figure;hold on;
        errorbar(in_betaPrevStim3,sim_betaPrevStimVis3med,betaPrevStim3MAD,betaPrevStim3MAD,'o','color',[.5 .5 .5],'linewidth',linewidth);
        scatter(in_betaPrevStim3,sim_betaPrevStimVis3med,circsize,'markerfacecolor',[0 .8 .8],'markeredgecolor','k');
        scatter(in_betaPrevStim3(subjectsNum),sim_betaPrevStimVis3med(subjectsNum),circsize,'filled','red','markeredgecolor','k'); % fictitious subject
        title('\bf\it \beta \rm\bf previous stimulus (t-3)','fontsize',fontaxis);
        if revPrevStim==0 % regular order 
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==1 % both in reverse order
            xlabel('Original (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==2 % only recovered in reverse order
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        end
        xlim([-.6 .6]);ylim([-.6 .6]);
        set(gca,'xtick',[-.6:.3:.6],'xticklabel',[-.6:.3:.6],'fontsize',fontaxis);
        set(gca,'ytick',[-.6:.3:.6],'yticklabel',[-.6:.3:.6]);
        line([-.6 .6],[-.6 .6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([-.6 .6],[0 0],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([0 0],[-.6 .6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        [rho,pvr]=corr(sim_betaPrevStimVis3med,in_betaPrevStim3')
        text(-.5,.5,sprintf('R^{2} = %.2f',rho^2));%, pv = %.3f',rho,pvr));
        set(f,'Units','normalized','color','w','PaperPositionMode','auto');
        set(f,'name','betaPrevStim3','NumberTitle','off');
        filename_fig=strcat(SIM_DIR,'\scatter - betaPrevStim3',revString);
        print(f,'-dtiff','-r300', filename_fig);
        saveas(f,strcat(filename_fig,'.eps'));
        saveas(f,strcat(filename_fig,'.tif'));
        saveas(f,strcat(filename_fig,'.fig'));
        
        % betaStim
        sim_betaStimmean=nanmean(sim_betaStim,2);
        sim_betaStimmed=nanmedian(sim_betaStim,2);
        betaStimMAD=mad(sim_betaStim,1,2); % 1 for median absolute deviation, 2-dimension
        f=figure;hold on;
        errorbar(in_betaStim,sim_betaStimmed,betaStimMAD,betaStimMAD,'o','color',[.5 .5 .5],'linewidth',linewidth);
        scatter(in_betaStim,sim_betaStimmed,circsize,'markerfacecolor',[0 .8 .8],'markeredgecolor','k');
        scatter(in_betaStim(subjectsNum),sim_betaStimmed(subjectsNum),circsize,'filled','red','markeredgecolor','k'); % fictitious subject
        title('\bf\it \beta \rm\bf current stimulus','fontsize',fontaxis);
        if revPrevStim==0 % regular order 
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==1 % both in reverse order
            xlabel('Original (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==2 % only recovered in reverse order
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        end
        xlim([-3 6]);ylim([-3 6]);
        set(gca,'xtick',[-3:3:6],'xticklabel',[-3:3:6],'fontsize',fontaxis);
        set(gca,'ytick',[-3:3:6],'yticklabel',[-3:3:6]);
        line([-3 6],[-3 6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([-3 6],[0 0],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([0 0],[-3 6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        [rho,pvr]=corr(sim_betaStimmed,in_betaStim')
        text(.1,4.5,sprintf('R^{2} = %.2f',rho^2));%, pv = %.3f',rho,pvr));
        set(f,'Units','normalized','color','w','PaperPositionMode','auto');
        set(f,'name','betaCurrStim','NumberTitle','off');
        filename_fig=strcat(SIM_DIR,'\scatter - betaCurrStim',revString);
        print(f,'-dtiff','-r300', filename_fig);
        saveas(f,strcat(filename_fig,'.eps'));
        saveas(f,strcat(filename_fig,'.tif'));
        saveas(f,strcat(filename_fig,'.fig'));
        
        % betaPrevC
        sim_betaPrevCmean=nanmean(sim_betaPrevC,2);
        sim_betaPrevCmed=nanmedian(sim_betaPrevC,2);
        betaPrevCMAD=mad(sim_betaPrevC,1,2); % 1 for median absolute deviation, 2-dimension
        f=figure;hold on;
        errorbar(in_betaPrevC,sim_betaPrevCmed,betaPrevCMAD,betaPrevCMAD,'o','color',[.5 .5 .5],'linewidth',linewidth);
        scatter(in_betaPrevC,sim_betaPrevCmed,circsize,'markerfacecolor',[0 .8 .8],'markeredgecolor','k');
        scatter(in_betaPrevC(subjectsNum),sim_betaPrevCmed(subjectsNum),circsize,'filled','red','markeredgecolor','k'); % fictitious subject
        title('\bf\it \beta \rm\bf previous choice','fontsize',fontaxis);
        if revPrevStim==0 % regular order 
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==1 % both in reverse order
            xlabel('Original (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==2 % only recovered in reverse order
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        end
        xlim([-3 6]);ylim([-3 6]);
        set(gca,'xtick',[-3:3:6],'xticklabel',[-3:3:6],'fontsize',fontaxis);
        set(gca,'ytick',[-3:3:6],'yticklabel',[-3:3:6]);
        line([-3 6],[-3 6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([-3 6],[0 0],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([0 0],[-3 6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        [rho,pvr]=corr(sim_betaPrevCmed,in_betaPrevC')
        text(-2,2,sprintf('R^{2} = %.2f',rho^2));%, pv = %.3f',rho,pvr));
        set(f,'Units','normalized','color','w','PaperPositionMode','auto');
        set(f,'name','betaPrevC','NumberTitle','off');
        filename_fig=strcat(SIM_DIR,'\scatter - betaPrevC',revString);
        print(f,'-dtiff','-r300', filename_fig);
        saveas(f,strcat(filename_fig,'.eps'));
        saveas(f,strcat(filename_fig,'.tif'));
        saveas(f,strcat(filename_fig,'.fig'));
        
        % beta0
        sim_beta0mean=nanmean(sim_beta0,2);
        sim_beta0med=nanmedian(sim_beta0,2);
        beta0MAD=mad(sim_beta0,1,2); % 1 for median absolute deviation, 2-dimension
        f=figure;hold on;
        errorbar(in_beta0,sim_beta0med,beta0MAD,beta0MAD,'o','color',[.5 .5 .5],'linewidth',linewidth);
        scatter(in_beta0,sim_beta0med,circsize,'markerfacecolor',[0 .8 .8],'markeredgecolor','k');
        scatter(in_beta0(subjectsNum),sim_beta0med(subjectsNum),circsize,'filled','red','markeredgecolor','k'); % fictitious subject
        title('\bf\it\beta\rm\bf0 - Bias','fontsize',fontaxis);
        if revPrevStim==0 % regular order 
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==1 % both in reverse order
            xlabel('Original (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        elseif revPrevStim==2 % only recovered in reverse order
            xlabel('Original (simulated) \it\beta \rm[A.U.]','fontsize',fontaxis);
            ylabel('Recovered (simulated in reverse order) \it\beta \rm[A.U.]','fontsize',fontaxis);
        end
        xlim([-3 6]);ylim([-3 6]);
        set(gca,'xtick',[-3:3:6],'xticklabel',[-3:3:6],'fontsize',fontaxis);
        set(gca,'ytick',[-3:3:6],'yticklabel',[-3:3:6]);
        line([-3 6],[-3 6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([-3 6],[0 0],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        line([0 0],[-3 6],'LineStyle',':','color',[.5 .5 .5],'linewidth',linewidth);
        [rho,pvr]=corr(sim_beta0med,in_beta0')
        text(-2,5,sprintf('R^{2} = %.2f',rho^2));%, pv = %.3f',rho,pvr));
        set(f,'Units','normalized','color','w','PaperPositionMode','auto');
        set(f,'name','beta0','NumberTitle','off');
        filename_fig=strcat(SIM_DIR,'\scatter - beta0',revString);
        print(f,'-dtiff','-r300', filename_fig);
        saveas(f,strcat(filename_fig,'.eps'));
        saveas(f,strcat(filename_fig,'.tif'));
        saveas(f,strcat(filename_fig,'.fig'));
        
    end
end