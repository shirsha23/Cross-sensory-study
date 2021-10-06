% This function arranges all data in one continuous record containing all
% trials in the right order of trialCount
% can help when analyzing prior trials etc.
% input: data in its original structure (for example: a few records for
%        each stimulus type)
% output: 1 continuous record containing all data in the right trials order
% fields included: trialCount, stimType, coherence, dir, response

% Shir Shalom 05/2018

function [dataRec]=oneDataRec(data)

priorsFlag=1; % Mark this flag as 1 if you have also priorResp in your structure
isPriorTrial=0; % Leave zero - it will indicate if a current trial is a prior trial

maxTrials=zeros(1,length(data));
for i=1:length(data)
    maxTrials(i)=max(data(i).Resp.trialCount);
end

lastTrial=max(maxTrials);

for i=1:lastTrial
    % find the current trial: record and index
    trialInd=[]; %trialInd is empty now
    trialRec=0;
    while isempty(trialInd) && trialRec<length(data)
        trialRec=trialRec+1;
        trialInd=find(data(trialRec).Resp.trialCount==i);
    end
    if isempty(trialInd) && priorsFlag==1 % then current trial is a prior trial
        trialRec=0;
        while isempty(trialInd) 
            trialRec=trialRec+1;
            trialInd=find(data(trialRec).PriorResp.trialCount==i);
        end
        isPriorTrial=1;
    end
    
    if ~isPriorTrial
        dataRec.trialCount(i)=data(trialRec).Resp.trialCount(trialInd);
        stimValueIdx = strmatch('STIMULUS_TYPE',{char(data(1).Rep(1,1).Trial(1).Param.name)},'exact');
        if ~isempty(stimValueIdx)
            dataRec.stimType(i)=data(trialRec).Rep.Trial(trialInd).Param(stimValueIdx).value;
        else
            dataRec.stimType(i)=NaN;
        end
        cohValueIdx = strmatch('STAR_MOTION_COHERENCE',{char(data(1).Rep(1,1).Trial(1).Param.name)},'exact');
        if ~isempty(cohValueIdx)
            dataRec.coherence(i)=data(trialRec).Rep.Trial(trialInd).Param(cohValueIdx).value;
        else
            dataRec.coherence(i)=NaN;
        end
        dataRec.dir(i)=data(trialRec).Resp.dir(trialInd);
        dataRec.response(i)=data(trialRec).Resp.response(trialInd);
        dataRec.isTestTrial(i)=1;
    else
        dataRec.trialCount(i)=data(trialRec).PriorResp.trialCount(trialInd);
        stimValueIdx = strmatch('STIMULUS_TYPE',{char(data(1).PriorRep(1,1).Trial(1).Param.name)},'exact');
        if ~isempty(stimValueIdx)
            dataRec.stimType(i)=data(trialRec).PriorRep.Trial(trialInd).Param(stimValueIdx).value;
        else
            dataRec.stimType(i)=NaN;
        end
        cohValueIdx = strmatch('STAR_MOTION_COHERENCE',{char(data(1).PriorRep(1,1).Trial(1).Param.name)},'exact');
        if ~isempty(cohValueIdx)
            dataRec.coherence(i)=data(trialRec).PriorRep.Trial(trialInd).Param(cohValueIdx).value;
        else
            dataRec.coherence(i)=NaN;
        end
        dataRec.dir(i)=data(trialRec).PriorResp.dir(trialInd);
        dataRec.response(i)=data(trialRec).PriorResp.response(trialInd);
        dataRec.isTestTrial(i)=0;
    end
    isPriorTrial=0;
end
end