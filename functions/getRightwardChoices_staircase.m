function [anaData]= getRightwardChoices_staircase(data,priorAnaData,specific_coherence)

% I previously used "getRightwardChoices which analyzed the variable "data.SavedInfo..."
% etc. For some reason the data structure changed and there is no longer a SavedInfo
% field? Rather the data itself is a structure of "# of stim" size. Hence
% the new function "getRightwardChoices_staircase"

%AZ updated to also take CombineDeltaPos and CombineDeltaNeg and VisualpriorL and VisualpriorR and VespriorL and VespriorR (priors still need to be done for combined) 
%AZ updated to take input specific_coherence. It will only return anaData for data with the specific coherence (vestibular is always returned) 

Ncrit=0; %only use headings with >Ncrit reps

iDirVes = 0;
dirArrayVes = [];
dirRepNumVes = [];
rightChoiceVes = [];
dirSeqVes=[];
NharderVes=0; NnotharderVes=0;
NeasierVes=0; NnoteasierVes=0;

iDirVisual = 0;
dirArrayVisual = [];
dirRepNumVisual = [];
rightChoiceVisual = [];
dirSeqVisual=[];
NharderVisual=0; NnotharderVisual=0;
NeasierVisual=0; NnoteasierVisual=0;

iDirVisualpriorL = 0;
dirArrayVisualpriorL = [];
dirRepNumVisualpriorL = [];
rightChoiceVisualpriorL = [];
dirSeqVisualpriorL=[];
NharderVisualpriorL=0; NnotharderVisualpriorL=0;
NeasierVisualpriorL=0; NnoteasierVisualpriorL=0;

iDirVisualpriorR = 0;
dirArrayVisualpriorR = [];
dirRepNumVisualpriorR = [];
rightChoiceVisualpriorR = [];
dirSeqVisualpriorR=[];
NharderVisualpriorR=0; NnotharderVisualpriorR=0;
NeasierVisualpriorR=0; NnoteasierVisualpriorR=0;

iDirVespriorL = 0;
dirArrayVespriorL = [];
dirRepNumVespriorL = [];
rightChoiceVespriorL = [];
dirSeqVespriorL=[];
NharderVespriorL=0; NnotharderVespriorL=0;
NeasierVespriorL=0; NnoteasierVespriorL=0;

iDirVespriorR = 0;
dirArrayVespriorR = [];
dirRepNumVespriorR = [];
rightChoiceVespriorR = [];
dirSeqVespriorR=[];
NharderVespriorR=0; NnotharderVespriorR=0;
NeasierVespriorR=0; NnoteasierVespriorR=0;

iDirCombine = 0;
dirArrayCombine = [];
dirRepNumCombine = [];
rightChoiceCombine = [];
dirSeqCombine=[];
NharderCombine=0; NnotharderCombine=0;
NeasierCombine=0; NnoteasierCombine=0;

iDirCombineDeltaPos = 0;
dirArrayCombineDeltaPos = [];
dirRepNumCombineDeltaPos = [];
rightChoiceCombineDeltaPos = [];
dirSeqCombineDeltaPos=[];
NharderCombineDeltaPos=0; NnotharderCombineDeltaPos=0;
NeasierCombineDeltaPos=0; NnoteasierCombineDeltaPos=0;

iDirCombineDeltaNeg = 0;
dirArrayCombineDeltaNeg = [];
dirRepNumCombineDeltaNeg = [];
rightChoiceCombineDeltaNeg = [];
dirSeqCombineDeltaNeg=[];
NharderCombineDeltaNeg=0; NnotharderCombineDeltaNeg=0;
NeasierCombineDeltaNeg=0; NnoteasierCombineDeltaNeg=0;

%this assumes that adapation values can be extracted from the first... 
idx = strmatch('ADAPTATION_ANGLE',{char(data(1).Rep(1,1).Trial(1).Param.name)},'exact');
if ~isempty(idx), anaData.ADAPTATION_ANGLE=data(1).Rep(1,1).Trial(1).Param(idx).value; else anaData.ADAPTATION_ANGLE=NaN; end
idx = strmatch('ADAPTATION_TYPE',{char(data(1).Rep(1,1).Trial(1).Param.name)},'exact');
if ~isempty(idx), anaData.ADAPTATION_TYPE=data(1).Rep(1,1).Trial(1).Param(idx).value; else anaData.ADAPTATION_TYPE=NaN; end

%added for specific_coherence. i.e. if we specify specific_coherence then only take the data with that coherence
if nargin>=3
    data_do=[];
    for s=1:length(data) %one for each stimulus type ves,vis & comb
        idx = strmatch('STAR_MOTION_COHERENCE',{char(data(s).Rep(1,1).Trial(1).Param.name)},'exact');
        if ~isempty(idx), coherence(s)=data(s).Rep(1,1).Trial(1).Param(idx).value; else coherence(s)=NaN; end
        idx = strmatch('STIMULUS_TYPE',{char(data(s).Rep(1,1).Trial(1).Param.name)},'exact');
        if ~isempty(idx), stim_type(s)=data(s).Rep(1,1).Trial(1).Param(idx).value; else stim_type(s)=NaN; end
        if (coherence(s)==specific_coherence) || (stim_type(s)==1) || (stim_type(s)==6) || (stim_type(s)==9); %include only that coherence. Always incl
            data_do = [data_do data(s)];
        end
    end
    anaData.STAR_MOTION_COHERENCE = specific_coherence;
else
    data_do = data;
    %just take coherence of the first block (NB note that this assumes constant coherence in the block)
    idx = strmatch('STAR_MOTION_COHERENCE',{char(data(1).Rep(1,1).Trial(1).Param.name)},'exact');
    if ~isempty(idx), anaData.STAR_MOTION_COHERENCE=data(1).Rep(1,1).Trial(1).Param(idx).value; else anaData.STAR_MOTION_COHERENCE=NaN; end
end


for s=1:length(data_do) %one for each stimulus type ves,vis & comb
    RepNum = length(data_do(s).Resp);   %repetition number
    for i=1:RepNum
        TrialNum = length(data_do(s).Resp(1,i).dir);   %trial number per repetition AZ: changed Rep(1,1) to Rep(1,RepNum)
        iStim = strmatch('STIMULUS_TYPE',{char(data_do(s).Rep(1,i).Trial(1).Param.name)},'exact');
        iJoy = strmatch('JOYBAR_DEGREE',{char(data_do(1).Rep(1,1).Trial(1).Param.name)},'exact');
        st=1;
        %st=2; %exclude first st-1 trials
        for j=st:TrialNum
            stim_type = data_do(s).Rep(1,i).Trial(1,j).Param(iStim).value;
            if ~isempty(iJoy), joy_val = data_do(s).Rep(1,i).Trial(1,j).Param(iJoy).value;
            else joy_val=NaN; end
            dir = data_do(s).Resp(1,i).dir(j);
            min_dir=min(abs(data_do(s).Resp(1,i).dir));
            max_dir=max(abs(data_do(s).Resp(1,i).dir));
            response = data_do(s).Resp(1,i).response(j);
            if response > 0 %AZ added  30-8-17 - Shir found bug that No choice was counted as left
            if stim_type == 1 %vestibular only
                dirSeqVes=[dirSeqVes dir];
                iInd = find(dirArrayVes == dir);
                if isempty(iInd)
                    iDirVes = iDirVes+1;
                    dirArrayVes(iDirVes) = dir;
                    dirRepNumVes(iDirVes) = 1;
                    rightChoiceVes(iDirVes) = 0;
                    joyVes{iDirVes}=joy_val;
                    iInd = iDirVes;
                else
                    dirRepNumVes(iInd)=dirRepNumVes(iInd)+1;
                    joyVes{iInd}=[joyVes{iInd} joy_val]; %add the value of the joystick selection to the array for that heading
                end
                if response == 2
                    right=1;
                else
                    right=0;
                end
                rightChoiceVes(iInd)=((dirRepNumVes(iInd)-1)*rightChoiceVes(iInd)+right)/dirRepNumVes(iInd);
                if j>1 && data_do(s).Resp(1,i).corr(j-1) ==1 && abs(data_do(s).Resp(1,i).dir(j-1))~=min_dir %Note: if heading is the minimum it can't get harder
                    if abs(data_do(s).Resp(1,i).dir(j))<abs(data_do(s).Resp(1,i).dir(j-1))
                        NharderVes=NharderVes+1;
                    else
                        NnotharderVes=NnotharderVes+1;
                    end
                elseif j>1 && data_do(s).Resp(1,i).corr(j-1) ==0 && abs(data_do(s).Resp(1,i).dir(j-1))~=max_dir %Note: if heading is the maximum it can't get easier
                    if abs(data_do(s).Resp(1,i).dir(j))>abs(data_do(s).Resp(1,i).dir(j-1))
                        NeasierVes=NeasierVes+1;
                    else
                        NnoteasierVes=NnoteasierVes+1;
                    end
                end
                
            elseif stim_type == 2  %visual
                dirSeqVisual=[dirSeqVisual dir];
                iInd = find(dirArrayVisual == dir);
                if isempty(iInd)
                    iDirVisual = iDirVisual+1;
                    dirArrayVisual(iDirVisual) = dir;
                    dirRepNumVisual(iDirVisual) = 1;
                    rightChoiceVisual(iDirVisual) = 0;
                    joyVisual{iDirVisual}=joy_val;
                    iInd = iDirVisual;
                else
                    dirRepNumVisual(iInd)=dirRepNumVisual(iInd)+1;
                    joyVisual{iInd}=[joyVisual{iInd} joy_val]; %add the value of the joystick selection to the array for that heading
                end
                if response == 2
                    right=1;
                else
                    right=0;
                end
                rightChoiceVisual(iInd)=((dirRepNumVisual(iInd)-1)*rightChoiceVisual(iInd)+right)/dirRepNumVisual(iInd);
                if j>1 && data_do(s).Resp(1,i).corr(j-1) ==1 && abs(data_do(s).Resp(1,i).dir(j-1))~=min_dir %Note: if heading is the minimum it can't get harder
                    if abs(data_do(s).Resp(1,i).dir(j))<abs(data_do(s).Resp(1,i).dir(j-1))
                        NharderVisual=NharderVisual+1;
                    else
                        NnotharderVisual=NnotharderVisual+1;
                    end
                elseif j>1 && data_do(s).Resp(1,i).corr(j-1) ==0 && abs(data_do(s).Resp(1,i).dir(j-1))~=max_dir %Note: if heading is the maximum it can't get easier
                    if abs(data_do(s).Resp(1,i).dir(j))>abs(data_do(s).Resp(1,i).dir(j-1))
                        NeasierVisual=NeasierVisual+1;
                    else
                        NnoteasierVisual=NnoteasierVisual+1;
                    end
                end
                
            elseif stim_type == 3 %combine
                dirSeqCombine=[dirSeqCombine dir];
                iInd = find(dirArrayCombine == dir);
                if isempty(iInd)
                    iDirCombine = iDirCombine+1;
                    dirArrayCombine(iDirCombine) = dir;
                    dirRepNumCombine(iDirCombine) = 1;
                    rightChoiceCombine(iDirCombine) = 0;
                    joyCombine{iDirCombine}=joy_val;
                    iInd = iDirCombine;
                else
                    dirRepNumCombine(iInd)=dirRepNumCombine(iInd)+1;
                    joyCombine{iInd}=[joyCombine{iInd} joy_val]; %add the value of the joystick selection to the array for that heading
                end
                if response == 2
                    right=1;
                else
                    right=0;
                end
                rightChoiceCombine(iInd)=((dirRepNumCombine(iInd)-1)*rightChoiceCombine(iInd)+right)/dirRepNumCombine(iInd);
                if j>1 && data_do(s).Resp(1,i).corr(j-1) ==1 && abs(data_do(s).Resp(1,i).dir(j-1))~=min_dir %Note: if heading is the minimum it can't get harder
                    if abs(data_do(s).Resp(1,i).dir(j))<abs(data_do(s).Resp(1,i).dir(j-1))
                        NharderCombine=NharderCombine+1;
                    else
                        NnotharderCombine=NnotharderCombine+1;
                    end
                elseif j>1 && data_do(s).Resp(1,i).corr(j-1) ==0 && abs(data_do(s).Resp(1,i).dir(j-1))~=max_dir %Note: if heading is the maximum it can't get easier
                    if abs(data_do(s).Resp(1,i).dir(j))>abs(data_do(s).Resp(1,i).dir(j-1))
                        NeasierCombine=NeasierCombine+1;
                    else
                        NnoteasierCombine=NnoteasierCombine+1;
                    end
                end
                
                
            elseif stim_type == 4 %combine delta pos
                dirSeqCombineDeltaPos=[dirSeqCombineDeltaPos dir];
                iInd = find(dirArrayCombineDeltaPos == dir);
                if isempty(iInd)
                    iDirCombineDeltaPos = iDirCombineDeltaPos+1;
                    dirArrayCombineDeltaPos(iDirCombineDeltaPos) = dir;
                    dirRepNumCombineDeltaPos(iDirCombineDeltaPos) = 1;
                    rightChoiceCombineDeltaPos(iDirCombineDeltaPos) = 0;
                    joyCombineDeltaPos{iDirCombineDeltaPos}=joy_val;
                    iInd = iDirCombineDeltaPos;
                else
                    dirRepNumCombineDeltaPos(iInd)=dirRepNumCombineDeltaPos(iInd)+1;
                    joyCombineDeltaPos{iInd}=[joyCombineDeltaPos{iInd} joy_val]; %add the value of the joystick selection to the array for that heading
                end
                if response == 2
                    right=1;
                else
                    right=0;
                end
                rightChoiceCombineDeltaPos(iInd)=((dirRepNumCombineDeltaPos(iInd)-1)*rightChoiceCombineDeltaPos(iInd)+right)/dirRepNumCombineDeltaPos(iInd);
                if j>1 && data_do(s).Resp(1,i).corr(j-1) ==1 && abs(data_do(s).Resp(1,i).dir(j-1))~=min_dir %Note: if heading is the minimum it can't get harder
                    if abs(data_do(s).Resp(1,i).dir(j))<abs(data_do(s).Resp(1,i).dir(j-1))
                        NharderCombineDeltaPos=NharderCombineDeltaPos+1;
                    else
                        NnotharderCombineDeltaPos=NnotharderCombineDeltaPos+1;
                    end
                elseif j>1 && data_do(s).Resp(1,i).corr(j-1) ==0 && abs(data_do(s).Resp(1,i).dir(j-1))~=max_dir %Note: if heading is the maximum it can't get easier
                    if abs(data_do(s).Resp(1,i).dir(j))>abs(data_do(s).Resp(1,i).dir(j-1))
                        NeasierCombineDeltaPos=NeasierCombineDeltaPos+1;
                    else
                        NnoteasierCombineDeltaPos=NnoteasierCombineDeltaPos+1;
                    end
                end
                
            elseif stim_type == 5 %combine delta neg
                dirSeqCombineDeltaNeg=[dirSeqCombineDeltaNeg dir];
                iInd = find(dirArrayCombineDeltaNeg == dir);
                if isempty(iInd)
                    iDirCombineDeltaNeg = iDirCombineDeltaNeg+1;
                    dirArrayCombineDeltaNeg(iDirCombineDeltaNeg) = dir;
                    dirRepNumCombineDeltaNeg(iDirCombineDeltaNeg) = 1;
                    rightChoiceCombineDeltaNeg(iDirCombineDeltaNeg) = 0;
                    joyCombineDeltaNeg{iDirCombineDeltaNeg}=joy_val;
                    iInd = iDirCombineDeltaNeg;
                else
                    dirRepNumCombineDeltaNeg(iInd)=dirRepNumCombineDeltaNeg(iInd)+1;
                    joyCombineDeltaNeg{iInd}=[joyCombineDeltaNeg{iInd} joy_val]; %add the value of the joystick selection to the array for that heading
                end
                if response == 2
                    right=1;
                else
                    right=0;
                end
                rightChoiceCombineDeltaNeg(iInd)=((dirRepNumCombineDeltaNeg(iInd)-1)*rightChoiceCombineDeltaNeg(iInd)+right)/dirRepNumCombineDeltaNeg(iInd);
                if j>1 && data_do(s).Resp(1,i).corr(j-1) ==1 && abs(data_do(s).Resp(1,i).dir(j-1))~=min_dir %Note: if heading is the minimum it can't get harder
                    if abs(data_do(s).Resp(1,i).dir(j))<abs(data_do(s).Resp(1,i).dir(j-1))
                        NharderCombineDeltaNeg=NharderCombineDeltaNeg+1;
                    else
                        NnotharderCombineDeltaNeg=NnotharderCombineDeltaNeg+1;
                    end
                elseif j>1 && data_do(s).Resp(1,i).corr(j-1) ==0 && abs(data_do(s).Resp(1,i).dir(j-1))~=max_dir %Note: if heading is the maximum it can't get easier
                    if abs(data_do(s).Resp(1,i).dir(j))>abs(data_do(s).Resp(1,i).dir(j-1))
                        NeasierCombineDeltaNeg=NeasierCombineDeltaNeg+1;
                    else
                        NnoteasierCombineDeltaNeg=NnoteasierCombineDeltaNeg+1;
                    end
                end
            
            elseif stim_type == 6  %Ves prior L
                dirSeqVespriorL=[dirSeqVespriorL dir];
                iInd = find(dirArrayVespriorL == dir);
                if isempty(iInd)
                    iDirVespriorL = iDirVespriorL+1;
                    dirArrayVespriorL(iDirVespriorL) = dir;
                    dirRepNumVespriorL(iDirVespriorL) = 1;
                    rightChoiceVespriorL(iDirVespriorL) = 0;
                    joyVespriorL{iDirVespriorL}=joy_val;
                    iInd = iDirVespriorL;
                else
                    dirRepNumVespriorL(iInd)=dirRepNumVespriorL(iInd)+1;
                    joyVespriorL{iInd}=[joyVespriorL{iInd} joy_val]; %add the value of the joystick selection to the array for that heading
                end
                if response == 2
                    right=1;
                else
                    right=0;
                end
                rightChoiceVespriorL(iInd)=((dirRepNumVespriorL(iInd)-1)*rightChoiceVespriorL(iInd)+right)/dirRepNumVespriorL(iInd);
                if j>1 && data_do(s).Resp(1,i).corr(j-1) ==1 && abs(data_do(s).Resp(1,i).dir(j-1))~=min_dir %Note: if heading is the minimum it can't get harder
                    if abs(data_do(s).Resp(1,i).dir(j))<abs(data_do(s).Resp(1,i).dir(j-1))
                        NharderVespriorL=NharderVespriorL+1;
                    else
                        NnotharderVespriorL=NnotharderVespriorL+1;
                    end
                elseif j>1 && data_do(s).Resp(1,i).corr(j-1) ==0 && abs(data_do(s).Resp(1,i).dir(j-1))~=max_dir %Note: if heading is the maximum it can't get easier
                    if abs(data_do(s).Resp(1,i).dir(j))>abs(data_do(s).Resp(1,i).dir(j-1))
                        NeasierVespriorL=NeasierVespriorL+1;
                    else
                        NnoteasierVespriorL=NnoteasierVespriorL+1;
                    end
                end
                
            elseif stim_type == 7  %visual prior L
                dirSeqVisualpriorL=[dirSeqVisualpriorL dir];
                iInd = find(dirArrayVisualpriorL == dir);
                if isempty(iInd)
                    iDirVisualpriorL = iDirVisualpriorL+1;
                    dirArrayVisualpriorL(iDirVisualpriorL) = dir;
                    dirRepNumVisualpriorL(iDirVisualpriorL) = 1;
                    rightChoiceVisualpriorL(iDirVisualpriorL) = 0;
                    joyVisualpriorL{iDirVisualpriorL}=joy_val;
                    iInd = iDirVisualpriorL;
                else
                    dirRepNumVisualpriorL(iInd)=dirRepNumVisualpriorL(iInd)+1;
                    joyVisualpriorL{iInd}=[joyVisualpriorL{iInd} joy_val]; %add the value of the joystick selection to the array for that heading
                end
                if response == 2
                    right=1;
                else
                    right=0;
                end
                rightChoiceVisualpriorL(iInd)=((dirRepNumVisualpriorL(iInd)-1)*rightChoiceVisualpriorL(iInd)+right)/dirRepNumVisualpriorL(iInd);
                if j>1 && data_do(s).Resp(1,i).corr(j-1) ==1 && abs(data_do(s).Resp(1,i).dir(j-1))~=min_dir %Note: if heading is the minimum it can't get harder
                    if abs(data_do(s).Resp(1,i).dir(j))<abs(data_do(s).Resp(1,i).dir(j-1))
                        NharderVisualpriorL=NharderVisualpriorL+1;
                    else
                        NnotharderVisualpriorL=NnotharderVisualpriorL+1;
                    end
                elseif j>1 && data_do(s).Resp(1,i).corr(j-1) ==0 && abs(data_do(s).Resp(1,i).dir(j-1))~=max_dir %Note: if heading is the maximum it can't get easier
                    if abs(data_do(s).Resp(1,i).dir(j))>abs(data_do(s).Resp(1,i).dir(j-1))
                        NeasierVisualpriorL=NeasierVisualpriorL+1;
                    else
                        NnoteasierVisualpriorL=NnoteasierVisualpriorL+1;
                    end
                end
                
            elseif stim_type == 9  %ves prior R
                dirSeqVespriorR=[dirSeqVespriorR dir];
                iInd = find(dirArrayVespriorR == dir);
                if isempty(iInd)
                    iDirVespriorR = iDirVespriorR+1;
                    dirArrayVespriorR(iDirVespriorR) = dir;
                    dirRepNumVespriorR(iDirVespriorR) = 1;
                    rightChoiceVespriorR(iDirVespriorR) = 0;
                    joyVespriorR{iDirVespriorR}=joy_val;
                    iInd = iDirVespriorR;
                else
                    dirRepNumVespriorR(iInd)=dirRepNumVespriorR(iInd)+1;
                    joyVespriorR{iInd}=[joyVespriorR{iInd} joy_val]; %add the value of the joystick selection to the array for that heading
                end
                if response == 2
                    right=1;
                else
                    right=0;
                end
                rightChoiceVespriorR(iInd)=((dirRepNumVespriorR(iInd)-1)*rightChoiceVespriorR(iInd)+right)/dirRepNumVespriorR(iInd);
                if j>1 && data_do(s).Resp(1,i).corr(j-1) ==1 && abs(data_do(s).Resp(1,i).dir(j-1))~=min_dir %Note: if heading is the minimum it can't get harder
                    if abs(data_do(s).Resp(1,i).dir(j))<abs(data_do(s).Resp(1,i).dir(j-1))
                        NharderVespriorR=NharderVespriorR+1;
                    else
                        NnotharderVespriorR=NnotharderVespriorR+1;
                    end
                elseif j>1 && data_do(s).Resp(1,i).corr(j-1) ==0 && abs(data_do(s).Resp(1,i).dir(j-1))~=max_dir %Note: if heading is the maximum it can't get easier
                    if abs(data_do(s).Resp(1,i).dir(j))>abs(data_do(s).Resp(1,i).dir(j-1))
                        NeasierVespriorR=NeasierVespriorR+1;
                    else
                        NnoteasierVespriorR=NnoteasierVespriorR+1;
                    end
                end
                
            elseif stim_type == 10  %visual prior R
                dirSeqVisualpriorR=[dirSeqVisualpriorR dir];
                iInd = find(dirArrayVisualpriorR == dir);
                if isempty(iInd)
                    iDirVisualpriorR = iDirVisualpriorR+1;
                    dirArrayVisualpriorR(iDirVisualpriorR) = dir;
                    dirRepNumVisualpriorR(iDirVisualpriorR) = 1;
                    rightChoiceVisualpriorR(iDirVisualpriorR) = 0;
                    joyVisualpriorR{iDirVisualpriorR}=joy_val;
                    iInd = iDirVisualpriorR;
                else
                    dirRepNumVisualpriorR(iInd)=dirRepNumVisualpriorR(iInd)+1;
                    joyVisualpriorR{iInd}=[joyVisualpriorR{iInd} joy_val]; %add the value of the joystick selection to the array for that heading
                end
                if response == 2
                    right=1;
                else
                    right=0;
                end
                rightChoiceVisualpriorR(iInd)=((dirRepNumVisualpriorR(iInd)-1)*rightChoiceVisualpriorR(iInd)+right)/dirRepNumVisualpriorR(iInd);
                if j>1 && data_do(s).Resp(1,i).corr(j-1) ==1 && abs(data_do(s).Resp(1,i).dir(j-1))~=min_dir %Note: if heading is the minimum it can't get harder
                    if abs(data_do(s).Resp(1,i).dir(j))<abs(data_do(s).Resp(1,i).dir(j-1))
                        NharderVisualpriorR=NharderVisualpriorR+1;
                    else
                        NnotharderVisualpriorR=NnotharderVisualpriorR+1;
                    end
                elseif j>1 && data_do(s).Resp(1,i).corr(j-1) ==0 && abs(data_do(s).Resp(1,i).dir(j-1))~=max_dir %Note: if heading is the maximum it can't get easier
                    if abs(data_do(s).Resp(1,i).dir(j))>abs(data_do(s).Resp(1,i).dir(j-1))
                        NeasierVisualpriorR=NeasierVisualpriorR+1;
                    else
                        NnoteasierVisualpriorR=NnoteasierVisualpriorR+1;
                    end
                end
                
                
                
            end
            end
        end
    end
end
if ~isempty(dirRepNumVes) %vest
    [sortDirVes, sortIndVes] = sort(dirArrayVes, 2);
    sortRightVes = rightChoiceVes(sortIndVes);
    sortdirRepNumVes=dirRepNumVes(sortIndVes);
    anaData.dirVes = sortDirVes;
    anaData.rightVes = sortRightVes;
    anaData.dirSeqVes = dirSeqVes;
    tempVes=cell2mat_nan(joyVes');
    anaData.joyVes=tempVes(sortIndVes,:);
    if ~isempty(sortRightVes)
        pfit_input = cat(2, sortDirVes', sortRightVes', sortdirRepNumVes');
        selVes=sortdirRepNumVes>Ncrit; %only take headings with Ncrit+1 reps
        IVes=pfit_input(selVes,2);
        if exist('priorAnaData') && ~isempty(priorAnaData), pfitDataVes = getPfitData(pfit_input(selVes,:),sortDirVes(selVes),priorAnaData.thresh95CIVes); %if a thresh 'prior' is given
        elseif all([IVes(1:min([find(IVes==1); length(IVes)])-1)==0; IVes(min([find(IVes==1); length(IVes)]):end)==1]) %if there are no intermediate values cannot calculate thresh.
            pfitDataVes.bias = cum_gaussfit_max1(pfit_input(selVes,:));
            anaData.threshVes=NaN;
        else
            pfit_input(pfit_input(1:end,2)==0,2)=1e-5; %AZ pfit program doesn't like pefect 0s (performs badly)
            pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-5); %AZ pfit program doesn't like perfect 1s
            pfitDataVes = getPfitData(pfit_input(selVes,:),sortDirVes(selVes));
            anaData.threshVes = pfitDataVes.thresh;
            anaData.xiVes = pfitDataVes.xi;
            anaData.pfitcurveVes = pfitDataVes.pfitcurve;
            anaData.bias95CIVes = pfitDataVes.bias95CI;
            anaData.biasSDVes = pfitDataVes.biasSD;
            anaData.thresh95CIVes = pfitDataVes.thresh95CI;
            anaData.RsquareVes = pfitDataVes.Rsquare;
        end
        anaData.biasVes = pfitDataVes.bias;
        anaData.NVes=[sum(sortdirRepNumVes(sortDirVes<0)) sum(sortdirRepNumVes(sortDirVes==0)) sum(sortdirRepNumVes(sortDirVes>0))];
        anaData.RepNumVes=sortdirRepNumVes;
        anaData.PerCorrectVes = 100*(sum(sortRightVes(sortDirVes>0) .* sortdirRepNumVes(sortDirVes>0)) + sum((1-sortRightVes(sortDirVes<0)) .* sortdirRepNumVes(sortDirVes<0)))/sum(sortdirRepNumVes);
        anaData.PerHarderVes = NharderVes/(NharderVes+NnotharderVes);
        anaData.PerEasierVes = NeasierVes/(NeasierVes+NnoteasierVes);
    end
else
    anaData.threshVes=NaN;
    anaData.biasVes=NaN;
end
if ~isempty(dirRepNumVisual) %vis
    [sortDirVisual, sortIndVisual] = sort(dirArrayVisual, 2);
    sortRightVisual = rightChoiceVisual(sortIndVisual);
    sortdirRepNumVisual=dirRepNumVisual(sortIndVisual);
    anaData.dirVisual = sortDirVisual;
    anaData.rightVisual = sortRightVisual;
    anaData.dirSeqVisual = dirSeqVisual;
    tempVisual=cell2mat_nan(joyVisual');
    anaData.joyVisual=tempVisual(sortIndVisual,:);
    if ~isempty(sortRightVisual)
        pfit_input = cat(2, sortDirVisual', sortRightVisual', sortdirRepNumVisual');
        selVisual=sortdirRepNumVisual>Ncrit; %only take headings with Ncrit+1 reps
        IVisual=pfit_input(selVisual,2);
        if exist('priorAnaData') && ~isempty(priorAnaData), pfitDataVisual = getPfitData(pfit_input(selVisual,:),sortDirVisual(selVisual),priorAnaData.thresh95CIVisual); %if a thresh 'prior' is given
        elseif all([IVisual(1:min([find(IVisual==1); length(IVisual)])-1)==0; IVisual(min([find(IVisual==1); length(IVisual)]):end)==1]) %if there are no intermediate values cannot calculate thresh.
            pfitDataVisual.bias = cum_gaussfit_max1(pfit_input(selVisual,:));
            anaData.threshVisual=NaN;
        else
            pfit_input(pfit_input(1:end,2)==0,2)=1e-5; %AZ pfit program doesn't like pefect 0s (performs badly)
            pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-5); %AZ pfit program doesn't like perfect 1s
            pfitDataVisual = getPfitData(pfit_input(selVisual,:),sortDirVisual(selVisual));
            anaData.threshVisual = pfitDataVisual.thresh;
            anaData.xiVisual = pfitDataVisual.xi;
            anaData.pfitcurveVisual = pfitDataVisual.pfitcurve;
            anaData.bias95CIVisual = pfitDataVisual.bias95CI;
            anaData.biasSDVisual = pfitDataVisual.biasSD;
            anaData.thresh95CIVisual = pfitDataVisual.thresh95CI;
            anaData.RsquareVisual = pfitDataVisual.Rsquare;
        end
        anaData.biasVisual = pfitDataVisual.bias;
        anaData.NVisual=[sum(sortdirRepNumVisual(sortDirVisual<0)) sum(sortdirRepNumVisual(sortDirVisual==0)) sum(sortdirRepNumVisual(sortDirVisual>0))];
        anaData.RepNumVisual=sortdirRepNumVisual;
        anaData.PerCorrectVisual = 100*(sum(sortRightVisual(sortDirVisual>0) .* sortdirRepNumVisual(sortDirVisual>0)) + sum((1-sortRightVisual(sortDirVisual<0)) .* sortdirRepNumVisual(sortDirVisual<0)))/sum(sortdirRepNumVisual);
        anaData.PerHarderVisual = NharderVisual/(NharderVisual+NnotharderVisual);
        anaData.PerEasierVisual = NeasierVisual/(NeasierVisual+NnoteasierVisual);
    end
else
    anaData.threshVisual=NaN;
    anaData.biasVisual=NaN;
end



if ~isempty(dirRepNumCombine) %combined
    [sortDirCombine, sortIndCombine] = sort(dirArrayCombine, 2);
    sortRightCombine = rightChoiceCombine(sortIndCombine);
    sortdirRepNumCombine=dirRepNumCombine(sortIndCombine);
    anaData.dirCombine = sortDirCombine;
    anaData.rightCombine = sortRightCombine;
    anaData.dirSeqCombine = dirSeqCombine;
    tempCombine=cell2mat_nan(joyCombine');
    anaData.joyCombine=tempCombine(sortIndCombine,:);
    if ~isempty(sortRightCombine)
        pfit_input = cat(2, sortDirCombine', sortRightCombine', sortdirRepNumCombine');
        selCombine=sortdirRepNumCombine>Ncrit;
        %why don't I have check "if all((pfit_input(selCombine,2)==0) | (pfit_input(selCombine,2)==1))" like for visual/vest? - add?
        pfit_input(pfit_input(1:end,2)==0,2)=1e-5; %AZ pfit program doesn't like pefect 0s (performs badly)
        pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-5); %AZ pfit program doesn't like perfect 1s
        pfitDataCombine = getPfitData(pfit_input(selCombine,:),sortDirCombine(selCombine));
        anaData.biasCombine = pfitDataCombine.bias;
        anaData.threshCombine = pfitDataCombine.thresh;
        anaData.xiCombine = pfitDataCombine.xi;
        anaData.pfitcurveCombine = pfitDataCombine.pfitcurve;
        anaData.bias95CICombine = pfitDataCombine.bias95CI;
        anaData.biasSDCombine = pfitDataCombine.biasSD;
        anaData.thresh95CICombine = pfitDataCombine.thresh95CI;
        anaData.RsquareCombine = pfitDataCombine.Rsquare;
        anaData.PerCorrectCombine = 100*(sum(sortRightCombine(sortDirCombine>0) .* sortdirRepNumCombine(sortDirCombine>0)) + sum((1-sortRightCombine(sortDirCombine<0)) .* sortdirRepNumCombine(sortDirCombine<0)))/sum(sortdirRepNumCombine);
    end
    anaData.biasCombine = pfitDataCombine.bias;
    anaData.NCombine=[sum(sortdirRepNumCombine(sortDirCombine<0)) sum(sortdirRepNumCombine(sortDirCombine==0)) sum(sortdirRepNumCombine(sortDirCombine>0))];
    anaData.RepNumCombine=sortdirRepNumCombine;
    anaData.PerCorrectCombine = 100*(sum(sortRightCombine(sortDirCombine>0) .* sortdirRepNumCombine(sortDirCombine>0)) + sum((1-sortRightCombine(sortDirCombine<0)) .* sortdirRepNumCombine(sortDirCombine<0)))/sum(sortdirRepNumCombine);
    anaData.PerHarderCombine = NharderCombine/(NharderCombine+NnotharderCombine);
    anaData.PerEasierCombine = NeasierCombine/(NeasierCombine+NnoteasierCombine);
else
    anaData.threshCombine=NaN;
    anaData.biasCombine=NaN;
end

if ~isempty(dirRepNumCombineDeltaPos) %CombineDeltaPos
    [sortDirCombineDeltaPos, sortIndCombineDeltaPos] = sort(dirArrayCombineDeltaPos, 2);
    sortRightCombineDeltaPos = rightChoiceCombineDeltaPos(sortIndCombineDeltaPos);
    sortdirRepNumCombineDeltaPos=dirRepNumCombineDeltaPos(sortIndCombineDeltaPos);
    anaData.dirCombineDeltaPos = sortDirCombineDeltaPos;
    anaData.rightCombineDeltaPos = sortRightCombineDeltaPos;
    anaData.dirSeqCombineDeltaPos = dirSeqCombineDeltaPos;
    tempCombineDeltaPos=cell2mat_nan(joyCombineDeltaPos');
    anaData.joyCombineDeltaPos=tempCombineDeltaPos(sortIndCombineDeltaPos,:);
    if ~isempty(sortRightCombineDeltaPos)
        pfit_input = cat(2, sortDirCombineDeltaPos', sortRightCombineDeltaPos', sortdirRepNumCombineDeltaPos');
        selCombineDeltaPos=sortdirRepNumCombineDeltaPos>Ncrit;
        %why don't I have check "if all((pfit_input(selCombineDeltaPos,2)==0) | (pfit_input(selCombineDeltaPos,2)==1))" like for visual/vest? - add?
        pfit_input(pfit_input(1:end,2)==0,2)=1e-5; %AZ pfit program doesn't like pefect 0s (performs badly)
        pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-5); %AZ pfit program doesn't like perfect 1s
        pfitDataCombineDeltaPos = getPfitData(pfit_input(selCombineDeltaPos,:),sortDirCombineDeltaPos(selCombineDeltaPos));
        anaData.biasCombineDeltaPos = pfitDataCombineDeltaPos.bias;
        anaData.threshCombineDeltaPos = pfitDataCombineDeltaPos.thresh;
        anaData.xiCombineDeltaPos = pfitDataCombineDeltaPos.xi;
        anaData.pfitcurveCombineDeltaPos = pfitDataCombineDeltaPos.pfitcurve;
        anaData.bias95CICombineDeltaPos = pfitDataCombineDeltaPos.bias95CI;
        anaData.biasSDCombineDeltaPos = pfitDataCombineDeltaPos.biasSD;
        anaData.thresh95CICombineDeltaPos = pfitDataCombineDeltaPos.thresh95CI;
        anaData.RsquareCombineDeltaPos = pfitDataCombineDeltaPos.Rsquare;
        anaData.PerCorrectCombineDeltaPos = 100*(sum(sortRightCombineDeltaPos(sortDirCombineDeltaPos>0) .* sortdirRepNumCombineDeltaPos(sortDirCombineDeltaPos>0)) + sum((1-sortRightCombineDeltaPos(sortDirCombineDeltaPos<0)) .* sortdirRepNumCombineDeltaPos(sortDirCombineDeltaPos<0)))/sum(sortdirRepNumCombineDeltaPos);
    end
    anaData.biasCombineDeltaPos = pfitDataCombineDeltaPos.bias;
    anaData.NCombineDeltaPos=[sum(sortdirRepNumCombineDeltaPos(sortDirCombineDeltaPos<0)) sum(sortdirRepNumCombineDeltaPos(sortDirCombineDeltaPos==0)) sum(sortdirRepNumCombineDeltaPos(sortDirCombineDeltaPos>0))];
    anaData.RepNumCombineDeltaPos=sortdirRepNumCombineDeltaPos;
    anaData.PerCorrectCombineDeltaPos = 100*(sum(sortRightCombineDeltaPos(sortDirCombineDeltaPos>0) .* sortdirRepNumCombineDeltaPos(sortDirCombineDeltaPos>0)) + sum((1-sortRightCombineDeltaPos(sortDirCombineDeltaPos<0)) .* sortdirRepNumCombineDeltaPos(sortDirCombineDeltaPos<0)))/sum(sortdirRepNumCombineDeltaPos);
    anaData.PerHarderCombineDeltaPos = NharderCombineDeltaPos/(NharderCombineDeltaPos+NnotharderCombineDeltaPos);
    anaData.PerEasierCombineDeltaPos = NeasierCombineDeltaPos/(NeasierCombineDeltaPos+NnoteasierCombineDeltaPos);
else
    anaData.threshCombineDeltaPos=NaN;
    anaData.biasCombineDeltaPos=NaN;
end

if ~isempty(dirRepNumCombineDeltaNeg) %CombineDeltaNeg
    [sortDirCombineDeltaNeg, sortIndCombineDeltaNeg] = sort(dirArrayCombineDeltaNeg, 2);
    sortRightCombineDeltaNeg = rightChoiceCombineDeltaNeg(sortIndCombineDeltaNeg);
    sortdirRepNumCombineDeltaNeg=dirRepNumCombineDeltaNeg(sortIndCombineDeltaNeg);
    anaData.dirCombineDeltaNeg = sortDirCombineDeltaNeg;
    anaData.rightCombineDeltaNeg = sortRightCombineDeltaNeg;
    anaData.dirSeqCombineDeltaNeg = dirSeqCombineDeltaNeg;
    tempCombineDeltaNeg=cell2mat_nan(joyCombineDeltaNeg');
    anaData.joyCombineDeltaNeg=tempCombineDeltaNeg(sortIndCombineDeltaNeg,:);
    if ~isempty(sortRightCombineDeltaNeg)
        pfit_input = cat(2, sortDirCombineDeltaNeg', sortRightCombineDeltaNeg', sortdirRepNumCombineDeltaNeg');
        selCombineDeltaNeg=sortdirRepNumCombineDeltaNeg>Ncrit;
        %why don't I have check "if all((pfit_input(selCombineDeltaNeg,2)==0) | (pfit_input(selCombineDeltaNeg,2)==1))" like for visual/vest? - add?
        pfit_input(pfit_input(1:end,2)==0,2)=1e-5; %AZ pfit program doesn't like pefect 0s (performs badly)
        pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-5); %AZ pfit program doesn't like perfect 1s
        pfitDataCombineDeltaNeg = getPfitData(pfit_input(selCombineDeltaNeg,:),sortDirCombineDeltaNeg(selCombineDeltaNeg));
        anaData.biasCombineDeltaNeg = pfitDataCombineDeltaNeg.bias;
        anaData.threshCombineDeltaNeg = pfitDataCombineDeltaNeg.thresh;
        anaData.xiCombineDeltaNeg = pfitDataCombineDeltaNeg.xi;
        anaData.pfitcurveCombineDeltaNeg = pfitDataCombineDeltaNeg.pfitcurve;
        anaData.bias95CICombineDeltaNeg = pfitDataCombineDeltaNeg.bias95CI;
        anaData.biasSDCombineDeltaNeg = pfitDataCombineDeltaNeg.biasSD;
        anaData.thresh95CICombineDeltaNeg = pfitDataCombineDeltaNeg.thresh95CI;
        anaData.RsquareCombineDeltaNeg = pfitDataCombineDeltaNeg.Rsquare;
        anaData.PerCorrectCombineDeltaNeg = 100*(sum(sortRightCombineDeltaNeg(sortDirCombineDeltaNeg>0) .* sortdirRepNumCombineDeltaNeg(sortDirCombineDeltaNeg>0)) + sum((1-sortRightCombineDeltaNeg(sortDirCombineDeltaNeg<0)) .* sortdirRepNumCombineDeltaNeg(sortDirCombineDeltaNeg<0)))/sum(sortdirRepNumCombineDeltaNeg);
    end
    anaData.biasCombineDeltaNeg = pfitDataCombineDeltaNeg.bias;
    anaData.NCombineDeltaNeg=[sum(sortdirRepNumCombineDeltaNeg(sortDirCombineDeltaNeg<0)) sum(sortdirRepNumCombineDeltaNeg(sortDirCombineDeltaNeg==0)) sum(sortdirRepNumCombineDeltaNeg(sortDirCombineDeltaNeg>0))];
    anaData.RepNumCombineDeltaNeg=sortdirRepNumCombineDeltaNeg;
    anaData.PerCorrectCombineDeltaNeg = 100*(sum(sortRightCombineDeltaNeg(sortDirCombineDeltaNeg>0) .* sortdirRepNumCombineDeltaNeg(sortDirCombineDeltaNeg>0)) + sum((1-sortRightCombineDeltaNeg(sortDirCombineDeltaNeg<0)) .* sortdirRepNumCombineDeltaNeg(sortDirCombineDeltaNeg<0)))/sum(sortdirRepNumCombineDeltaNeg);
    anaData.PerHarderCombineDeltaNeg = NharderCombineDeltaNeg/(NharderCombineDeltaNeg+NnotharderCombineDeltaNeg);
    anaData.PerEasierCombineDeltaNeg = NeasierCombineDeltaNeg/(NeasierCombineDeltaNeg+NnoteasierCombineDeltaNeg);
else
    anaData.threshCombineDeltaNeg=NaN;
    anaData.biasCombineDeltaNeg=NaN;
end




if ~isempty(dirRepNumVisualpriorL) %VisualpriorL prior L
    [sortDirVisualpriorL, sortIndVisualpriorL] = sort(dirArrayVisualpriorL, 2);
    sortRightVisualpriorL = rightChoiceVisualpriorL(sortIndVisualpriorL);
    sortdirRepNumVisualpriorL=dirRepNumVisualpriorL(sortIndVisualpriorL);
    anaData.dirVisualpriorL = sortDirVisualpriorL;
    anaData.rightVisualpriorL = sortRightVisualpriorL;
    anaData.dirSeqVisualpriorL = dirSeqVisualpriorL;
    tempVisualpriorL=cell2mat_nan(joyVisualpriorL');
    anaData.joyVisualpriorL=tempVisualpriorL(sortIndVisualpriorL,:);
    if ~isempty(sortRightVisualpriorL)
        pfit_input = cat(2, sortDirVisualpriorL', sortRightVisualpriorL', sortdirRepNumVisualpriorL');
        selVisualpriorL=sortdirRepNumVisualpriorL>Ncrit; %only take headings with Ncrit+1 reps
        IVisualpriorL=pfit_input(selVisualpriorL,2);
        if exist('priorAnaData') && ~isempty(priorAnaData), pfitDataVisualpriorL = getPfitData(pfit_input(selVisualpriorL,:),sortDirVisualpriorL(selVisualpriorL),priorAnaData.thresh95CIVisualpriorL); %if a thresh 'prior' is given
        elseif all([IVisualpriorL(1:min([find(IVisualpriorL==1); length(IVisualpriorL)])-1)==0; IVisualpriorL(min([find(IVisualpriorL==1); length(IVisualpriorL)]):end)==1]) %if there are no intermediate values cannot calculate thresh.
            pfitDataVisualpriorL.bias = cum_gaussfit_max1(pfit_input(selVisualpriorL,:));
            anaData.threshVisualpriorL=NaN;
        else
            pfit_input(pfit_input(1:end,2)==0,2)=1e-5; %AZ pfit program doesn't like pefect 0s (performs badly)
            pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-5); %AZ pfit program doesn't like perfect 1s
            pfitDataVisualpriorL = getPfitData(pfit_input(selVisualpriorL,:),sortDirVisualpriorL(selVisualpriorL));
            anaData.threshVisualpriorL = pfitDataVisualpriorL.thresh;
            anaData.xiVisualpriorL = pfitDataVisualpriorL.xi;
            anaData.pfitcurveVisualpriorL = pfitDataVisualpriorL.pfitcurve;
            anaData.bias95CIVisualpriorL = pfitDataVisualpriorL.bias95CI;
            anaData.biasSDVisualpriorL = pfitDataVisualpriorL.biasSD;
            anaData.thresh95CIVisualpriorL = pfitDataVisualpriorL.thresh95CI;
            anaData.RsquareVisualpriorL = pfitDataVisualpriorL.Rsquare;
        end
        anaData.biasVisualpriorL = pfitDataVisualpriorL.bias;
        anaData.NVisualpriorL=[sum(sortdirRepNumVisualpriorL(sortDirVisualpriorL<0)) sum(sortdirRepNumVisualpriorL(sortDirVisualpriorL==0)) sum(sortdirRepNumVisualpriorL(sortDirVisualpriorL>0))];
        anaData.RepNumVisualpriorL=sortdirRepNumVisualpriorL;
        anaData.PerCorrectVisualpriorL = 100*(sum(sortRightVisualpriorL(sortDirVisualpriorL>0) .* sortdirRepNumVisualpriorL(sortDirVisualpriorL>0)) + sum((1-sortRightVisualpriorL(sortDirVisualpriorL<0)) .* sortdirRepNumVisualpriorL(sortDirVisualpriorL<0)))/sum(sortdirRepNumVisualpriorL);
        anaData.PerHarderVisualpriorL = NharderVisualpriorL/(NharderVisualpriorL+NnotharderVisualpriorL);
        anaData.PerEasierVisualpriorL = NeasierVisualpriorL/(NeasierVisualpriorL+NnoteasierVisualpriorL);
    end
else
    anaData.threshVisualpriorL=NaN;
    anaData.biasVisualpriorL=NaN;
end
if ~isempty(dirRepNumVisualpriorR) %VisualpriorR prior R
    [sortDirVisualpriorR, sortIndVisualpriorR] = sort(dirArrayVisualpriorR, 2);
    sortRightVisualpriorR = rightChoiceVisualpriorR(sortIndVisualpriorR);
    sortdirRepNumVisualpriorR=dirRepNumVisualpriorR(sortIndVisualpriorR);
    anaData.dirVisualpriorR = sortDirVisualpriorR;
    anaData.rightVisualpriorR = sortRightVisualpriorR;
    anaData.dirSeqVisualpriorR = dirSeqVisualpriorR;
    tempVisualpriorR=cell2mat_nan(joyVisualpriorR');
    anaData.joyVisualpriorR=tempVisualpriorR(sortIndVisualpriorR,:);
    if ~isempty(sortRightVisualpriorR)
        pfit_input = cat(2, sortDirVisualpriorR', sortRightVisualpriorR', sortdirRepNumVisualpriorR');
        selVisualpriorR=sortdirRepNumVisualpriorR>Ncrit; %only take headings with Ncrit+1 reps
        IVisualpriorR=pfit_input(selVisualpriorR,2);
        if exist('priorAnaData') && ~isempty(priorAnaData), pfitDataVisualpriorR = getPfitData(pfit_input(selVisualpriorR,:),sortDirVisualpriorR(selVisualpriorR),priorAnaData.thresh95CIVisualpriorR); %if a thresh 'prior' is given
        elseif all([IVisualpriorR(1:min([find(IVisualpriorR==1); length(IVisualpriorR)])-1)==0; IVisualpriorR(min([find(IVisualpriorR==1); length(IVisualpriorR)]):end)==1]) %if there are no intermediate values cannot calculate thresh.
            pfitDataVisualpriorR.bias = cum_gaussfit_max1(pfit_input(selVisualpriorR,:));
            anaData.threshVisualpriorR=NaN;
        else
            pfit_input(pfit_input(1:end,2)==0,2)=1e-5; %AZ pfit program doesn't like pefect 0s (performs badly)
            pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-5); %AZ pfit program doesn't like perfect 1s
            pfitDataVisualpriorR = getPfitData(pfit_input(selVisualpriorR,:),sortDirVisualpriorR(selVisualpriorR));
            anaData.threshVisualpriorR = pfitDataVisualpriorR.thresh;
            anaData.xiVisualpriorR = pfitDataVisualpriorR.xi;
            anaData.pfitcurveVisualpriorR = pfitDataVisualpriorR.pfitcurve;
            anaData.bias95CIVisualpriorR = pfitDataVisualpriorR.bias95CI;
            anaData.biasSDVisualpriorR = pfitDataVisualpriorR.biasSD;
            anaData.thresh95CIVisualpriorR = pfitDataVisualpriorR.thresh95CI;
            anaData.RsquareVisualpriorR = pfitDataVisualpriorR.Rsquare;
        end
        anaData.biasVisualpriorR = pfitDataVisualpriorR.bias;
        anaData.NVisualpriorR=[sum(sortdirRepNumVisualpriorR(sortDirVisualpriorR<0)) sum(sortdirRepNumVisualpriorR(sortDirVisualpriorR==0)) sum(sortdirRepNumVisualpriorR(sortDirVisualpriorR>0))];
        anaData.RepNumVisualpriorR=sortdirRepNumVisualpriorR;
        anaData.PerCorrectVisualpriorR = 100*(sum(sortRightVisualpriorR(sortDirVisualpriorR>0) .* sortdirRepNumVisualpriorR(sortDirVisualpriorR>0)) + sum((1-sortRightVisualpriorR(sortDirVisualpriorR<0)) .* sortdirRepNumVisualpriorR(sortDirVisualpriorR<0)))/sum(sortdirRepNumVisualpriorR);
        anaData.PerHarderVisualpriorR = NharderVisualpriorR/(NharderVisualpriorR+NnotharderVisualpriorR);
        anaData.PerEasierVisualpriorR = NeasierVisualpriorR/(NeasierVisualpriorR+NnoteasierVisualpriorR);
    end
else
    anaData.threshVisualpriorR=NaN;
    anaData.biasVisualpriorR=NaN;
end

if ~isempty(dirRepNumVespriorL) %VespriorL prior L
    [sortDirVespriorL, sortIndVespriorL] = sort(dirArrayVespriorL, 2);
    sortRightVespriorL = rightChoiceVespriorL(sortIndVespriorL);
    sortdirRepNumVespriorL=dirRepNumVespriorL(sortIndVespriorL);
    anaData.dirVespriorL = sortDirVespriorL;
    anaData.rightVespriorL = sortRightVespriorL;
    anaData.dirSeqVespriorL = dirSeqVespriorL;
    tempVespriorL=cell2mat_nan(joyVespriorL');
    anaData.joyVespriorL=tempVespriorL(sortIndVespriorL,:);
    if ~isempty(sortRightVespriorL)
        pfit_input = cat(2, sortDirVespriorL', sortRightVespriorL', sortdirRepNumVespriorL');
        selVespriorL=sortdirRepNumVespriorL>Ncrit; %only take headings with Ncrit+1 reps
        IVespriorL=pfit_input(selVespriorL,2);
        if exist('priorAnaData') && ~isempty(priorAnaData), pfitDataVespriorL = getPfitData(pfit_input(selVespriorL,:),sortDirVespriorL(selVespriorL),priorAnaData.thresh95CIVespriorL); %if a thresh 'prior' is given
        elseif all([IVespriorL(1:min([find(IVespriorL==1); length(IVespriorL)])-1)==0; IVespriorL(min([find(IVespriorL==1); length(IVespriorL)]):end)==1]) %if there are no intermediate values cannot calculate thresh.
            pfitDataVespriorL.bias = cum_gaussfit_max1(pfit_input(selVespriorL,:));
            anaData.threshVespriorL=NaN;
        else
            pfit_input(pfit_input(1:end,2)==0,2)=1e-5; %AZ pfit program doesn't like pefect 0s (performs badly)
            pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-5); %AZ pfit program doesn't like perfect 1s
            pfitDataVespriorL = getPfitData(pfit_input(selVespriorL,:),sortDirVespriorL(selVespriorL));
            anaData.threshVespriorL = pfitDataVespriorL.thresh;
            anaData.xiVespriorL = pfitDataVespriorL.xi;
            anaData.pfitcurveVespriorL = pfitDataVespriorL.pfitcurve;
            anaData.bias95CIVespriorL = pfitDataVespriorL.bias95CI;
            anaData.biasSDVespriorL = pfitDataVespriorL.biasSD;
            anaData.thresh95CIVespriorL = pfitDataVespriorL.thresh95CI;
            anaData.RsquareVespriorL = pfitDataVespriorL.Rsquare;
        end
        anaData.biasVespriorL = pfitDataVespriorL.bias;
        anaData.NVespriorL=[sum(sortdirRepNumVespriorL(sortDirVespriorL<0)) sum(sortdirRepNumVespriorL(sortDirVespriorL==0)) sum(sortdirRepNumVespriorL(sortDirVespriorL>0))];
        anaData.RepNumVespriorL=sortdirRepNumVespriorL;
        anaData.PerCorrectVespriorL = 100*(sum(sortRightVespriorL(sortDirVespriorL>0) .* sortdirRepNumVespriorL(sortDirVespriorL>0)) + sum((1-sortRightVespriorL(sortDirVespriorL<0)) .* sortdirRepNumVespriorL(sortDirVespriorL<0)))/sum(sortdirRepNumVespriorL);
        anaData.PerHarderVespriorL = NharderVespriorL/(NharderVespriorL+NnotharderVespriorL);
        anaData.PerEasierVespriorL = NeasierVespriorL/(NeasierVespriorL+NnoteasierVespriorL);
    end
else
    anaData.threshVespriorL=NaN;
    anaData.biasVespriorL=NaN;
end
if ~isempty(dirRepNumVespriorR) %VespriorR prior R
    [sortDirVespriorR, sortIndVespriorR] = sort(dirArrayVespriorR, 2);
    sortRightVespriorR = rightChoiceVespriorR(sortIndVespriorR);
    sortdirRepNumVespriorR=dirRepNumVespriorR(sortIndVespriorR);
    anaData.dirVespriorR = sortDirVespriorR;
    anaData.rightVespriorR = sortRightVespriorR;
    anaData.dirSeqVespriorR = dirSeqVespriorR;
    tempVespriorR=cell2mat_nan(joyVespriorR');
    anaData.joyVespriorR=tempVespriorR(sortIndVespriorR,:);
    if ~isempty(sortRightVespriorR)
        pfit_input = cat(2, sortDirVespriorR', sortRightVespriorR', sortdirRepNumVespriorR');
        selVespriorR=sortdirRepNumVespriorR>Ncrit; %only take headings with Ncrit+1 reps
        IVespriorR=pfit_input(selVespriorR,2);
        if exist('priorAnaData') && ~isempty(priorAnaData), pfitDataVespriorR = getPfitData(pfit_input(selVespriorR,:),sortDirVespriorR(selVespriorR),priorAnaData.thresh95CIVespriorR); %if a thresh 'prior' is given
        elseif all([IVespriorR(1:min([find(IVespriorR==1); length(IVespriorR)])-1)==0; IVespriorR(min([find(IVespriorR==1); length(IVespriorR)]):end)==1]) %if there are no intermediate values cannot calculate thresh.
            pfitDataVespriorR.bias = cum_gaussfit_max1(pfit_input(selVespriorR,:));
            anaData.threshVespriorR=NaN;
        else
            pfit_input(pfit_input(1:end,2)==0,2)=1e-5; %AZ pfit program doesn't like pefect 0s (performs badly)
            pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-5); %AZ pfit program doesn't like perfect 1s
            pfitDataVespriorR = getPfitData(pfit_input(selVespriorR,:),sortDirVespriorR(selVespriorR));
            anaData.threshVespriorR = pfitDataVespriorR.thresh;
            anaData.xiVespriorR = pfitDataVespriorR.xi;
            anaData.pfitcurveVespriorR = pfitDataVespriorR.pfitcurve;
            anaData.bias95CIVespriorR = pfitDataVespriorR.bias95CI;
            anaData.biasSDVespriorR = pfitDataVespriorR.biasSD;
            anaData.thresh95CIVespriorR = pfitDataVespriorR.thresh95CI;
            anaData.RsquareVespriorR = pfitDataVespriorR.Rsquare;
        end
        anaData.biasVespriorR = pfitDataVespriorR.bias;
        anaData.NVespriorR=[sum(sortdirRepNumVespriorR(sortDirVespriorR<0)) sum(sortdirRepNumVespriorR(sortDirVespriorR==0)) sum(sortdirRepNumVespriorR(sortDirVespriorR>0))];
        anaData.RepNumVespriorR=sortdirRepNumVespriorR;
        anaData.PerCorrectVespriorR = 100*(sum(sortRightVespriorR(sortDirVespriorR>0) .* sortdirRepNumVespriorR(sortDirVespriorR>0)) + sum((1-sortRightVespriorR(sortDirVespriorR<0)) .* sortdirRepNumVespriorR(sortDirVespriorR<0)))/sum(sortdirRepNumVespriorR);
        anaData.PerHarderVespriorR = NharderVespriorR/(NharderVespriorR+NnotharderVespriorR);
        anaData.PerEasierVespriorR = NeasierVespriorR/(NeasierVespriorR+NnoteasierVespriorR);
    end
else
    anaData.threshVespriorR=NaN;
    anaData.biasVespriorR=NaN;
end