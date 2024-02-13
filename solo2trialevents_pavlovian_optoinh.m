function TE = solo2trialevents_pavlovian_optoinh(filepath,ifsave)
%SOLO2TRIALEVENTS_PAVLOVIAN_OPTOINH   Creates trial events structure.
%   SOLO2TRIALEVENTS2_PAVLOVIAN_OPTOINH(SESSPATH,IFSAVE) Creates TE trial events
%   structure. The events are extracted from the saved behavior file (FILEPATH).
%   Most of the events are given relative to the TrialStart event. IFSAVE
%   logical input argument controls whether to save the results.

%   Balazs Hangya, Panna Hegedus, Victoria Lyakhova
%   hangya.balazs@koki.hu
%   16/03/2023
           

if nargin < 2
    ifsave = 1;
end

% Determine filepath and load data
load(filepath);

% a dynamic calculation of the session name, required due to non-uniform
% naming format of the sessions
q = strfind(filepath,'_');
w = strfind(filepath,'.mat');

sessionID = filepath(q(end)+1:w);

% No. of trials (exclude trial 1 and last trial may not have been completed
ntrials = SessionData.nTrials;

% Preallocation of space
% Training session parameters
[TE.NTrials TE.TrialStart TE.TrialEnd TE.TrialType,TE.TrialType_AS,TE.TrialType_ANS] = deal(nan(1,ntrials));

% ITI
TE.ITIBeginning = nan(1,ntrials); % Begining of each ITI
TE.ITIDuration = nan(1,ntrials);
TE.ITIEnd = nan(1,ntrials);
TE.DeliverFeedback = nan(1,ntrials);
TE.DeliverAllFeedback = nan(1,ntrials);
TE.DeliverFreeWater = nan(1,ntrials);
TE.Feedback = nan(1,ntrials);
TE.Omission = nan(1,ntrials);
TE.Punishment = nan(1,ntrials);
TE.LEDOn = cell(1,ntrials);
TE.LEDOff = cell(1,ntrials);

% Stimulation
TE.StimulusOn = nan(1,ntrials);
TE.StimulusOff = nan(1,ntrials);
TE.StimulusDuration = nan(1,ntrials); % stimulus duration
TE.ToneType1 = nan(1,ntrials); % proportion of trialtype1 in the session settings
TE.ToneType2 = nan(1,ntrials); % proportion of trialtype2 in the session settings
TE.NoPunishment = nan(1,ntrials); % presence or absense of punishment in the session settings 1=no punishment, 0= punishment
TE.FreeWater = nan(1,ntrials); %how many portions of freewater the animal is given

% Response window
TE.Delay = nan(1,ntrials); % Delay period after the stimulus delivery
TE.ResponseWindow = nan(1,ntrials); % StimulusDuration + Delay
TE.ResponseWindowEnd = nan(1,ntrials);

% Reaction
TE.LickIn = cell(1,ntrials); % breaking the beam of the water port in headfixed protocol
TE.LickOut = cell(1,ntrials); % end of breaking the beam
TE.RT = nan(1,ntrials);
TE.GoRT = nan(1,ntrials);
TE.NoGoRT = nan(1,ntrials);
TE.Tup = cell(1,ntrials); % Timeups
TE.RewardValveTime = nan(1,ntrials);
TE.PunishValveTime = nan(1,ntrials);
TE.StartRewardValveTime = nan(1,ntrials);
TE.EndRewardValveTime = nan(1,ntrials);
TE.StartPunishValveTime = nan(1,ntrials);
TE.EndPunishValveTime = nan(1,ntrials);

% Outcomes
TE.Hit = nan(1,ntrials); % lick + reward
TE.AllReward = nan(1,ntrials); % lick + reward, partitioned to TrialTypes
TE.Reward = nan(1,ntrials); % lick + reward, partitioned to TrialTypes
TE.CorrectRejection = nan(1,ntrials); % no lick + punishment
TE.FalseAlarm = nan(1,ntrials); % lick + punishment
TE.Punishment = nan(1,ntrials);
TE.Omission = nan(1,ntrials);
TE.Miss = nan(1,ntrials); % no lick + missed reward
TE.LickOmission = nan(1,ntrials); % omission+ lick
TE.NoLickOmission = nan(1,ntrials); % omission+no lick
TE.FreeWaterTrial = nan(1,ntrials); % free water trials - deleted in the current version
TE.sessionID = cell(1,ntrials);
TE.NTrials = nan(1,ntrials);
TE.Animalname = cell(1,ntrials);

% Sync correction (state machine bug)
TE.NeedsSyncCorrection = nan(1,ntrials);

% Animal name
rawchar = length(filepath);
Animalname_pre = strfind(filepath,'PVINH');
Animalname_length = 7;
Animalname = filepath((Animalname_pre(1)):(Animalname_pre(1)+Animalname_length-1));

% Defining variables and parameters trialwise
for currentTrial = 1:ntrials
    TE.sessionID {1,currentTrial} = sessionID;
    TE.NTrials(1,currentTrial) = ntrials; %saving out the number of trials, just to be sure
    TE.Animalname(1,currentTrial) = cellstr(Animalname);
    TE.TrialStart = SessionData.TrialStartTimestamp; % Start of each trial
    
    %ITI
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'ITI')   % no ITI state in free water trials
        TE.ITIBeginning(1,currentTrial) = SessionData.RawEvents.Trial{1,currentTrial}.States.ITI(1);
        TE.ITIDuration(1,currentTrial) = SessionData.RawEvents.Trial{1,currentTrial}.States.ITI(2) - TE.ITIBeginning(1,currentTrial);
        TE.ITIEnd(1,currentTrial) = SessionData.RawEvents.Trial{1,currentTrial}.States.ITI(2);
    end
end

TE.TrialType = SessionData.TrialTypes; % 1 = Type1, 2 = Type2

% Outcomes
TE.Hit(SessionData.TrialOutcome==1) = 1; % lick + reward
TE.CorrectRejection(SessionData.TrialOutcome == 4) = 1; % no lick + punishment
TE.FalseAlarm(SessionData.TrialOutcome == 0) = 1; % lick + punishment
TE.Miss(SessionData.TrialOutcome == 2) = 1; % no lick + missed reward
TE.LickOmission(SessionData.TrialOutcome == 5) = 1; % omission + lick
TE.NoLickOmission(SessionData.TrialOutcome == 3) = 1 ; % omission + no lick

% For partitioning reinforcement
for currentTrial = 1:ntrials
    if SessionData.TrialOutcome(currentTrial) ==1 && SessionData.TrialTypes(currentTrial) ==1
        TE.Reward(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial)==1 && SessionData.TrialTypes(currentTrial) == 2
        TE.Reward(currentTrial) = 2;
    elseif SessionData.TrialOutcome(currentTrial) == 0  && SessionData.TrialTypes(currentTrial) ==1
        TE.Punishment(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial) == 4  && SessionData.TrialTypes(currentTrial) ==1
        TE.Punishment(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial) == 0  && SessionData.TrialTypes(currentTrial) ==2
        TE.Punishment(currentTrial) = 2;
    elseif SessionData.TrialOutcome(currentTrial) == 4  && SessionData.TrialTypes(currentTrial) ==2
        TE.Punishment(currentTrial) = 2;
    end
    
    if SessionData.TrialOutcome(currentTrial) == 3 || SessionData.TrialOutcome(currentTrial) == 5 && SessionData.TrialTypes(currentTrial) ==1
        TE.Omission(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial) == 3 || SessionData.TrialOutcome(currentTrial) == 5 && SessionData.TrialTypes(currentTrial) ==2
        TE.Omission(currentTrial) = 2;
    end
    
    if SessionData.TrialOutcome(currentTrial) ==1 && SessionData.TrialTypes(currentTrial) ==1
        TE.AllReward(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial)==1 && SessionData.TrialTypes(currentTrial) == 2
        TE.AllReward(currentTrial) = 2;
    elseif SessionData.TrialOutcome(currentTrial) ==2 && SessionData.TrialTypes(currentTrial) ==1
        TE.AllReward(currentTrial) = 1;
    elseif SessionData.TrialOutcome(currentTrial)==2 && SessionData.TrialTypes(currentTrial) == 2
        TE.AllReward(currentTrial) = 2;
    end
    
end

Sinx = find(~isnan(TE.Punishment)) + 1;

if ~isempty(Sinx)
    Sinx(end) = [];
    NSinx = find(isnan(TE.Punishment)) + 1;
    NSinx(end) = [];
else
    NSinx = [];
end

TE.TrialType_AS(Sinx) = TE.TrialType(Sinx);
TE.TrialType_ANS(NSinx) = TE.TrialType(NSinx);

for currentTrial = 1:ntrials
    if SessionData.TrialOutcome(currentTrial) ==1 || SessionData.TrialOutcome(currentTrial)==2
        TE.Feedback(currentTrial) = 1; % reward
    elseif SessionData.TrialOutcome(currentTrial)==0 || SessionData.TrialOutcome(currentTrial)==4
        TE.Feedback(currentTrial) = 2; % punishment
    elseif SessionData.TrialOutcome(currentTrial)==5 || SessionData.TrialOutcome(currentTrial)==3
        TE.Feedback(currentTrial) = 3; % omission
    end
end

% Stimulus
for currentTrial = 1:ntrials
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   % free water trials don't have these
        TE.StimulusOn(1,currentTrial) = SessionData.RawEvents.Trial{1,currentTrial}.States.StartStimulus(1:1);
        TE.Delay(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Delay(1:1)-SessionData.RawEvents.Trial{1,currentTrial}.States.Delay(2);
        
        TE.LEDOn{1,currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.States.NoLick(:,1);
        TE.LEDOff{1,currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.States.NoLick(:,2);
    end
    
    % DeliverFeedback
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'StartStimulus')   % free water trials don't have these
        if SessionData.TrialOutcome(currentTrial) == 4 || SessionData.TrialOutcome(currentTrial) == 0 % lick& no lick punishment
            TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
        elseif SessionData.TrialOutcome(currentTrial) == 1 % lick&reward
            TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
        elseif SessionData.TrialOutcome(currentTrial) == 2 && SessionData.TrialSettings(1).DirectDelivery == 1   %no lick& reward
            TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
        elseif SessionData.TrialOutcome(currentTrial) == 2 && SessionData.TrialSettings(1).DirectDelivery == 0   %no lick& no reward
            TE.DeliverFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
        end
        
        % DeliverAllFeedback
        if SessionData.TrialOutcome(currentTrial) == 4 || SessionData.TrialOutcome(currentTrial) == 0  %punishment
            TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1) ;
        elseif SessionData.TrialOutcome(currentTrial) == 1 %lick reward
            TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
        elseif SessionData.TrialOutcome(currentTrial) == 2 && SessionData.TrialSettings(1).DirectDelivery == 1
            TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
        elseif SessionData.TrialOutcome(currentTrial) == 2 && SessionData.TrialSettings(1).DirectDelivery ==0
            TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
        elseif SessionData.TrialOutcome(currentTrial) == 5 || SessionData.TrialOutcome(currentTrial) == 3 %omission
            TE.DeliverAllFeedback(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.PostUS(1:1);
        end
    else   % free water trials
        TE.DeliverFreeWater(currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
        TE.FreeWaterTrial(currentTrial) = 1;
    end
    
end

% Reaction
for currentTrial = 20:(ntrials-1)
    if SessionData.TrialOutcome(currentTrial) == 1
        TE.LickIn{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;% breaking the beam of the water port in headfixed protocol
    elseif SessionData.TrialOutcome(currentTrial) == 0
        TE.LickIn{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;% breaking the beam of the water port in headfixed protocol
    elseif SessionData.TrialOutcome(currentTrial) == 5
        TE.LickIn{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Port1In;% breaking the beam of the water port in headfixed protocol
    end
    
    rellicktimes = TE.LickIn{1,currentTrial}-TE.StimulusOn(1,currentTrial);
    firstlickinx = find(rellicktimes>0,1,'first');
    if ~isempty(firstlickinx)   % only if there were licks after stimulus
        TE.RT(1,currentTrial) = rellicktimes(firstlickinx);
    end
end
[TE.GoRT, TE.NoGoRT] = deal(TE.RT);
TE.GoRT(TE.TrialType~=1) = NaN;   % 'reaction time' in the likely rewarded trials
TE.NoGoRT(TE.TrialType~=2) = NaN;   % 'reaction time' in the likely punished trials

for currentTrial = 1:ntrials
    TE.Tup{1, currentTrial} = SessionData.RawEvents.Trial{1, currentTrial}.Events.Tup;
    TE.StartRewardValveTime(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(1:1);
    TE.EndRewardValveTime(1,currentTrial) = SessionData.RawEvents.Trial{1, currentTrial}.States.Reward(2);
    if isnan(TE.FreeWaterTrial(1,currentTrial)) % no Punish state in free water trials
        TE.StartPunishValveTime(1,currentTrial)= SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(1:1);
        TE.EndPunishValveTime(1,currentTrial)= SessionData.RawEvents.Trial{1, currentTrial}.States.Punish(2);
    end
end

% Sync correction
for currentTrial = 1:ntrials
    if isfield(SessionData.RawEvents.Trial{1,currentTrial}.States,'DeliverStimulus')   % TTL was sent with 25 ms delay by mistake
        TE.NeedsSyncCorrection(currentTrial) = 1;
    end
end
TE = shortenTE(TE,isnan(TE.FreeWaterTrial));  % skip free water trials for now - not trivial to sync (no TTLs)

% Stimulus settings
for currentTrial = 1:ntrials
    TE.ToneType1(1,currentTrial) = SessionData.TrialSettings(1, currentTrial).Type1; % proportion of trialtype1 in the session settings
    TE.ToneType2(1,currentTrial) = SessionData.TrialSettings(1, currentTrial).Type2; % proportion of trialtype2 in the session settings
    TE.NoPunishment(1,currentTrial) = SessionData.TrialSettings(1, currentTrial).NoPunish; % presence or absense of punishment in the session settings 1=no punishment, 0= punishment
    TE.FreeWater(1,currentTrial)=SessionData.TrialSettings(1, currentTrial).GUI.DirectWater; %how many portions of freewater the animal is given
end

% Save
if ifsave == 1
    
    % For each session, save the individual variables from the extracted TE
    savename = [Animalname '_' sessionID];
    savename2 = 'TE';
    savename3 = 'TrialEvents';
    save([fileparts(filepath) filesep savename],'-struct','TE')
    save([fileparts(filepath) filesep savename2],'-struct','TE')
    save([fileparts(filepath) filesep savename3],'-struct','TE')
end
disp('TE analysis complete')

% -------------------------------------------------------------------------
function TE2 = shortenTE(TE2,shinx)

% Eliminate behavioral trials
fnm = fieldnames(TE2);
for k = 1:length(fieldnames(TE2))
    TE2.(fnm{k}) = TE2.(fnm{k})(shinx);
end