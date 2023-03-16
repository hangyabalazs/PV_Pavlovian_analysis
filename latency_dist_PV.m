function lat = latency_dist_PV(pvcells,sourcedir, stimuli)
%LATENCY_DIST_PV   Load response latency statistics.
%   LAT = LATENCY_DIST_PV(CELLIDS,DIR) loads latency statistics for cue,
%   reward and punishment (STIMULI) responsive cells of CELLIDS from SOURCEDIR. The results
%   are returned in a struct containing one field for each of cue,
%   reward and punishment. These individual fields contain an array of
%   five fields for start, peak, duration, maxvalue and minvalue for activated and inhibited
%   neurons, respectively.
%
%   See also ULTIMATE_PSTH and PSTH_STATS.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   19-Apr-2019

%   Code review: BH 12/20/19, 1/3/20

% Input arguments
if nargin < 1 || isempty(pvcells)
    pvcells = vpselectcells(pvcells);
end

% Responsive cells
cueresp = getvalue('cueresponse', pvcells);
rewardresp = getvalue('rewardresponse', pvcells);
punishmentresp = getvalue('punishresponse', pvcells);

% Load latency statistics
lat = struct;

switch stimuli
    case 'cue'
        p1_vpcells = pvcells(cueresp == 1); % excited by cue
        p2_vpcells = pvcells(cueresp == -1); % inhibited by cue
        stimstr = '_StimulusOn_TrialType.mat';
        numConditions = 1;
    case 'reward'
        p1_vpcells = pvcells(rewardresp == 1); % excited by reward
        p2_vpcells = pvcells(rewardresp == -1); % inhibited by reward
        stimstr = '_DeliverAllFeedback_RewardedTrials.mat';
        numConditions = 1;
    case 'punishment'
        p1_vpcells = pvcells(punishmentresp == 1); % excited by punishment
        p2_vpcells = pvcells(punishmentresp == -1); % inhibited by punishment
        stimstr = '_DeliverAllFeedback_PunishedTrials.mat';
        numConditions = 1;
    case 'punishment_expectation'
        p1_vpcells = pvcells(punishmentresp == 1); % excited by punishment
        p2_vpcells = pvcells(punishmentresp == -1); % inhibited by punishment
        stimstr = '_DeliverAllFeedback_Punishment.mat';
        numConditions = 2;
end

cells = [{p1_vpcells}, {p2_vpcells}];
[latency_start, latency_peak, latency_duration, minvalue, maxvalue, tags] = deal(cell(1,numConditions));
for l = 1:2
    numCells = length(cells{l});
    if numCells ~= 0
        for j = 1:numConditions
            [latency_start{j}, latency_peak{j}, latency_duration{j}] = deal(nan(1,numCells));
            for k = 1:numCells   % loop through cells
                cellid = char(cells{l}(k));
                [animalID, sessionID, tt, u]  = cellid2tags(cellid);
                sttc = [animalID '_' sessionID '_' num2str(tt) '_' num2str(u)];
                load([sourcedir '\' sttc stimstr]);   % load response analysis file
                close all;
                
                latency_start{j}(k) = stats1{1,j}.activation_start(1);
                latency_peak{j}(k) = stats1{1,j}.activation_peak(1);
                latency_duration{j}(k) = stats1{1,j}.activation_time(1);
                maxvalue{j}(k) = stats1{1,j}.maxvalue(1);
                minvalue{j}(k) = stats1{1,j}.minvalue(1);
                tags{j} = stats1{1,j}.tag;
                latency_i_start{j}(k) = stats1{1,j}.inhibition_start(1);
                latency_i_peak{j}(k) = stats1{1,j}.inhibition_peak(1);
                latency_i_duration{j}(k) = stats1{1,j}.inhibition_time(1);
            end
            
            % Output
            lat.(stimuli).excited.latency_start = latency_start;
            lat.(stimuli).inhibited.latency_start = latency_i_start;
            lat.(stimuli).excited.latency_peak = latency_peak;
            lat.(stimuli).inhibited.latency_peak = latency_i_peak;
            lat.(stimuli).excited.latency_duration = latency_duration;
            lat.(stimuli).inhibited.latency_duration = latency_i_duration;
            lat.(stimuli).maxvalue = maxvalue;
            lat.(stimuli).minvalue = minvalue;
            lat.(stimuli).tags = tags;
        end
    end
end