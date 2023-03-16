function lick_psth_summary_allcond(cellids,resdir,issave)
%LICK_PSTH_SUMMARY_ALLCOND(CELLIDS, RESDIR, ISSAVE)   Average lick PETH.
%   ACH_LICK_PSTH_SUMMARY_ALLCOND plots average lick PETH (beam break time stamps
%   aligned to an event) for all conditions of a probabilistic pavlovian
%   condiitoning task (likely reward ,unlikely punishment, likely punishment
%   unlikely reward). All sessions from mice are used.
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   07-June-2021

%   Code review:

% Input arguments
if nargin < 1 || isempty(cellids)
    cellids = select_pv_cells(cellids);
end

if nargin < 2
    resdir = ('D:\auditory_pavlovian_PV_cellbase\_data_analysis\lick_psth_allcond');   % result directory
end
if nargin < 3
    issave = true;
end

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Select sessions (all sessions of all animals)
cellids = cellids'; % BF PV+ cell IDs
NumCells = length(cellids); % number of cells
Sessions =cell(1,NumCells); % sessions of recording
for i = 1:NumCells % extract session names from cellids
    [ratname,session,~,~] = cellid2tags(cellids{i});
    Sessions{i} = [ratname '_' session];
end
Sessions = unique(Sessions); % eliminate duplications
NumSessions = length(Sessions);

% Time window
wn = [-3 3];   % in seconds
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector

% PETH
[Hit_allpsth, Hit_allpsth2, FA_allpsth, FA_allpsth2] = deal([]);
for iS = 1:NumSessions
    cSession = Sessions{iS}; % extract sessionIDs in the form of 1x2 cell
    [rat,remain] = strtok(cSession,'_');
    [session, ~] = strtok(remain, '_');
    sessionid = [{rat}, {session}];
    [spsth_hit, spsth_hit2, spsth_fa, spsth_fa2] = main(sessionid,wn,dt);
    Hit_allpsth = [Hit_allpsth; spsth_hit];
    Hit_allpsth2 = [Hit_allpsth2; spsth_hit2];
    FA_allpsth = [FA_allpsth; spsth_fa]; 
    FA_allpsth2 = [FA_allpsth2; spsth_fa2];
    H = gcf;
    if issave
        cellidt = [sessionid{1} '_' sessionid{2}];
        fnm = fullfile(resdir, [cellidt '_LICK.fig']);   % save
        saveas(H,fnm)
        fnm = fullfile(resdir, [cellidt '_LICK.jpg']);
        saveas(H,fnm)
    end
    close(H)
end

% Plot & save
H = figure;
green = [51 204 51] / 255;   % colors for plotting
blue = [0 0 1];
red = [216 41 0] / 255;
dark_red = [0.5 0 0];
errorshade(time,nanmean(Hit_allpsth),nanstd(Hit_allpsth)/sqrt(size(Hit_allpsth,1)),...
    'LineColor',green,'ShadeColor',green)
hold on
errorshade(time,nanmean(Hit_allpsth2),nanstd(Hit_allpsth2)/sqrt(size(Hit_allpsth2,1)),...
    'LineColor',blue,'ShadeColor',blue)
hold on
errorshade(time,nanmean(FA_allpsth),nanstd(FA_allpsth)/sqrt(size(FA_allpsth,1)),...
    'LineColor',red,'ShadeColor',red)
hold on
errorshade(time,nanmean(FA_allpsth2),nanstd(FA_allpsth2)/sqrt(size(FA_allpsth2,1)),...
    'LineColor',dark_red,'ShadeColor',dark_red)
if issave
    fnm = fullfile(resdir,'average_lickPSTH_allcond.fig');
    saveas(H,fnm)
        set(H, 'renderer', 'painters')
    fnm2 = fullfile(resdir,'average_lickPSTH_allcond.eps');
    fnm3 = fullfile(resdir,'average_lickPSTH_allcond.jpg');
    saveas(H,fnm)
    saveas(H,fnm2)
    saveas(H,fnm3)
end

% -------------------------------------------------------------------------
function [spsth_hit, spsth_hit2, spsth_fa, spsth_fa2] = main(cellid,wn,dt)

% Filter input
filterinput_hit = 'TrialType==1&AllReward==1';
filterinput_hit2 = 'TrialType==1&Punishment==1';
filterinput_fa = 'TrialType==2&AllReward==2';
filterinput_fa2 = 'TrialType==2&Punishment==2';


% Calcualte lick PSTH
[~, spsth_hit] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_hit,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_fa] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_fa,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_hit2] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_hit2,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_fa2] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.05,'event_filter','custom','filterinput',filterinput_fa2,...
    'isadaptive',0,'maxtrialno',Inf);

% Lick raster
H = figure;
viewlick(cellid,'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','sigma', 0.05,...
    'ShowEvents',{{'StimulusOn' 'StimulusOff'}},...
    'Partitions','#TrialType','window',wn)
maximize_figure(H)