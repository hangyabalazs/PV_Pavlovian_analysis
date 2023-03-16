function lick_psth_summary(pvcells,resdir,issave)
%LICK_PSTH_SUMMARY   Average lick PETH.
%   LICK_PSTH_SUMMARY plots average lick PETH (beam break time stamps
%   aligned to an event). 10 last sessions from mice are used.
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   20-Dec-2019

%   Code review: BH 12/20/19, 1/2/20, 4/8/20

% Input arguments
if nargin < 1 || isempty(pvcells)
    pvcells = select_pv_cells(pvcells);
end
if nargin < 2
    resdir = ('E:\auditory_pavlovian_cellbase\_paper_figs\code_revirew\Fig1_behavior\');   % result directory
end
if nargin < 3
    issave = true;
end

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Select sessions: mice w no feedback delay that had enough training 
cellids = pvcells'; % BF PV+ cell IDs
NumCells = length(cellids); % number of cells
Sessions =cell(1,NumCells); % sessions of recording
for i = 1:NumCells % extract session names from cellids
    [ratname,session,~,~] = cellid2tags(cellids{i});
    Sessions{i} = [ratname '_' session];
end
Sessions = unique(Sessions); % eliminate duplications
NumSessions = length(Sessions);

% Time window
wn = [-5 5];   % in seconds
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector

% PETH
[Hit_allpsth, FA_allpsth] = deal([]);
for iS = 1:NumSessions
    cSession = Sessions{iS}; % extract sessionIDs in the form of 1x2 cell
    [rat,remain] = strtok(cSession,'_');
    [session, ~] = strtok(remain, '_');
    sessionid = [{rat}, {session}];
    [spsth_hit, spsth_fa] = main(sessionid,wn,dt);
    Hit_allpsth = [Hit_allpsth; spsth_hit]; %#ok<*AGROW>
    FA_allpsth = [FA_allpsth; spsth_fa]; %#ok<*AGROW>
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
red = [216 41 0] / 255;
errorshade(time,nanmean(Hit_allpsth),nanstd(Hit_allpsth)/sqrt(size(Hit_allpsth,1)),...
    'LineColor',green,'ShadeColor',green)
hold on
errorshade(time,nanmean(FA_allpsth),nanstd(FA_allpsth)/sqrt(size(FA_allpsth,1)),...
    'LineColor',red,'ShadeColor',red)
if issave
    fnm = fullfile(resdir,'average_lickPSTH.fig');
    set(H, 'renderer', 'painters')
    fnm2 = fullfile(resdir,'average_lickPSTH.eps');
    fnm3 = fullfile(resdir,'average_lickPSTH.jpg');
    saveas(H,fnm)
    saveas(H,fnm2)
    saveas(H,fnm3)
    
end

% -------------------------------------------------------------------------
function [spsth_hit, spsth_fa] = main(cellid,wn,dt)

% Filter input
filterinput_hit = 'TrialType==1';
filterinput_fa = 'TrialType==2';

% Calcualte lick PSTH
[~, spsth_hit] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.08,'event_filter','custom','filterinput',filterinput_hit,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_fa] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.08,'event_filter','custom','filterinput',filterinput_fa,...
    'isadaptive',0,'maxtrialno',Inf);

% Lick raster
H = figure;
viewlick(cellid,'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav',...
    'ShowEvents',{{'StimulusOn' 'StimulusOff'}},...
    'Partitions','#TrialType','window',wn)
maximize_figure(H)