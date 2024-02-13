function lick_psth_summary_FF(sessions,resdir,issave)
%LICK_PSTH_SUMMARY_FF Average lick PETH for fiber photometry recordings.
%   LICK_PSTH_SUMMARY_FF(SESSIONS,RESDIR,ISSAVE) plots average lick PETH
%   (beam break time stamps aligned to an event). Session indicated in
%   SESSIONS as pairs of animal and folder names in a cell are used. PETHS 
%   are saved to RESDIR if indicated in the ISSAVE logalical variable. 
%
%   See also ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   20-Dec-2019

%   Code review: BH 12/20/19, 1/2/20, 4/8/20

% Input arguments
if nargin < 1 || isempty(sessions)
    sessions = [{{'SOM-HDB3' '230807a'}} {{'SOM-HDB3' '230808a'}} {{'SOM-HDB3' '230809a'}} ...
     {{'SOM-HDB3' '230810a'}} {{'SOM-HDB3' '230811a'}} {{'SOM-HDB4' '230802a'}} ...
    {{'SOM-HDB4' '230803a'}} {{'SOM-HDB4' '230804a'}} {{'SOM-HDB4' '230807a'}}...
    {{'SOM-HDB4' '230808a'}} {{'SOM-HDB5' '230802a'}} {{'SOM-HDB5' '230803a'}}...
    {{'SOM-HDB5' '230804a'}} {{'SOM-HDB5' '230807a'}} {{'SOM-HDB5' '230808a'}}...
    {{'SOM-HDB6' '230802a'}} {{'SOM-HDB6' '230803a'}} {{'SOM-HDB6' '230804a'}}...
    {{'SOM-HDB6' '230807a'}} {{'SOM-HDB6' '230808a'}}];
end
if nargin < 2
    resdir = 'F:\SOM\Plots\';   % result directory
end
if nargin < 3
    issave = true;
end

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Select sessions: mice w no feedback delay that had enough training 
cellids = sessions'; % BF PV+ cell IDs
NumCells = length(cellids); % number of cells
Session_list =cell(1,NumCells); % sessions of recording
for i = 1:NumCells % extract session names from cellids
    c_session = sessions{i};
    ratname = c_session{1};
    session = c_session{2};
    Session_list{i} = [ratname '_' session];
end
Session_list = unique(Session_list); % eliminate duplications
NumSessions = length(Session_list);

% Time window
wn = [-5 5];   % in seconds
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector

% PETH
[Hit_allpsth, FA_allpsth] = deal([]);
for iS = 1:NumSessions
    cSession = Session_list{iS}; % extract sessionIDs in the form of 1x2 cell
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