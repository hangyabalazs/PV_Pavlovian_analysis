function quickanalysis_pavlovian_PV(cellids, resdir)
%QUICKANALYSIS_PAVLOVIAN_PV  Analysis of in vivo neuronal recordings.
%   QUICKANALYSIS_PAVLOVIAN_PV(CELLIDS, RESDIR) is designed as an offline
%   analysis tool for tetrode data and behavior, which can be executed in
%   an unsupervised manner on a daily bases. It gives a quick overview of
%   the experiment including response profiles of clustered neurons,
%   light-triggered PSTH and psychometric plots of behavior.
%   It relies on CellBase data handling system.

%   See also ADDNEWCELLS, PREALIGNSPIKES and VIEWCELL2B.
%   Balazs Hangya, Panna Hegedus
%   hangya.balazs@koki.hu
%   12/03/2023

% Input argument check
narginchk(0,2)

% Choose cellids
if nargin < 1 | isempty(cellids)
    cellids = select_pv_cells(cellids);
end

% Stop if error
dbstop if error

% Make results directory if not exist
if nargin < 2
    resdir = getpref('cellbase','datapath');
end

if ~isdir(resdir)
    mkdir(resdir)
end

% Prealign spikes for trial events
problem_behav_cellid = [];
for iC = 1:length(cellids)
    cellid = cellids(iC);
    try
        prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_pavlovian,'filetype','event','ifsave',1,'ifappend',0)
    catch
        disp('Error in prealignSpikes.');
        problem_behav_cellid = [problem_behav_cellid cellid];
    end
    
    
    % Cue response
    G = figure;
    pause(0.01)
    viewcell2b(cellid,'TriggerName','StimulusOn','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'DeliverAllFeedback'}},'Partitions','#TrialType','window',[-3 3])
    maximize_figure(G)
    
    cellidt = char(cellid);
    cellidt(cellidt=='.') = '_';
    fnm = fullfile(resdir,[cellidt '_StimulusOn.jpg']);   % save
    saveas(G,fnm)
    close(G)
    
    % Reward response
    H = figure;
    pause(0.01)
    viewcell2b(cellid,'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#AllReward','window',[-3 3])
    maximize_figure(H)
    
    cellidt = char(cellid);
    cellidt(cellidt=='.') = '_';
    fnm = fullfile(resdir,[cellidt '_Reward.jpg']);   % save
    saveas(H,fnm)
    close(H)
    
    % Punishment response
    I = figure;
    pause(0.01)
    viewcell2b(cellid,'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#Punishment','window',[-3 3])
    maximize_figure(I)
    
    cellidt = char(cellid);
    cellidt(cellidt=='.') = '_';
    fnm = fullfile(resdir,[cellidt '_Punishment.jpg']);   % save
    saveas(I,fnm)
    close(I)
    
     % Median split
    J = figure;
    pause(0.01)
    viewcell2b(cellid,'TriggerName','DeliverAllFeedback','SortEvent','TrialStart','sigma', 0.07,'eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#MedialSplit','window',[-3 3])
    maximize_figure(J)
    
    cellidt = char(cellid);
    cellidt(cellidt=='.') = '_';
    fnm = fullfile(resdir,[cellidt '_MedianSplit.jpg']);   % save
    saveas(J,fnm)
    close(J)
end




% View light-triggered raster and PSTH
TrigEvent = 'BurstOn';
SEvent = 'BurstOff';
win = [-0.2 0.5];
% parts = 'all';
parts = '#BurstNPulse';
dt = 0.001;
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
for iCell = 1:length(cellids)
    cellid = cellids(iCell);
    H = figure;
    viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',H,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','off')
    maximize_figure(H)
    
    cellidt = cellid{1};
    cellidt(cellidt=='.') = '_';
    fnm = fullfile(resdir,[cellidt '_LS.jpg']);   % save
    saveas(H,fnm)
    close(H)
end