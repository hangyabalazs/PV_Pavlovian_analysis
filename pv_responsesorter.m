function pv_responsesorter(cellids,issave,response_resdir,responsespec)
%PV_RESPONSESORTER   Peri-event time histogram.
%   PV_RESPONSESORTER(CELLIDS,ISSAVE) calculates non-adaptive PSTHs for a
%   set of cells (see ULTIMATE_PSTH) aligned to stimulus onset and feedback
%   delivery. Statistical tests are performed to probe significant firing
%   rate changes of responses (see PSTH_STATS). Indicators for significant
%   responses (p < 0.001, one-sided Mann-Whitney U-test) are added to
%   CellBase as properties.
%   Input parameters:
%       CELLIDS - list of cell IDs; if empty or not specified, all
%           well-separated cells are selected (ID>20, L-ratio<0.15; see
%           LRATIO) from ventral pallidum
%       ISSAVE - controls saving
%       RESDIR - results directory
%       RESPONSESPEC - 'cue' for cue response, 'rew' for reward response
%           and 'pun' for punishment response
%
%   See also ULTIMATE_PSTH and LRATIO.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   11-Nov-2018

%   Code review: BH 7/24/19, 12/2/19, 4/8/20

% Directories
if ~isfolder(response_resdir)
    mkdir(response_resdir)
end

% Input argument check
narginchk(0,4);
if nargin < 2
    issave = true;   % default saving behavior
end
if nargin < 1
    pv_cells = select_pv_cells([]);   % all well-isolated units
else
    pv_cells = cellids;
end
numCells = length(pv_cells);

% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

% Preprocess trial events file if the needed variables do no exist
% for p = 1:numCells
%     cell2check = pv_cells{p};   % current cell
%     preproc_TE(cell2check)
% end

% Raster + PSTH
switch responsespec
    case 'cue'
        
        % Raster + PSTH aligned to stimulus onset
        alignevent = 'StimulusOn';   % trigger event
        shevent = 'DeliverAllFeedback';  % show-event
        partition = '#TrialType';   % partition trials
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-2 0];   % baseline window for MW-test
        twin = [0 1];   % test-window for MW-test
        
        % Add propoerty for grouping
        propname = 'cueresponse';
        if ~ismember(propname,listtag('prop'))
            insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname)
        end
        
        % PSTH
        for iC = 1:numCells
            cellid = pv_cells{iC};   % current cell
            stats1 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir);
            
            % Add property to CellBase - inhibition taken over activation!
            if length(stats1) > 1  % if there are two sounds
                if stats1{1}.Wpi < 0.01 || stats1{2}.Wpi < 0.01
                    st = setvalue(cellid,propname,-1);   % inhibited
                elseif stats1{1}.Wpa < 0.01 || stats1{2}.Wpa < 0.01
                    st = setvalue(cellid,propname,1);   % activated
                else
                    st = setvalue(cellid,propname,0);   % non-responsive
                end
            else   % if there was one cue tone
                if stats1{1}.Wpi < 0.01
                    st = setvalue(cellid,propname,-1);   % inhibited
                elseif stats1{1}.Wpa < 0.01
                    st = setvalue(cellid,propname,1);   % activated
                else
                    st = setvalue(cellid,propname,0);   % non-responsive
                end
            end
            if ~st
                error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
            end
        end
        
    case 'rew'
        
        % Raster + PSTH aligned to feedback delivery
        alignevent = 'DeliverAllFeedback';   % trigger event
        shevent = 'StimulusOn';  % show-event
        partition = '#RewardedTrials';   % partition trials
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-1.1 0];   % baseline window for MW-test (before cue) - -3.4 was the limit
        twin = [0 0.5];   % test-window for MW-test
        
        % Add propoerty for grouping
        propname1 = 'rewardresponse2';
        if ~ismember(propname1,listtag('prop'))
            insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname1)
        end
        
        % PSTH
        for iC = 1:numCells
            cellid = pv_cells{iC};   % current cell
            stats2 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir);
            
            % Add property to CellBase
            if ~isnan(stats2{1}.Wpi) && ~isnan(stats2{1}.Wpa)
                if stats2{1}.Wpi < 0.01 && stats2{1}.Wpa > 0.01   % reward
                    st = setvalue(cellid,propname1,-1);   % inhibited
                elseif stats2{1}.Wpa < 0.01 && stats2{1}.Wpi > 0.01
                    st = setvalue(cellid,propname1,1);   % activated
                elseif stats2{1}.Wpa < 0.01 && stats2{1}.Wpi < 0.01
                    st = setvalue(cellid,propname1,2);   % activated and inhibited
                else
                    st = setvalue(cellid,propname1,0);   % non-responsive
                end
            elseif ~isnan(stats2{1}.Wpi) && isnan(stats2{1}.Wpa)
                if stats2{1}.Wpi < 0.01  % reward
                    st = setvalue(cellid,propname1,-1);   % inhibited
                else
                    st = setvalue(cellid,propname1,0);   % non-responsive
                end
            elseif isnan(stats2{1}.Wpi) && ~isnan(stats2{1}.Wpa)
                if stats2{1}.Wpa < 0.01  % reward
                    st = setvalue(cellid,propname1,1);   % activated
                else
                    st = setvalue(cellid,propname1,0);   % non-responsive
                end
            else
                st = setvalue(cellid,propname1,0);   % non-responsive
            end
            if ~st
                error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
            end
        end
        
    case 'pun'
        
        % Raster + PSTH aligned to feedback delivery
        alignevent = 'DeliverAllFeedback';   % trigger event
        shevent = 'StimulusOn';  % show-event
        partition = '#PunishedTrials';   % partition trials
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-3 -2];   % baseline window for MW-test (before cue) - -3.4 was the limit
        twin = [0 0.2];   % test-window for MW-test
        
        % Add propoerty for grouping
        propname2 = 'punishresponse';
        if ~ismember(propname2,listtag('prop'))
            insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname2)
        end
        
        % PSTH
        for iC = 1:numCells
            cellid = pv_cells{iC};   % current cell
            TE = loadcb(cellid,'TrialEvents');
            if ~all(isnan(TE.Punishment)) && isfield(TE, 'Punishment')  % if punishment was part of the session
                
                stats2 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir);
                
                % Add property to CellBase
                if ~isnan(stats2{1}.Wpi) && ~isnan(stats2{1}.Wpa)
                    if stats2{1}.Wpi < 0.01 && stats2{1}.Wpa > 0.01   % punishment
                        st = setvalue(cellid,propname2,-1);   % inhibited
                    elseif stats2{1}.Wpa < 0.01 && stats2{1}.Wpi > 0.01
                        st = setvalue(cellid,propname2,1);   % activated
                    elseif stats2{1}.Wpa < 0.01 && stats2{1}.Wpi < 0.01
                        st = setvalue(cellid,propname2,2);   % activated and inhibited
                    else
                        st = setvalue(cellid,propname2,0);   % non-responsive
                    end
                elseif ~isnan(stats2{1}.Wpi) && isnan(stats2{1}.Wpa)
                    if stats2{1}.Wpi < 0.01  % punishment
                        st = setvalue(cellid,propname2,-1);   % inhibited
                    else
                        st = setvalue(cellid,propname2,0);   % non-responsive
                    end
                elseif isnan(stats2{1}.Wpi) && ~isnan(stats2{1}.Wpa)
                    if stats2{1}.Wpa < 0.01  % punishment
                        st = setvalue(cellid,propname2,1);   % activated
                    else
                        st = setvalue(cellid,propname2,0);   % non-responsive
                    end
                else
                    st = setvalue(cellid,propname2,0);   % non-responsive
                end
                
                if ~st
                    error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
                end
            end
        end
        
    case 'pun_expectation'
        
        % Raster + PSTH aligned to feedback delivery
        alignevent = 'DeliverAllFeedback';   % trigger event
        shevent = 'StimulusOn';  % show-event
        partition = '#Punishment';   % partition trials
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-3 -2];   % baseline window for MW-test (before cue) - -3.4 was the limit
        twin = [0 0.5];   % test-window for MW-test
        
        % Add propoerty for grouping
        propname2 = 'punishresponse_expectation';
        if ~ismember(propname2,listtag('prop'))
            insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname2)
        end
        
        % PSTH
        for iC = 1:numCells
            cellid = pv_cells{iC};   % current cell
            TE = loadcb(cellid,'TrialEvents');
            if ~all(isnan(TE.Punishment)) && isfield(TE, 'Punishment')  % if punishment was part of the session
                
                stats2 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir);
                
                % Add property to CellBase
                if ~isnan(stats2{1}.Wpi) && ~isnan(stats2{1}.Wpa)
                    if stats2{1}.Wpi < 0.01 && stats2{1}.Wpa > 0.01   % punishment
                        st = setvalue(cellid,propname2,-1);   % inhibited
                    elseif stats2{1}.Wpa < 0.01 && stats2{1}.Wpi > 0.01
                        st = setvalue(cellid,propname2,1);   % activated
                    elseif stats2{1}.Wpa < 0.01 && stats2{1}.Wpi < 0.01
                        st = setvalue(cellid,propname2,2);   % activated and inhibited
                    else
                        st = setvalue(cellid,propname2,0);   % non-responsive
                    end
                elseif ~isnan(stats2{1}.Wpi) && isnan(stats2{1}.Wpa)
                    if stats2{1}.Wpi < 0.01  % punishment
                        st = setvalue(cellid,propname2,-1);   % inhibited
                    else
                        st = setvalue(cellid,propname2,0);   % non-responsive
                    end
                elseif isnan(stats2{1}.Wpi) && ~isnan(stats2{1}.Wpa)
                    if stats2{1}.Wpa < 0.01  % punishment
                        st = setvalue(cellid,propname2,1);   % activated
                    else
                        st = setvalue(cellid,propname2,0);   % non-responsive
                    end
                else
                    st = setvalue(cellid,propname2,0);   % non-responsive
                end
                
                if ~st
                    error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
                end
            end
        end
        
    case 'om'
        
        % Raster + PSTH aligned to feedback delivery
        alignevent = 'DeliverAllFeedback';   % trigger event
        shevent = 'StimulusOn';  % show-event
        partition = '#Omission';   % partition trials
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-3 -2];   % baseline window for MW-test (before cue) - -3.4 was the limit
        twin = [0 0.2];   % test-window for MW-test
        
        % Add propoerty for grouping
        propname1 = 'omissionresponse';
        if ~ismember(propname1,listtag('prop'))
            insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname1)
        end
        
        % PSTH
        for iC = 1:numCells
            cellid = pv_cells{iC};   % current cell
            stats2 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir);
            
            % Add property to CellBase
            if ~isnan(stats2{1}.Wpi) && ~isnan(stats2{1}.Wpa)
                if stats2{1}.Wpi < 0.01 && stats2{1}.Wpa > 0.01   % reward
                    st = setvalue(cellid,propname1,-1);   % inhibited
                elseif stats2{1}.Wpa < 0.01 && stats2{1}.Wpi > 0.01
                    st = setvalue(cellid,propname1,1);   % activated
                elseif stats2{1}.Wpa < 0.01 && stats2{1}.Wpi < 0.01
                    st = setvalue(cellid,propname1,2);   % activated and inhibited
                else
                    st = setvalue(cellid,propname1,0);   % non-responsive
                end
            elseif ~isnan(stats2{1}.Wpi) && isnan(stats2{1}.Wpa)
                if stats2{1}.Wpi < 0.01  % reward
                    st = setvalue(cellid,propname1,-1);   % inhibited
                else
                    st = setvalue(cellid,propname1,0);   % non-responsive
                end
            elseif isnan(stats2{1}.Wpi) && ~isnan(stats2{1}.Wpa)
                if stats2{1}.Wpa < 0.01  % reward
                    st = setvalue(cellid,propname1,1);   % activated
                else
                    st = setvalue(cellid,propname1,0);   % non-responsive
                end
            else
                st = setvalue(cellid,propname1,0);   % non-responsive
            end
            if ~st
                error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
            end
        end
end

% -------------------------------------------------------------------------
function stats1 = rasterPSTH(cellid,alignevent,shevent,partition,wn,dt,sigma,bwin,twin,issave,response_resdir)

% Raster plot and PSTH
TE = loadcb(cellid,'TrialEvents');   % load trial events
SP = loadcb(cellid,'EVENTSPIKES');
fld = fieldnames(TE);
if ~isequal(length(SP.event_stimes{1}),length(TE.(fld{1})))
    error('MATLAB:vpresponsesorter:rasterPSTH:trialMismatch',...
        'Trial number mismatch between TrialEvents and EVENTSPIKES.')
end

Raster plot
viewcell2b(cellid,'TriggerName',alignevent,'SortEvent','TrialStart','sigma',sigma,...
    'eventtype','behav','ShowEvents',{{shevent}},'Partitions',partition,'window',wn,'PSTHPlot',false);
V_handle = gcf;
maximize_figure

% Peri-event time histogram
% PSTHaxis_handle = findobj(allchild(V_handle),'type','axes','XLim',[0 1],'YLim',[0 1],'Tag','');   % handle for the empty PSTH axes
[~, ~, ~, ~, ~, stats1] = ...
    ultimate_psth(cellid,'trial',alignevent,wn,...
    'dt',dt,'sigma',sigma,'parts',partition,'isadaptive',0,...
    'maxtrialno',Inf,'baselinewin',bwin,'testwin',twin,'relative_threshold',0.1,'display',true); % calculate psth

% Plot & save
if ~iscell(stats1)
    stats1 = {stats1};   % only one PSTH
end
[~, tags] = partition_trials(TE,partition);
[mylabels, mycolors] = makeColorsLabels(@defineLabelsColors_Balazs,tags);
NumStats = length(stats1);
close_handle = nan(1,NumStats);
for iP =  1:NumStats
    Ls = findobj(allchild(stats1{iP}.axis_handle),'Type','line');   % all lines in the plot
    clr = findobj(Ls,'Type','line','Color','black');   % re-color the black PSTH
    set(clr,'Color',mycolors{iP})
    tx = findobj(allchild(stats1{iP}.axis_handle),'Type','text');   % re-position text object
    x_lim = xlim;
    y_lim = ylim;
    set(tx(1),'Position',[x_lim(1)+diff(x_lim)*0.1 (iP-1)*diff(y_lim)*0.4+y_lim(1)+diff(y_lim)*0.6 0])
    tx(1).String = regexprep(tx(1).String,'MW test',mylabels{iP});   % more meaningful labels
    if length(tx) > 1
        set(tx(2),'Position',[x_lim(1)+diff(x_lim)*0.1 (iP-1)*diff(y_lim)*0.4+y_lim(1)+diff(y_lim)*0.4  0])
        tx(2).String = regexprep(tx(2).String,'MW test',mylabels{iP});   % more meaningful labels
    end
    copyobj(tx,PSTHaxis_handle)   % copy PSTH and stat text to the raster figure
    copyobj(Ls,PSTHaxis_handle)
    xlabel(PSTHaxis_handle,['Time from ' alignevent])
    close_handle(iP) = stats1{iP}.figure_handle;
end

% Save figure
if issave
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(response_resdir,[cellidt '_' alignevent '_' partition(2:end) '.jpg']);
    saveas(V_handle,fnm)
    fnm2 = fullfile(response_resdir,[cellidt '_' alignevent '_' partition(2:end) '.fig']);
    saveas(V_handle,fnm2)
    fnm = fullfile(response_resdir,[cellidt '_' alignevent '_' partition(2:end) '.mat']);
    warning('off','MATLAB:Figure:FigureSavedToMATFile')
    save(fnm,'stats1')
end
close all
% close([V_handle, close_handle])

% -------------------------------------------------------------------------
function preproc_TE(cellid)
TE = loadcb(cellid,'TrialEvents');   % load trial events
NumTrials = TE.NTrials(1);
if ~isfield(TE, 'RewardedTrials')
    TE.RewardedTrials = nan(1,NumTrials);
    TE.RewardedTrials(~isnan(TE.Reward))=1;
    
end

if ~isfield(TE, 'PunishedTrials')
    TE.PunishedTrials = nan(1,NumTrials);
    TE.PunishedTrials(~isnan(TE.Punishment))=1;
end

[mouse, session, ~, ~] = cellid2tags(cellid);
save(fullfile(getpref('cellbase', 'datapath'), mouse, session, 'TrialEvents.mat'), '-struct', 'TE')
