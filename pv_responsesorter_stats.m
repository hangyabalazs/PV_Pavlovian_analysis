function pv_responsesorter_stats(cellids,issave,response_resdir,responsespec)
%PV_RESPONSESORTER_STATS   Peri-event time histogram.
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
for p = 1:numCells
    cell2check = pv_cells{p};   % current cell
    preproc_TE(cell2check)
end

% Raster + PSTH
switch responsespec
    case 'cue'
        
        % Raster + PSTH aligned to stimulus onset
        alignevent = 'StimulusOn';   % trigger event
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        bwin = [-2 0];   % baseline window for MW-test
        twin = [0.5 1];   % test-window for MW-test
        
        % Add propoerty for grouping
        propname = 'cueresponse_stat';
        if ~ismember(propname,listtag('prop'))
            insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname)
        end
        
        % PSTH
        for iC = 1:numCells
            cellid = pv_cells{iC};   % current cell
            
            [~, ~, ~, ~, spt ,~] = ultimate_psth(cellid, 'trials', alignevent, wn,'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',0,...
                'event_filter','none','maxtrialno',Inf,'data_type','real');
            
            % Convert window to indices
            st = abs(wn(1)) / dt;	% in ms
            nullindex = st + 1;	% index for time 0
            
            WNb = [bwin(1)/dt+nullindex bwin(2)/dt+nullindex-1];   % test window; convert to indices
            WNb = round(WNb);
            lWNb = WNb(2) - WNb(1) + 1;   % length of test window
            
            WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
            WNt = round(WNt);
            lWNt = WNt(2) - WNt(1) + 1;   % length of test window
            
            
            % Extract spikes from binraster
            b_spikes = spt(:,WNb(1):WNb(2)); % baseline spike counts
            avg_b_spikes = sum(b_spikes,2)/(lWNb)*1000; % trial-by-trial average FR
            median_b = median(avg_b_spikes); % median spike count is baseline win to decide activation or inhibition
            
            t_spikes = spt(:,WNt(1):WNt(2)); % spike counts after cue
            avg_t_spikes = sum(t_spikes,2)/(lWNt)*1000; % trial_by_trial average FR
            median_t = median(avg_t_spikes);
            
            % Statistics
            [H,p] = boxstat(avg_b_spikes, avg_t_spikes, 'baseline', 'test', 0.01, 'paired')
            set(H, 'renderer', 'painters')
            cellidt = char(cellid);
            cellidt(cellidt=='.') = '_';
            saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.fig']))
            saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.eps']))
            saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.jpg']))
            close(H)
            
            % Add property to CellBase - inhibition taken over activation!
            if (median_b < median_t) && p < 0.01
                st = setvalue(cellid,propname,1);   % activation
            elseif (median_b > median_t) && p < 0.01
                st = setvalue(cellid,propname,-1);   % inhibition
            else
                st = setvalue(cellid,propname,0);   % non-responsive
            end
        end
        if ~st
            error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
        end
        
    case 'rew'
        
        % Raster + PSTH aligned to feedback delivery
        alignevent = 'DeliverAllFeedback';   % trigger event
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        bwin = [-3 -2];   % baseline window for MW-test (before cue) - -3.4 was the limit
        twin = [0 0.5];   % test-window for MW-test
        
        % Add propoerty for grouping
        propname1 = 'rewardresponse_stat';
        if ~ismember(propname1,listtag('prop'))
            insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname1)
        end
        
        % PSTH
        for iC = 1:numCells
            cellid = pv_cells{iC};   % current cell
            [~, ~, ~, ~, spt ,~] = ultimate_psth(cellid, 'trials', alignevent, wn,'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',0,...
                'event_filter','custom','filterinput','RewardedTrials==1','maxtrialno',Inf,'data_type','real');
            
            % Convert window to indices
            st = abs(wn(1)) / dt;	% in ms
            nullindex = st + 1;	% index for time 0
            
            WNb = [bwin(1)/dt+nullindex bwin(2)/dt+nullindex-1];   % test window; convert to indices
            WNb = round(WNb);
            lWNb = WNb(2) - WNb(1) + 1;   % length of test window
            
            WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
            WNt = round(WNt);
            lWNt = WNt(2) - WNt(1) + 1;   % length of test window
            
            % Extract spikes from binraster
            b_spikes = spt(:,WNb(1):WNb(2)); % baseline spike counts
            avg_b_spikes = sum(b_spikes,2)/(lWNb)*1000; % trial-by-trial average FR
            median_b = median(avg_b_spikes); % median spike count is baseline win to decide activation or inhibition
            
            t_spikes = spt(:,WNt(1):WNt(2)); % spike counts after cue
            avg_t_spikes = sum(t_spikes,2)/(lWNt)*1000; % trial_by_trial average FR
            median_t = median(avg_t_spikes);
            
            % Statistics
            [H,p] = boxstat(avg_b_spikes, avg_t_spikes, 'baseline', 'test', 0.01,'paired')
            set(H, 'renderer', 'painters')
            cellidt = char(cellid);
            cellidt(cellidt=='.') = '_';
            saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.fig']))
            saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.eps']))
            saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.jpg']))
            close(H)
            
            % Add property to CellBase - inhibition taken over activation!
            if (median_b < median_t) && p < 0.01
                st = setvalue(cellid,propname1,1);   % activation
            elseif (median_b > median_t) && p < 0.01
                st = setvalue(cellid,propname1,-1);   % inhibition
            else
                st = setvalue(cellid,propname1,0);   % non-responsive
            end
        end
        if ~st
            error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
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
        propname2 = 'punishresponse_stat';
        if ~ismember(propname2,listtag('prop'))
            insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname2)
        end
        
        % PSTH
        for iC = 1:numCells
            cellid = pv_cells{iC};   % current cell
            TE = loadcb(cellid,'TrialEvents');
            if ~all(isnan(TE.Punishment)) && isfield(TE, 'Punishment')  % if punishment was part of the session
                
                [~, ~, ~, ~, spt ,~] = ultimate_psth(cellid, 'trials', alignevent, wn,'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',0,...
                    'event_filter','custom','filterinput','PunishedTrials==1','maxtrialno',Inf,'data_type','real');
                
                % Convert window to indices
                st = abs(wn(1)) / dt;	% in ms
                nullindex = st + 1;	% index for time 0
                
                WNb = [bwin(1)/dt+nullindex bwin(2)/dt+nullindex-1];   % test window; convert to indices
                WNb = round(WNb);
                lWNb = WNb(2) - WNb(1) + 1;   % length of test window
                
                WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
                WNt = round(WNt);
                lWNt = WNt(2) - WNt(1) + 1;   % length of test window
                
                % Extract spikes from binraster
                b_spikes = spt(:,WNb(1):WNb(2)); % baseline spike counts
                avg_b_spikes = sum(b_spikes,2)/(lWNb)*1000; % trial-by-trial average FR
                median_b = median(avg_b_spikes); % median spike count is baseline win to decide activation or inhibition
                
                t_spikes = spt(:,WNt(1):WNt(2)); % spike counts after cue
                avg_t_spikes = sum(t_spikes,2)/(lWNt)*1000; % trial_by_trial average FR
                median_t = median(avg_t_spikes);
                % Statistics
                [H,p] = boxstat(avg_b_spikes, avg_t_spikes, 'baseline', 'test', 0.01, 'paired')
                set(H, 'renderer', 'painters')
                cellidt = char(cellid);
                cellidt(cellidt=='.') = '_';
                saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.fig']))
                saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.eps']))
                saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.jpg']))
                close(H)
                
                % Add property to CellBase - inhibition taken over activation!
                if (median_b < median_t) && p < 0.01
                    st = setvalue(cellid,propname2,1);   % activation
                elseif (median_b > median_t) && p < 0.01
                    st = setvalue(cellid,propname2,-1);   % inhibition
                else
                    st = setvalue(cellid,propname2,0);   % non-responsive
                end
            end
            if ~st
                error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
            end
        end
        
    case 'om'
        
        % Raster + PSTH aligned to feedback delivery
        alignevent = 'DeliverAllFeedback';   % trigger event
        wn = [-3 3];   % full raster window in seconds
        dt = 0.001;   % raster resolution, in seconds
        sigma = 0.02;   % controls smoothing for 'spsth'
        bwin = [-3 -2];   % baseline window for MW-test (before cue) - -3.4 was the limit
        twin = [0 0.2];   % test-window for MW-test
        
        % Add propoerty for grouping
        propname3 = 'omissionresponse_stat';
        if ~ismember(propname3,listtag('prop'))
            insertdata([CELLIDLIST' num2cell(nan(size(CELLIDLIST')))],'type','property','name',propname3)
        end
        
        % PSTH
        for iC = 1:numCells
            cellid = pv_cells{iC};   % current cell
            [~, ~, ~, ~, spt ,~] = ultimate_psth(cellid, 'trials', alignevent, wn,'dt',dt,'display',false,'sigma',0.08,'parts','all','isadaptive',0,...
                'event_filter','custom','filterinput','OmittedTrials==1','maxtrialno',Inf,'data_type','real');
            
            % Convert window to indices
            st = abs(wn(1)) / dt;	% in ms
            nullindex = st + 1;	% index for time 0
            
            WNb = [bwin(1)/dt+nullindex bwin(2)/dt+nullindex-1];   % test window; convert to indices
            WNb = round(WNb);
            lWNb = WNb(2) - WNb(1) + 1;   % length of test window
            
            WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
            WNt = round(WNt);
            lWNt = WNt(2) - WNt(1) + 1;   % length of test window
            
            % Extract spikes from binraster
            b_spikes = spt(:,WNb(1):WNb(2)); % baseline spike counts
            avg_b_spikes = sum(b_spikes,2)/(lWNb)*1000; % trial-by-trial average FR
            median_b = median(avg_b_spikes); % median spike count is baseline win to decide activation or inhibition
            
            t_spikes = spt(:,WNt(1):WNt(2)); % spike counts after cue
            avg_t_spikes = sum(t_spikes,2)/(lWNt)*1000; % trial_by_trial average FR
            median_t = median(avg_t_spikes);
            
            % Statistics
            [H,p] = boxstat(avg_b_spikes, avg_t_spikes, 'baseline', 'test', 0.01, 'paired')
            set(H, 'renderer', 'painters')
            cellidt = char(cellid);
            cellidt(cellidt=='.') = '_';
            saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.fig']))
            saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.eps']))
            saveas(H, fullfile(response_resdir, [responsespec '_' cellidt '.jpg']))
            close(H)
            
            % Add property to CellBase - inhibition taken over activation!
            if (median_b < median_t) && p < 0.01
                st = setvalue(cellid,propname3,1);   % activation
            elseif (median_b > median_t) && p < 0.01
                st = setvalue(cellid,propname3,-1);   % inhibition
            else
                st = setvalue(cellid,propname3,0);   % non-responsive
            end
            if ~st
                error('MATLAB:vpresponsesorter:setvalueStatus','Error using setvalue.m')
            end
        end
end


% -------------------------------------------------------------------------
function preproc_TE(cellid)
TE = loadcb(cellid,'TrialEvents');   % load trial events
NumTrials = TE.NTrials(1);
if ~isfield(TE, 'RewardedTrials')
    TE.RewardedTrials = nan(1,NumTrials);
    TE.RewardedTrials(~isnan(TE.AllReward))=1;
    
end

if ~isfield(TE, 'PunishedTrials')
    TE.PunishedTrials = nan(1,NumTrials);
    TE.PunishedTrials(~isnan(TE.Punishment))=1;
end

if ~isfield(TE, 'OmittedTrials')
    TE.OmittedTrials = nan(1,NumTrials);
    TE.OmittedTrials(~isnan(TE.Omission))=1;
end
[mouse, session, ~, ~] = cellid2tags(cellid);
save(fullfile(getpref('cellbase', 'datapath'), mouse, session, 'TrialEvents.mat'), '-struct', 'TE')
