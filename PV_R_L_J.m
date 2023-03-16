function PV_R_L_J(cellids, resdir, spec, issave)
%PV_R_L_J(CELLIDS, RESDIR, ISSAVE) Latency and jitter CDF of event evoked
%spiking.
%   PV_R_L_J(CELLIDS,RESDIR, ISSAVE) calculates latency and jitter for
%   event evoked spiking of neurons listed in CELLIDS. If CELLIDS is and
%   empy cell, parvalbumin+ neurons of the current cellbase are
%   automatically selected. SPEC defines the need for light ('tag'), cue
%   ('cue') or reinforcement ('rew' or 'pun') aligned spiking. Results,
%   latency CDF and jitter CDF are saved to RESDIR.
%
%   See also RELIABILITY_LATENCY_JITTER

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   16-Apr-2020

%   Code review: BH 6/7/20

% Input check
narginchk(0,4)
if nargin < 1
    cellids = [];
end

% Choose ChAT cells
if isempty(cellids)
    cellids = select_pv_cells(cellids); % select PV+ cells from CELLIDLIST
end

if nargin < 2
    resdir = fullfile(getpref('cellbase', 'datapath'), 'R_L_J\'); % results directory
end

% Results directory
if ~isfolder(resdir)
    mkdir(resdir)
end

% Choose filters
switch spec
    case 'tag'
        evtype = 'stim';
        alignevent = 'PulseOn';
        evfilter = 'none';
        finput = {[]};
        win = [-0.005 0.01];
        baselinewin = [-0.005 0];
        testwin = [0 0.01];
        reltr = 0.5;
    case 'cue'
        evtype = 'trial';
        alignevent = 'StimulusOn';
        evfilter = 'custom';
        finput = {'TrialType==1' 'TrialType==2'};
        win = [-1 1];
        baselinewin = [-1 0];
        testwin = [0 1];
        reltr = 0.05;
    case 'rew'
        evtype = 'trial';
        alignevent = 'DeliverAllFeedback';
        evfilter = 'custom';
        finput = {'AllReward==1' 'AllReward==2'};
        win = [-0.02 0.15];
        baselinewin = [-0.02 0];
        testwin = [0 0.15];
        reltr = 0.05;
    case 'pun'
        evtype = 'trial';
        alignevent = 'DeliverAllFeedback';
        evfilter = 'custom';
        finput = {'Punishment==1' 'Punishment==2'};
        win = [-0.02 0.15];
        baselinewin = [-0.02 0];
        testwin = [0 0.15];
        reltr = 0.05;
end

% Reliability, latency, jitter
NumCell = length(cellids); % number of cells
[Reliability Latency Jitter] = deal(nan(1,NumCell)); % preallocate space
SpikeNumberDistribution = cell(1,NumCell);

for k = 1:length(finput)
    for iC = 1:NumCell
        cellid = cellids{iC};
        
        % Reliability, latency, jitter
        [~, latency, jitter, ~, ~, ~, ~, ~, ~] = ...
        reliability_latency_jitter(cellid,...
        'event_type',evtype,'event',alignevent,'window',win,...
        'event_filter',evfilter,'filterinput',finput{k},'isadaptive',0,...
        'baselinewin',baselinewin,'testwin',testwin,'relative_threshold',reltr,...  % 'jitterdefinition','burst',...
        'display',true);
        Latency(iC) = latency;
        Jitter(iC) = jitter;
        close all
    end
    
    % Plot and save CDF
    if issave
        
        % Latency CDF
        blue = [0 153 255] / 255;
        figure
        edges = 0:0.0001:max(Latency)+0.01;   % bin edges
        dist_latency = histc(Latency,edges);   % latency histogram
        dist_latency = [0 dist_latency(1:end-1)];   % values corresponding to the edges
        dist_latency = dist_latency / sum(dist_latency);   % normalize
        stairs(edges,cumsum(dist_latency),'Color',blue,'LineWidth',2)
        axis tight
        set(gca,'box','off','FontSize',12,'TickDir','out')
        xlabel('Latency')
        setmyplot_balazs
        set(gcf, 'Renderer', 'painters')
        if ~isempty(finput{1})
            [filt, trialtype] = strtok(finput{k}, '==');
        else
            filt = 'lightstim';
            trialtype = '_';
        end
        saveas(gcf, fullfile(resdir, [filt '_' trialtype(end) '_Latency_CDF.fig'])) % save
        saveas(gcf, fullfile(resdir, [filt '_' trialtype(end) '_Latency_CDF.eps'])) % save
        
        % Jitter CDF
        figure
        edges = 0:0.0001:max(Jitter)+0.01;   % bin edges
        dist_jitter = histc(Jitter,edges);   % jitter histogram
        dist_jitter = [0 dist_jitter(1:end-1)];   % values corresponding to the edges
        dist_jitter = dist_jitter / sum(dist_jitter);   % normalize
        stairs(edges,cumsum(dist_jitter),'Color',blue,'LineWidth',2)
        axis tight
        set(gca,'box','off','FontSize',12,'TickDir','out')
        xlabel('Jitter')
        setmyplot_balazs
        set(gcf, 'Renderer', 'painters')
        saveas(gcf, fullfile(resdir, [filt '_' trialtype(end) '_Jitter_CDF.fig'])) % save
        saveas(gcf, fullfile(resdir, [filt '_' trialtype(end) '_Jitter_CDF.eps'])) % save
        close all
    end
end