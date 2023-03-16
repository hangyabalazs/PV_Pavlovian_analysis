function PV_lightpsth_colored(cellids, impdir, resdir)
%   PV_LIGHTPSTH_COLORED Population PSTH of light-aligned PSTH.
%   PV_LIGHTPSTH_COLORED(CELLIDS, IMPDIR, RESDIR) calculates adaptive PSTH 
%   (see DAPSTH) aligned to light pulse onset using all pulses of the most
%   efficient stimulation frequency. PSTHs are normalized by peak value.
%
%   See also ULTIMATE_PSTH and DAPSTH.

% Reliability data
if isempty(impdir)
    impdir = fullfile(getpref('cellbase', 'datapath'), 'taggingsummary');
end

% Results directory
if ~isfolder(resdir)
    mkdir(resdir)
end

NumCells = length(cellids);   % number of cholinergic cells

% PSTH
wn = [-0.005 0.01];   % window
dt = 0.0002;   % resolution (1 ms)
time = wn(1):dt:wn(2);   % time vector
[psths hsts psths_norm hsts_norm] = deal(nan(NumCells,length(time)));
for iC = 1:NumCells
    cellid = cellids{iC};   % current cell
    
    % Light stim. frequencies
    SE = loadcb(cellid,'StimEvents');
    burst_types = sort(unique(SE.BurstNPulse),'ascend');
    disp(burst_types)
    
    % Reliability data
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(impdir,['taggingsummary' cellidt '_TAGSUM.mat']);
    load(fnm)
    [~, m2] = max(reliability);
    bnp = burst_types(m2);
    disp(bnp);
    
    % PSTH
    [psth, spsth, spsth_se, ~, spt] = ...
        ultimate_psth(cellid,'stim','PulseOn',wn,...
        'dt',dt,'display',true,'sigma',0.001,'parts','all','isadaptive',2,...
        'event_filter','custom','filterinput',['BurstNPulse==' num2str(bnp)],'maxtrialno',Inf);   % PSTH
    sspt = sum(spt);
    psths(iC,:) = psth;  % adaptive psth
    hsts(iC,:) = sspt;  % summed biraster
    psths_norm(iC,:) = zscore(psth);  % normalized adaptive psth
    hsts_norm(iC,:) = sspt / max(sspt);  % normalized summed biraster
    
    % Plot
    figure
    bar(time,sspt,'FaceColor','k','EdgeColor','k','BarWidth',1)
    xlim(wn)
end

% Plot summary
[m1 m2] = max(hsts_norm(:,25:50),[],2);
[srt Ia] = sort(m2,'ascend');   % sort based on FA response
figure   % plot all FA PSTHs, sorted
imagesc(time,1:NumCells,hsts_norm(Ia,:))
colormap(hot)
saveas(gcf,fullfile(resdir,'lightpoppsth_sorted.fig'))