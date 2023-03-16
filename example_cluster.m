function [out1, out2, out3] = example_cluster(cellid)
%EXAMPLE_CLUSTER   Example cluster of a tagged neuron.
%   [O1, O2, O3] = EXAMPLE_CLUSTER(CELLID) plots MClust projections and
%   waveforms, calculates ID and L-ratio for CELLID. Spike shape correlation
%   between light-evoked and spontaneous spikes is calculated. Waveforms of
%   other clusters in the session are also plotted. Waveforms of the
%   example cell, the tetrode neighbors and light-triggered PSTH are
%   returned in the three output arguments.
%
%   See also PLOTWAVEFORMS, PLOTALLWAVEFORMS and VIEWCELL2B.

%   Balazs Hangya
%   hangya.balazs@koki.hu
%   15-Dec-2020

% Plot MClust projections
[lim1, lim2] = findStimPeriod(cellid);   % find putative stimulated period
% HM = plot_mclust_projections2(cellid,'feature_names',{'Amplitude','Energy'},...
%     'stim_period',[lim1 lim2],'usefastplot',false);  % save the best projections individually

% Plot light-evoked and spont. spikes
out1 = plotwaveforms(cellid,'correlation',true,'evoked',true,'maxnum',5000);

% All cells of the session
[animalID, sessionID, tetrode] = cellid2tags(cellid);
cellids = findcell('rat',animalID,'session',sessionID,'tetrode',tetrode);
ID_PC = getvalue('ID_PC',cellids);
Lr_PC = getvalue('Lr_PC',cellids);
cellids = cellids(ID_PC>20&Lr_PC<0.15);

% All waveforms
NumCells = length(cellids);
out2 = struct('ID_PC',{},'Lr_PC',{},'waveforms',{},'figure',{});
for iC = 1:NumCells
    [walls, mean_all, H_all] = plotallwaveforms(cellids{iC});
    out2(iC).ID_PC = getvalue('ID_PC',cellids{iC});
    out2(iC).Lr_PC = getvalue('Lr_PC',cellids{iC});
    out2(iC).waveforms = walls;
    out2(iC).mean_waveforms = mean_all;
    out2(iC).figure = H_all;
end

% Light-triggered raster and PSTH
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
H = figure;
viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
    'FigureNum',H,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
    'EventMarkerWidth',0,'PlotZeroLine','off')
maximize_figure(H)
out3.figure = H;