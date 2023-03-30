function viewphotometry_avg_PV(cellids, root, resdir, varargin)
%VIEWPHOTOMETRY_AVG   Average PSTH of fiber photometry data.
%   VIEWPHOTOMETRY_AVG(CELLIDS, RESDIR, VARARGIN) creates an average PSTH
%   of fiber photometry data recorded from sessions listed in CELLIDS and located
%   in ROOT. PETHs are saved to RESDIR.

%   See also VIEWPHOTOMETRY.M, PERIEVENT.M

%   Balint Kiraly, Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu

% Default arguments
default_args={...
    'window',               [-3 3];...
    'dt',                   0.01;...
    'sigma',                10;...
    'isadaptive'            false;...
    'FigureNum',            1;...
    'Signal'                'dff' %s465, s405, dff
    'TriggerEvent',         'DeliverAllFeedback';...
    'SortEvent',            'TrialStart';...
    'ShowEvents',           {{'StimulusOn'}};...
    'ShowEventsColors',     {{[0 0.8 0] [0.8 0.8 0] [0 0.8 0.8]}};...
    'Num2Plot',             'all';...
    'PlotDashedEvent',      '';...
    'PlotDashedCondition',  'min';...
    'PSTHPlot',             1;...
    'PSTHlinewidth',        1.5;...
    'DashedLineStyle',      ':';...
    'LastEvents',           '';...
    'Partitions',           'all';...
    'PrintCellID',          'on';...
    'PrintCellIDPos',       'bottom-right';...
    'BurstPSTH'             'off';......
    };
[g,error] = parse_args(default_args,varargin{:});


% Input argument check
if nargin < 1
    error('Not enough input arguments. Please provide, at least, an animal and a session.');
end

isnorm = 0;

if ~isfolder(resdir)
    mkdir(resdir)
end

if strcmp(g.Partitions, '#AllReward') % baseline windows used for normalization
    bwin = [-3 -2];
    twin = 0.6;
    color1 = [0.22 0.78 0.62];
    color2 = [0 0.4 0.4];
elseif strcmp(g.Partitions, '#Reward') % baseline windows used for normalization
    bwin = [-3 -2];
    twin = 0.7;
    color1 = [0.22 0.78 0.62];
elseif strcmp(g.Partitions, '#PunishedTrials')
    bwin = [-3 -2];
    twin = 1;
    color1 = [0.2 0 0];
    color2 = [1 0.6 0.8];
elseif strcmp(g.Partitions, '#TrialType')
    bwin = [-1 0];
    twin = 0.5;
    color1 = [0.2 0.6 0.6];
    color2 = [0.4 0 0.2];
end

% Preallocate
[psths1, psths1_st, psths2, psths2_st]=deal(nan(length(cellids),72289)); 
[psths_isos1, psths_isos2]=deal(nan(length(cellids),72289)); 
[mean1, mean2, mean1_st, mean2_st, maxval1, maxval2]=deal(nan(1,size(cellids,2)));

% Loop through sessions
for c = 1:length(cellids)
    current_cellid = cellids{c};
    animalID = current_cellid{1};
    sessionID = current_cellid{2};
    channel = current_cellid{3};
    
    if strcmp(channel, 'Ch1')
        g.Signal = 'dff_D';
        g.isosSignal1 = 's465_D';
        g.isosSignal2 = 's405_D';
    elseif strcmp(channel, 'Ch2')
        g.Signal = 'dff_A';
        g.isosSignal1 = 's465_A';
        g.isosSignal2 = 's405_A';
    end
    
    % Load data and FiberEvents
    path = [root filesep animalID filesep sessionID filesep];
    TE = load([path 'FiberEvents.mat']);
    TE.PunishedTrials = nan(1,length(TE.NTrials));
    TE.PunishedTrials(~isnan(TE.Punishment))=1;
    DATA = load([path filesep 'proF.mat']);
    TE.Blocknum(TE.Hit~=1)=NaN;
    
    % Extracting valid trials
    [COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);
    vinx = cellfun(@(s)(~isempty(s)),COMPTRIALS);
    COMPTRIALS = COMPTRIALS(vinx);
    
    TAGS = TAGS(vinx);
    trigev = TE.(g.TriggerEvent);
    if ~iscell(trigev)
        valid_trials = find(~isnan(trigev));
    else
        valid_trials = find(cellfun(@(s)~isempty(s),trigev));
    end
    
    % Creating time vector
    time = g.window(1):(1/DATA.sr):g.window(2);
    b_inx=find((time>=bwin(1))&time<bwin(2)); % baseline window indices
    
    % Sort trials
    NUMevents = length(g.SortEvent);
    if iscellstr(g.SortEvent)
        sort_var = nan(NUMevents,NUMtrials);
        for iS = 1:NUMevents
            sort_var(iS,:) = TE.(g.SortEvent{iS}) - TE.(g.TriggerEvent);
        end
        sort_var = min(sort_var);
    elseif ~isempty(g.SortEvent)
        if ~iscell(TE.(g.TriggerEvent))
            sort_var = TE.(g.SortEvent) - TE.(g.TriggerEvent);
        else
            gte = nan(1,NUMtrials);
            inx = ~cellfun(@isempty,TE.(g.TriggerEvent));
            gte(inx) = cell2mat(cellfun(@(s)s(1),TE.(g.TriggerEvent)(inx),'UniformOutput',false));
            sort_var = TE.(g.SortEvent) - gte;
        end
    else
        sort_var = NaN;
    end
    
    [mylabels, mycolors, mycolors2,mylinestyle] = makeColorsLabels(@defineLabelsColors_Balazs,TAGS);
    
    
    % Peri-event matrix for the given contingencies
    for iPAR = 1:length(TAGS)
        valT = intersect(valid_trials,COMPTRIALS{1,iPAR});
        sortT = TE.(g.TriggerEvent)(valT);
        vdisc = [];
        if strcmp(g.TriggerEvent,'TrialStart')
            [~,vdisc] = min(abs((sortT) - (DATA.tss-DATA.tss(1))));
        else
            startT = TE.TrialStart(valT);
            %vdisc = round((startT + sortT) * DATA.sr);
            for i= 1:length(startT)
                [~,vdisc(i)] = min(abs((startT(i) + sortT(i)) - (DATA.tss-DATA.tss(1))));
            end
        end
        
        nT = [1:1:length(vdisc)];
        
        [fibmean,fibcv,spmat] = perievent(vdisc,DATA.(g.Signal),DATA.sr,g.window,isnorm,g.dt,g.sigma); % peri-event
        [~,~,spmat_isos1] = perievent(vdisc,DATA.(g.isosSignal1),DATA.sr,g.window,isnorm,g.dt,g.sigma); % peri-event
        [~,~,spmat_isos2] = perievent(vdisc,DATA.(g.isosSignal2),DATA.sr,g.window,isnorm,g.dt,g.sigma); % peri-event

    
        [~,inx]=sort(sort_var(valT));
        
        if max(inx)>size(spmat,1)
            delinx = find(inx>size(spmat,1))
            inx(delinx)=[];
        end
        
        % Remove big artifacts
        spmat_avg = mean(spmat, 1);
        spmat_std = std(spmat);
        boundary_up = spmat_avg+2*spmat_std;
        boundary_low = spmat_avg-2*spmat_std;
        
        for j = 1:size(inx,2)
            all_points = length(boundary_up);
            ten_pct = all_points/10;
            deviation_up=find(spmat(inx(j),:)>boundary_up);
            deviation_low = find(spmat(inx(j),:)<boundary_low);
            if (length(deviation_up)+length(deviation_low))>ten_pct
                inx(j)=nan;
            end
        end
        inx(isnan(inx))=[]; % Only keep low noise trials
        
        if iPAR==1
            psths1(c,:) = nanmean(spmat(inx,:),1);
            mn = nanmean(psths1(c,b_inx));
            sd = nanstd(psths1(c,b_inx));
            psths1_st(c,:) = (psths1(c,:)-mn)/sd;
            
            psths_isos1(c,:)= nanmean(spmat_isos1(inx,:),1);
            psths_isos2(c,:)= nanmean(spmat_isos2(inx,:),1);

        elseif iPAR==2
            psths2(c,:) = nanmean(spmat(inx,:),1);
            psths2_st(c,:) = (psths2(c,:)-mn)/sd; % use the same mean and sd for normalization
        end
        
        w1_inx = min(find(time>0));
        w2_inx = min(find(time>twin));
    end
end

tag = g.Partitions;

% Plot
figure;
errorshade(time, nanmean(psths1_st,1), nanstd(psths1_st,1)/sqrt(size(psths1_st,1)),...
    'LineColor',color1,'ShadeColor',color1)

% Save figure
fnm = fullfile(resdir,[tag '.fig']);
saveas(gcf,fnm)
set(gcf,'renderer', 'painters')
fnm = fullfile(resdir,[tag '.eps']);
saveas(gcf,fnm)
fnm = fullfile(resdir,[tag '.bmp']);
saveas(gcf,fnm)
set(gca,'Xlim', [-1 1])

close(gcf)

% Plot isosbestic signals
figure;
errorshade(time, nanmean(psths_isos1,1), nanstd(psths_isos1,1)/sqrt(size(psths_isos1,1)),...
    'LineColor',color1,'ShadeColor',color1)
hold on
errorshade(time, nanmean(psths_isos2,1), nanstd(psths_isos2,1)/sqrt(size(psths_isos2,1)),...
    'LineColor',color2,'ShadeColor',color2)
hold off
set(gca,'Xlim', [-1 1])

% Save plots
fnm = fullfile(resdir,[tag 'isosbestic.fig']);
saveas(gcf,fnm)
set(gcf,'renderer', 'painters')
fnm = fullfile(resdir,[tag 'isosbestic.eps']);
saveas(gcf,fnm)
fnm = fullfile(resdir,[tag 'isosbestic.bmp']);
saveas(gcf,fnm)
