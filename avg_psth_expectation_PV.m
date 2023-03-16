function avg_psth_expectation_PV(pvcells, responsetype, ispart, resdir, data_type)
%AVG_PSTH_EXPECTATION_PV  Average, normalized peri-event time-histograms
%comparing neuronal response to expected and unexpected stimuli.
%
%   AVG_PSTH_EXPECTATION_PV(PVCELLS, RESPONSETYPE, ISLICK, RESDIR) calculates
%   normalized average PSTHs  of cells defined by PVCELLS, aligned to
%   RESPONSETYPE (cue, reward or punishment) and compares neuronal response
%   to expected and unexpected stimuli (timestamps extracted from
%   TrialEvents). If ISPART is true, trials are partitioned by the defined
%   variable and plotted separately. Cells with >5 trials are used for
%   the PSTH. Average firing rate and maximal response in the response
%   window are compared as well. The PSTHs and box-whiskers plots are saved
%   to RESDIR.
%
%   AVG_PSTH_EXPECTATION_PV(PVCELLS, RESPONSETYPE, ISLICK, RESDIR, DATA_TYPE)
%   calls ULTIMATE_PSTH for virtual (simulated) spikes if DATA_TYPE input
%   argument is 'virtual'. This can be used to analyze simulated data the
%   same way as real spiking.
%
%   See also ULTIMATE_PSTH and BOXSTAT.

%   Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   01-May-2020

%   Code review: 5/11/20, 6/9/20, 7/30/20, 11/26/20

% Input arguments
if nargin < 1 || isempty(pvcells)
    pvcells = select_pv_cells(pvcells);
end
NumCell = length(pvcells);
if nargin < 5
    data_type = 'real';  % control simulated data
end

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Define filterinput
switch responsetype
    case 'cue' % cure response
        filters = {'TrialType==1' 'TrialType==2'};
        align = 'StimulusOn';
        twin = [0 0.5];
        bwin = [-1 0];   % baseline window
    case 'lickcue' % cue reaponse with anticipatory licking
        filters = {'Lick==1&TrialType==1' 'Lick==1&TrialType==2'};
        align = 'StimulusOn';
        twin = [0 0.5];
        bwin = [-1 0];   % baseline window
    case 'reward' % reward response
        filters = {'AllReward==1' 'AllReward==2'};
        align = 'DeliverAllFeedback';
        twin = [0 0.2];
        bwin = [-3 -2];   % baseline window
    case 'lickreward' % reward response after anticipatory licking
        filters = {'Reward==1' 'Reward==2'};
        align = 'DeliverAllFeedback';
        twin = [0 0.2];
        bwin = [-3 -2];   % baseline window
    case 'lickpunish' % lick & punish response
        filters = {'FalseAlarm==1&TrialType==1' 'FalseAlarm==1&TrialType==2'};
        align = 'DeliverAllFeedback';
        twin = [0 0.2];
        bwin = [-3 -2];   % baseline window
    case 'punish' % punishment response
        filters = {'Punishment==1' 'Punishment==2'};
        align = 'DeliverAllFeedback';
        twin = [0 0.2];
        bwin = [-3 -2];   % baseline window
    case 'omission'  % omission responses
        filters = {'Omission==1' 'Omission==2'};  % omission after likely reward (1) or likely punishment cue (2)
        align = 'DeliverAllFeedback';
        twin = [0 0.5];
        bwin = [-3 -2];   % baseline window
    case 'mediansplit'  % omission responses
        filters = {'MedialSplit==1' 'MedialSplit==2'};  % omission after likely reward (1) or likely punishment cue (2)
        align = 'DeliverAllFeedback';
        twin = [0 0.5];
        bwin = [-3 -2];   % baseline window
end

% Filter trials for anticipatory lick response
switch ispart
    case 'lick'
        part = '#Lick';
        evtype = 'trial';
    case 'cue_corr'
        part = '#Cueresp';
        evtype = 'lick';
        twin = [0 1];  % licks can be counted in the entire period before reinforcement
        bwin = [-2 0];   % baseline window
    case 'cue'
        part = '#Cueresp';
        evtype = 'trial';
    case 'previous_outcome'
        part = '#PreviousOutcome';
        evtype = 'trial';
    case 'none'
        part = 'all';
        evtype = 'trial';
end

% Time window
wn = [-4 4];
dt = 0.001;
time = wn(1)*1000:dt*1000:wn(2)*1000;   % time vector

% Index for time 0
st = abs(wn(1)) / dt;   % in ms
nullindex = st + 1;

% PSTH
R1 = runanalysis(@ultimate_psth,...
    evtype, align, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts',part,'isadaptive',0,...
    'event_filter','custom', 'filterinput',filters{1},'maxtrialno',Inf,'data_type',data_type,...
    'baselinewin',bwin,'testwin',twin, 'forcesmoothedstat', true, 'relative_threshold',0.1,'cellids',pvcells);
R2 = runanalysis(@ultimate_psth,...
    evtype, align, wn,...
    'dt',dt,'display',false,'sigma',0.08,'parts',part,'isadaptive',0,...
    'event_filter','custom', 'filterinput',filters{2},'maxtrialno',Inf,'data_type',data_type,...
    'baselinewin',bwin,'testwin',twin, 'forcesmoothedstat', true, 'relative_threshold',0.1,'cellids',pvcells);

% Window for testing the potential effect
WNt = [twin(1)/dt+nullindex twin(2)/dt+nullindex-1];   % test window; convert to indices
WNt = round(WNt);
lWNt = WNt(2) - WNt(1) + 1;   % length of test window

WNb = [bwin(1)/dt+nullindex bwin(2)/dt+nullindex-1];   % test window; convert to indices
WNb = round(WNb);
lWNb = WNb(2) - WNb(1) + 1;   % length of test window

BR1 = R1(:,5); % binraster for cue1
BR2 = R2(:,5); % binraster for cue2
STAT1 = R1(:,6);
STAT2 = R2(:,6);

avg_spikecount1 = nan(size(R1{1,1},1), size(BR1 ,1)); % preallocate space for average response and maxvalue
avg_spikecount2 = nan(size(R2{1,1},1), size(BR2 ,1));
bsl_spikecount1 = nan(size(R1{1,1},1), size(BR1 ,1)); % preallocate space for average baseline firing
bsl_spikecount2 = nan(size(R2{1,1},1), size(BR2 ,1));
maxvalue1 = nan(size(R1{1,1},1), size(BR1 ,1));
maxvalue2 = nan(size(R2{1,1},1), size(BR2 ,1));

H = figure;
if size(R1{1,1},1) > 1 % if there are more than one partitions, plot them on different figures
    G = figure;
end

NumParts = size(R1{1,1},1); % number of partitions
for p = 1:NumParts % if there are more than one partitions
    
    % Set colors for plotting
    if size(R1{1,1},1) == 1 % no partition
        linecolor1 = [0 0.9 0] ;
        shadecolor1 = [0 0.9 0] ;
        linecolor2 = [0.5 0 0];
        shadecolor2 = [0.5 0 0];
    else
        if p == 1
            linecolor1 = [0 0.9 1] ;
            shadecolor1 = [0 0.9 1] ;
            linecolor2 = [0 0.9 1];
            shadecolor2 = [0 0.9 1];
        elseif p == 2
            linecolor1 = [0.5 0 1];
            shadecolor1 = [0.5 0 1];
            linecolor2 = [0.5 0 1];
            shadecolor2 = [0.5 0 1];
        elseif p == 3
            linecolor1 = [0.1 0 1];
            shadecolor1 = [0.1 0 1];
            linecolor2 = [0.1 0 1];
            shadecolor2 = [0.1 0 1];
        end
    end
    
    [psth_R1, psth_R2] = deal(nan(NumCell,size(R1{1,1},2)));
    baseline_inx = find(time==bwin(1)*1000):find(time==bwin(2)*1000); % baseline indices
    Rinx = 2;
    for iC = 1:NumCell % Z-score by baseline
        
        % set trial selection
        if NumParts ==1 % if there is one partition only
            T1size = size(R1{iC,5},1);
            T2size = size(R2{iC,5},1);
        else % two partitions
            T1size = size(R1{iC,5}{2},1); % 'no lick' trials for reward cue
            T2size = size(R2{iC,5}{1},1); % 'lick' trials for punish cue
        end
        
        if (T1size>5) && (T2size>5) %  only use cells with >5 trials in each partition
            
            mn = mean(R1{iC,Rinx}(1,baseline_inx)); % means and sds calculated for the average PSTH
            sd = std(R1{iC,Rinx}(1,baseline_inx));
            
            if NumParts == 1 % if there is one partition - compared plots are normalized with sd and mean calculated from R1
                mn2 = mn;
                sd2 = sd;
            else
                mn2 = mean(R2{iC,Rinx}(1,baseline_inx));  % p == 1 for the normalizing constants, to preserve partition differences
                sd2 = std(R2{iC,Rinx}(1,baseline_inx));
            end
            if strcmp(evtype,'lick')
                psth_R1(iC,:) = R1{iC,Rinx}(p,:);   % do not normalize lick raster
                psth_R2(iC,:) = R2{iC,Rinx}(p,:);
            else
                psth_R1(iC,:) = (R1{iC,Rinx}(p,:) - mn) / sd;   % the same mean and SD used within a cell,
                psth_R2(iC,:) = (R2{iC,Rinx}(p,:) - mn2) / sd2;
            end
            
            % Binraster and statistics
            if NumParts > 1 % if more than one partition
                binrast1 = BR1{iC}{p};
                binrast2 = BR2{iC}{p};
                
                maxvalue1(p,iC) = STAT1{iC}{p}.maxvalue; % maxvalues
                maxvalue2(p,iC) = STAT2{iC}{p}.maxvalue;
            else
                binrast1 = BR1{iC};
                binrast2 = BR2{iC};
                
                maxvalue1(p,iC) = STAT1{iC}.maxvalue; % maxvalues
                maxvalue2(p,iC) = STAT2{iC}.maxvalue;
            end
            
            avg_spikecount1_cell = nan(1,size(binrast1 ,1)); % calculate average spike counts from binraster
            avg_spikecount2_cell = nan(1,size(binrast2 ,1));
            
            bsl_spikecount1_cell = nan(1,size(binrast1 ,1)); % calculate average spike counts from binraster
            bsl_spikecount2_cell = nan(1,size(binrast2 ,1));
            
            for j = 1:size(binrast1 ,1)   % loop through trials
                avg_spikecount1_cell(j) = sum(binrast1(j, WNt(1):WNt(2)), 2) / lWNt / dt;
                bsl_spikecount1_cell(j) = sum(binrast1(j, WNb(1):WNb(2)), 2) / lWNb / dt;
            end
            avg_spikecount1(p,iC) = mean(avg_spikecount1_cell);
            bsl_spikecount1(p,iC) = mean(bsl_spikecount1_cell);
            
            for k = 1:size(binrast2 ,1)   % loop through trials
                avg_spikecount2_cell(k) = sum(binrast2(k, WNt(1):WNt(2)), 2) / lWNt / dt;
                bsl_spikecount2_cell(k) = sum(binrast2(k, WNb(1):WNb(2)), 2) / lWNb / dt;
            end
            avg_spikecount2(p,iC) = mean(avg_spikecount2_cell);
            bsl_spikecount2(p,iC) = mean(bsl_spikecount2_cell);
        end
    end
    
    psth_R1(any(isnan(psth_R1), 2), :) = []; % remove NaNs
    psth_R2(any(isnan(psth_R2), 2), :) = [];
%     
%     maxvalue1 = maxvalue1-bsl_spikecount1;
%     maxvalue2 = maxvalue2-bsl_spikecount2;

    % Plot
    inx = 1:size((psth_R1),1);  % all
    inx2 = 1:size((psth_R2),1);  % all
    disp(length(inx))
    
    figure(H) % plot average PSTHs
    if NumParts > 1 && isvarname('filters')
        title(filters{1})
    end
    errorshade(time,mean(psth_R1(inx,:)),std(psth_R1(inx,:))/sqrt(size(psth_R1(inx,:),1)),...
        'LineColor',linecolor1,'ShadeColor',shadecolor1)
    set(gca,'XLim',[-3000 3000]);
    set(gca,'YLim',[-(min(mean(psth_R1(inx,:)))+5) max(mean(psth_R1(inx,:)))+5]);
    
    if size(R1{1,1},1) > 1 % if there are more than one partitions, plot them on different figures
        figure(G)
        title(filters{2})
    end
    errorshade(time,mean(psth_R2(inx2,:)),std(psth_R2(inx2,:))/sqrt(size(psth_R2(inx2,:),1)),...
        'LineColor',linecolor2,'ShadeColor',shadecolor2)
    set(gca,'XLim',[-3000 3000]);
    set(gca,'YLim',[-(min(mean(psth_R2(inx,:)))+10) max(mean(psth_R2(inx,:)))+10]);
end

% Save figure
filename = fullfile(resdir,['compare_expectation_' responsetype '.fig']);
saveas(H,filename);
set(H,'renderer', 'painters')
filename2 = fullfile(resdir,['compare_expectation_' responsetype '.eps']);
saveas(H,filename2);
filename3 = fullfile(resdir,['compare_expectation_' responsetype '.jpg']);
saveas(H,filename3);
close(H)

if size(R1{1,1},1) > 1 % if there are more than one partitions, plot them on different figures
    filename = fullfile(resdir,['compare_expectation_' responsetype '_2.fig']);
    saveas(G,filename);
    set(G, 'renderer', 'painters')
    filename = fullfile(resdir,['compare_expectation_' responsetype '_2.eps']);
    saveas(G,filename);
    filename = fullfile(resdir,['compare_expectation_' responsetype '_2.jpg']);
    saveas(G,filename);
end
close all

avg_spikecount1(:,any(isnan(avg_spikecount1(p,:)), 1)) = []; % remove NaNs
avg_spikecount2(:,any(isnan(avg_spikecount2(p,:)), 1)) = [];


if size(R1{1,1},1) <= 1 % when trials are not partitioned
    
    % Box-whisker plot
    [H, Wp1] = boxstat(avg_spikecount1,avg_spikecount2,filters{1},filters{2},0.05,'paired');
    [G, Wp2] = boxstat(maxvalue1,maxvalue2,filters{1},filters{2},0.05,'paired');
    [L, Wp3] = boxstat(bsl_spikecount1,bsl_spikecount2,filters{1},filters{2},0.05,'paired');
    
    % Boxplot for max value difference
    figure;
    boxplot(maxvalue1-maxvalue2);
    hold on;
    %     plot(ones(size(maxvalue1))+(rand(size(maxvalue1))-0.5)/2,maxvalue1-maxvalue2,'.')
    y = maxvalue1 - maxvalue2;
    x = ones(length(y),1);
    scatter(x(:),y(:),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.05);
    I = gcf;
    
    figure;
    boxplot(avg_spikecount1-avg_spikecount2);
    hold on;
    %     plot(ones(size(maxvalue1))+(rand(size(maxvalue1))-0.5)/2,maxvalue1-maxvalue2,'.')
    y = avg_spikecount1 - avg_spikecount2;
    x = ones(length(y),1);
    scatter(x(:),y(:),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.05);
    J = gcf;
    
    figure;
    boxplot(bsl_spikecount1-bsl_spikecount2);
    hold on;
    y = bsl_spikecount1 - bsl_spikecount2;
    x = ones(length(y),1);
    scatter(x(:),y(:),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.05);
    K = gcf;
    
    % Save figure
    filename = fullfile(resdir,['compare_expectation_' responsetype  '_boxstat.fig']);
    saveas(H,filename);
    
    filename = fullfile(resdir,['compare_maxvalue_' responsetype  '_boxstat.fig']);
    saveas(G,filename);
    
    filename = fullfile(resdir,['compare_baseline_' responsetype  '_boxstat.fig']);
    saveas(L,filename);
    
    filename = fullfile(resdir,['compare_maxdifference_' responsetype  '_boxstat.fig']);
    saveas(I,filename);
    
    filename = fullfile(resdir,['compare_avgdifference_' responsetype  '_boxstat.fig']);
    saveas(J,filename);
    
    filename = fullfile(resdir,['compare_bsldifference_' responsetype  '_boxstat.fig']);
    saveas(K,filename);
    
    % Save variables
    filename = fullfile(resdir,['compare_FR_maxvalue_' responsetype  '.mat']);
    save(filename,'avg_spikecount1','avg_spikecount2','maxvalue1','maxvalue2');
    
else
    
    % Box-whisker plot
    [H1, Wp1] = boxstat(avg_spikecount1(1,:),avg_spikecount1(2,:),[filters{1} ' ' part ' 1'],[filters{1} ' ' part ' 2'],0.05,'paired');
    [H2, Wp1] = boxstat(avg_spikecount2(1,:),avg_spikecount2(2,:),[filters{2} ' ' part ' 1'],[filters{2} ' ' part ' 2'],0.05,'paired');
    
    % Save figure
    filename = fullfile(resdir,['compare_expectation_' responsetype  '_boxstat_T1.fig']);
    filename2 = fullfile(resdir,['compare_expectation_' responsetype  '_boxstat_T2.fig']);
    
    saveas(H1,filename);
    saveas(H2,filename2);
    
    [G1, Wp2] = boxstat(maxvalue1(1,:),maxvalue1(2,:),[filters{1} ' ' part ' 1'],[filters{1} ' ' part ' 2'],0.05,'paired');
    [G2, Wp2] = boxstat(maxvalue2(1,:),maxvalue2(2,:),[filters{2} ' ' part ' 1'],[filters{2} ' ' part ' 2'],0.05,'paired');
    figure;
    boxplot(maxvalue1(1,:)-maxvalue1(2,:));
    I1 = gcf;
    figure;
    boxplot(maxvalue2(1,:)-maxvalue2(2,:));
    I2 = gcf;
    
    % Save figure
    filename3 = fullfile(resdir,['compare_maxvalue_' responsetype  '_boxstat_T1.fig']);
    filename4 = fullfile(resdir,['compare_maxvalue_' responsetype  '_boxstat_T2.fig']);
    filename5 = fullfile(resdir,['compare_maxdifference_' responsetype  '_boxstat_T1.fig']);
    filename6 = fullfile(resdir,['compare_maxdifference_' responsetype  '_boxstat_T2.fig']);
    
    saveas(G1,filename3);
    saveas(G2,filename4);
    saveas(I1,filename5);
    saveas(I2,filename6);
    
    % Save variables
    filename = fullfile(resdir,['compare_FR_maxvalue_' responsetype  '.mat']);
    save(filename,'avg_spikecount1','avg_spikecount2','maxvalue1','maxvalue2');
    
    %     close all
end