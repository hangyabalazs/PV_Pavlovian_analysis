function Pop_PETH_Cluster(varargin)
%POP_PETH_CLUSTER Searhcing for clusters of similarly behaving neurons.
%   POP_PETH_CLUSTER Finds groups of similary behaving cells based on the
%   k-means clustering performed in the space spanned by the first
%   PCA_dim(input) principal components of the normalized PETHs of the
%   neurons. See viewcell2b (CELLBASE) for descripton of the optional
%   inpusts. Optionally (defined in the input usedoubletrigger, default is
%   yes) 2 trigger events can be used, in this case trial by trial
%   variabilites in the time between the two triggers is equalized by
%   dropping the above average time segments.
%
%   See also VIEWCELL2b

% Bálint Király, 08/02/24

% Load CellBase
loadcb
Lratio = getvalue('Lr_PC');
ID = getvalue('ID_PC');
Area = getvalue('Area1');
Mousenum = getvalue('RatID');
Tetrodenum = getvalue('Tetrode');
isArea = strcmp(Area,'HDB');

% Find valid cellids
ptinx = ID > 20 & Lratio < 0.15 & isArea;
I = find(ptinx);
IDS = [CELLIDLIST(I) 'HDB30_181002a_2.1' 'HDB30_181002a_2.2'];

% Find putative tagged cells
cellids = [];
cellids_all = select_pv_cells(cellids); % select PV+ cells from CELLIDLIST of chosen CellBase
cellids = setdiff(cellids_all, {'HDB23_180221a_5.2' 'HDB17_170810a_4.1'}); % remove cells with poor tagging
taggedIDS = [cellids 'HDB30_181002a_2.1' 'HDB30_181002a_2.2']; % add HDB30 neurons
skipinx=[];

default = {...
    'window',               [-3 3];... % time window around the 1st trigger event: calculated  // % !!! do not set it too narrow: longer increased or decsreased block can pull away the mean value of Z-score !!!
    'pwindow'               [-0.5,2];...% time window around the 1st trigger: plotted
    'cwindow',              [0.0,1.7;0.0,1.7];... % time window around the 1st trigger: considerd for the clustering
    'dt',                   0.005;... % time resolution (s)
    'sigma',                0.01;...  % smoothing kernel
    'Trigger1Name',         'StimulusOn';... % 1st Trigger --> time(0:Trigger2)  //e.g. cue
    'Trigger2Name'          'DeliverAllFeedback';... % 2nd Trigger --> time(Trigger2:end) //e.g. feedback
    'Partitions',           '#TrialType';...  % trial partitions
    'ClusterPartitions'     [1,2];... % partition # used for clustering e.g. [1,2]; the first one wil be used for sorting
    'PlotPartitions'        [1,2];... % partition # to be plotted
    'baselinelength'        1;... % length of the baseline period (s) from the begining of the 'window' //only for auROC
    'ROCWindow'             10;... % length of the moving window (in dt points), used for the auROC histogram comparisson
    'avgtimediff'           1.4;... % ~avg time difference between triggers
    'integrallength'        0.6;... % length of the window (from the first trigger event) used for sorting clusters and cells
    'maxtimediff'           2;...  % maximum time difference between triggers
    'ClusterNum'            5;... % required number of clusters
    'PCA_dim'               [3,3];...
    'clustering_method'     'complete';...
    'clustering_metric'     'euclidean';...
    'LastEvents',           '';...
    'BurstPSTH'             'off';...
    'showtagged'            '*';... % '*':mark tagged cells with * / 'IDS': write tagged IDS  / 0: don't mark tagged neurons / 'showall' write all IDS
    'norm_method'           'Z-score';... % 'auROC' / 'Z-score'
    'usedoubletrigger'      1;... % 1:yes; 0:no
    'savefig'               1;... %1:yes, 0:no
    'savematrix'            1;... %1:yes, 0:no
    'avgplotcolors'         ['g','r','y','b'];...
    'savepath'              'F:\HDB_PV\clustering',...
    };

%--------------------------------------------------------------------------
% Calculate the Z-score and AUROC psth of the cells
%--------------------------------------------------------------------------

[g,error] = parse_args(default,varargin{:});
%[g,error] = parse_args(default);
mkdir(g.savepath)

% Defining time windows
g.baselinelength = g.baselinelength/g.dt;
margin = g.sigma*3;
time = g.window(1):g.dt:g.window(2);
time_m = g.window(1)-margin:g.dt:g.window(2)+margin; %plotted time window + margin for smoothing
time_plot = g.pwindow(1):g.dt:g.pwindow(2); % plotted time window
time_calc = g.window(1)-g.maxtimediff-margin:g.dt:g.window(2)+g.maxtimediff+margin; % extended time window for calculation
for partnum = 1:length(g.ClusterPartitions)
    ctime = g.cwindow(partnum,1):g.dt:g.cwindow(partnum,2); % time window considered for clustering
    [~,c1] = (min(abs(time-ctime(1))));
    [~,c2] = (min(abs(time-ctime(end))));
    cinx{partnum,:} = c1:c2;
end

i = 0;
for cellinx = 1:length(IDS)
    
    cellid = IDS(cellinx);
    
    % load event files
    TE = loadcb(cellid,'TrialEvents');
    SP = loadcb(cellid,'EVENTSPIKES');
    
    % Partition trials, and skip cells without at least 1 trial from all required partitions
    try
        [COMPTRIALS2, ~] = partition_trials(TE,g.Partitions);
        COMPTRIALS = COMPTRIALS2(g.ClusterPartitions);
        if min(cellfun('length',COMPTRIALS))<5
            sprintf('%s skipped because of a missing partiton',cellid{1})
            skipinx = [skipinx cellinx];
            continue
        end
    catch
        sprintf('%s skipped because of a missing partiton',cellid{1})
        skipinx = [skipinx cellinx];
        continue
    end
    
    try
        COMPTRIALS = COMPTRIALS2(g.PlotPartitions);
    catch
        sprintf('%s missing partition',cellid{1})
    end
    
    
    % find valid trials in which Trigger2 happened
    trigev1 = TE.(g.Trigger1Name);
    trigev2 = TE.(g.Trigger2Name);
    if ~iscell(trigev2)
        valid_trials = find(~isnan(trigev1) & ~isnan(trigev2));
    else
        valid_trials = find(cellfun(@(s)~isempty(s),trigev2));
    end
    
    % find Trigger2 events
    trigger_pos = findcellstr(SP.events(:,1),g.Trigger2Name);
    
    % Handle TriggerName mismatch
    if trigger_pos == 0
        error('Trigger name not found');
    else
        TriggerEvent = SP.events{trigger_pos,2};
    end
    if ~isfield(TE,TriggerEvent)
        error('TriggerEvent mismatch: supply correct Events structure')
    end
    
    % Spike times
    alltrials = 1:size(SP.event_stimes{1},2);
    stimes  = SP.event_stimes{trigger_pos}(alltrials);
    if strcmp(g.BurstPSTH,'on')
        stimes = detect_bursts(stimes);
    end
    
    % Event windows
    if ~iscellstr(g.LastEvents) && (strcmpi(g.LastEvents,'none') || isempty(g.LastEvents))
        ev_windows = SP.event_windows{trigger_pos};
    else
        ev_windows = get_last_evtime(TE,TriggerEvent,g.LastEvents);
    end
    
    % Calculate binraster
    NUMtrials = length(TE.(g.Trigger1Name));
    if iscell(stimes{1})   % deal with lick-aligned raster
        stimes2 = [stimes{1:end}];
        binraster0 = stimes2binraster(stimes2,time_calc,g.dt);
        binraster = nan(NUMtrials,size(binraster0,2));
        %     binraster2 = nan(NUMtrials,size(binraster0,2));
        for k = 1:NUMtrials   % calculate sum of rows for each trial, which will be used for the PSTH
            sind = sum(cellfun(@length,stimes(1:k-1))) + 1;
            eind = sind + length(stimes{k}) - 1;
            %         disp([sind eind])
            binraster(k,:) = mean(binraster0(sind:eind,:),1);
            %         binraster2(k,:) = sum(stimes2binraster(stimes{k},time,g.dt),1);
        end
    else
        binraster = stimes2binraster(stimes,time_calc,g.dt);
    end
    
    % For variable windows, change padding to NaN to ensure correct averaging - BH
    if ~isempty(g.LastEvents)
        for iT = 1:NUMtrials    % loop through trials
            inx = time_calc > ev_windows(iT,2);
            binraster(iT,inx) = NaN;
        end
    end
    
    % Calculate double triggered rasters
    if g.usedoubletrigger == 1
        trigtime = find(time_calc==0);
        EventTimes = trialevents2relativetime(TE,TriggerEvent,g.Trigger1Name); % time of trgi1 relative to trig2
        Braster = NaN(size(binraster,1),length(time_m)); % inicialize double triggered raster
        trigtimediff = floor((TE.(TriggerEvent)-TE.(g.Trigger1Name))/g.dt); % timedifference between triggers
        trigtimediff(trigtimediff>g.avgtimediff/g.dt) = g.avgtimediff/g.dt; % cut trials with longer time diff than avg
        for j = valid_trials
            if abs(EventTimes(j)) < g.maxtimediff  % check if the time diff is too long because of a recordign mistake and skip problematic trials
                B1 = NaN(1,floor(g.avgtimediff/g.dt)); % inicialize raster between triggers
                % calculate raster between triggers; in case of shorter
                % then avg time diff between triggers bins reamin filled with nans
                B1(1:trigtimediff(j)) = binraster(j,trigtime+floor(EventTimes(j)/g.dt):trigtime+floor(EventTimes(j)/g.dt)+trigtimediff(j)-1);
                startinx = trigtime+floor((EventTimes(j)/g.dt))+floor(time_m(1)/g.dt):trigtime+floor(EventTimes(j)/g.dt)-1; %time before trig1
                %fill rester [before trig1; between triggers; after trig2 ]
                Braster(j,:) = [binraster(j,startinx),B1,binraster(j,trigtime:trigtime+ceil(time_m(end)/g.dt-g.avgtimediff/g.dt))];
            else
                sprintf('Omitted trial in %s (time difference between the trigger events bigger than maxtimediff) \n', cellid{1})
            end
        end
        
    end
    
    % Calculate PSTH
    if g.usedoubletrigger == 1
        [psth, spsth, ~] = binraster2psth(Braster,g.dt,g.sigma,COMPTRIALS,valid_trials);
        psth = psth(:,find(abs((time_m-g.window(1)))<g.dt*0.99):find(abs((time_m-g.window(2)))<g.dt*0.99));
        spsth = spsth(:,find(abs((time_m-g.window(1)))<g.dt*0.99):find(abs((time_m-g.window(2)))<g.dt*0.99));
    else
        g.Trigger1Name = g.Trigger2Name;
        [psth, spsth, ~] = binraster2psth(binraster,g.dt,g.sigma,COMPTRIALS,valid_trials);
        psth = psth(:,floor(g.maxtimediff/g.dt):end-floor(g.maxtimediff/g.dt));
        spsth = spsth(:,floor(g.maxtimediff/g.dt):end-floor(g.maxtimediff/g.dt));
        g.avgtimediff = 0;
    end
    
    %--------------------------------------------------------------------------
    % Normalize PETHs
    %--------------------------------------------------------------------------
    test = sum(isnan(zscore(squeeze(spsth(:,:)),0,2))')>1;
    if sum(test(g.ClusterPartitions))>0
        sprintf('%s skipped because there was not enough data from a partiton',cellid{1})
        skipinx = [skipinx cellinx];
        continue
    end
    i = i + 1;
    
    trigtime = find(time_m==0);
    for k = 1:length(g.PlotPartitions)
        if size(spsth,1)<length(g.PlotPartitions) & ismember(k,setdiff(g.PlotPartitions,g.ClusterPartitions))
            continue
        end
        % Z-score map
        Zpsth(i,k,:) = zscore(squeeze(spsth(k,:)),0,2); % Z-score spsth
        integralrespons(i,k) = sum(Zpsth(i,k,trigtime:trigtime+floor(g.integrallength/g.dt))); %integrated respons from trigger1 for sorting
        
        % auROC map
        maximum = floor(max(psth(k,:)))+1; %maximal length of the histogram
        b = histcounts(psth(k,1:g.baselinelength),'BinWidth',1,'BinLimits',[0 , maximum]); % computing baselinehistorgram
        baselinehist = b/g.baselinelength; % baseline histogram normalization
        z = 0;
        for q = g.ROCWindow:g.ROCWindow:length(psth)
            z = z + 1;
            a = histcounts(psth(k,q-(g.ROCWindow-1):q),'BinWidth',1,'BinLimits',[0 , maximum]); %computing moving window histogram
            ahist = a / g.ROCWindow;         %moving window histogram normalization
            for qq = 1:1:maximum                  % moving the criteria
                p_a(qq) = sum(ahist(qq:end));     % probability that the activity during the window is bigger than the criteria
                p_b(qq) = sum(baselinehist(qq:end)); % probability that the activity during the baseline is bigger than the criteria
                
            end
            auROC(i,k,z) = trapz(p_b,p_a) *- 1; % the area under the p_a-p_b curve (quantifies the degree of overlap between the two spike count distirbuitons)
            clear p_b
            clear p_a
        end
        
    end
    
end


%--------------------------------------------------------------------------
% Clustering based on the first three principal components of the of the
% clusterpartions of the normalized PETHs
%--------------------------------------------------------------------------


% Handle skipped cellids
cellnum = 1:size(Zpsth,1);
IDS_orig = IDS;
IDS(skipinx) = [];
taggednum = find(ismember(IDS,taggedIDS)==1);


% Calcualte pricinpal components based on  the selected normalization
% method
PCA2 = [];
VAR = [];
for f = 1:length(g.ClusterPartitions)
    if g.norm_method == 'Z-score'
        [~,PCA1,~,~,VAR1] = pca(((squeeze(Zpsth(:,f,cinx{f}))))); % calculate principal copmponents of the Z-scored psth in the clustering time window
    elseif g.norm_method == 'auROC'
        PCA1 = pca(((squeeze(auROC(:,f,floor(cinx{f}(1)/g.ROCWindow):floor(cinx{f}(end)/g.ROCWindow))))')); % calculate principal copmponents of the auROC in the clustering time window
    end
    PCA2 = [PCA2,PCA1(:,1:g.PCA_dim(f))]; % take the 1st x=PCA_dim principal components
    VAR = [VAR;VAR1(1:g.PCA_dim(f))]; % take the 1st x=PCA_dim principal components
end
figure
Clusters = kmeans(PCA2,g.ClusterNum,'Replicates',100);
SORT_ = zeros(1,g.ClusterNum);
for j = 1:g.ClusterNum
    SORT_(j) = mean(PCA2(Clusters==j,1));
end
[~,rule] = sort(SORT_);
TEMP = zeros(length(Clusters),1);
for j = 1:g.ClusterNum
    TEMP(Clusters==rule(j))=j;
end
scatter3(PCA2(:,1),PCA2(:,2),PCA2(:,3),50,TEMP*100,'.'); % 3D PCA cluster plot

xlabel('PCA_1');
ylabel('PCA_2');
zlabel('PCA_3');

mkdir(g.savepath)
if g.savefig == 1
    saveas(gca,[g.savepath '3DPCA.fig'])
end

[~,inx] = sortrows([TEMP,PCA2(:,3)]); %inx: order of the sorted cells (original index)
clusterborders = find(diff(Clusters(inx))~=0);

%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------


clustered auROC maps

figure;
suptitle(sprintf('auROC maps, clusters based on %s',g.norm_method));
for k = 1:size(auROC,2) %rows -- corresponding to different partitions
    subplot(size(auROC,2),8,[(k-1)*8+1:(k-1)*8+6])
    imagesc(time,cellnum,squeeze(auROC(inx,k,:)));
    colormapdefiner(1,0,0.5,100,size(auROC,2),8,[(k-1)*8+1:(k-1)*8+6],g.Trigger1Name,g.Trigger2Name,-1.5);
    hold on
    additionalplots(time,time_plot,clusterborders,cellnum(end),g,inx,taggednum,IDS);
    hold off
    subplot(size(auROC,2),8,(k-1)*8+7)
    imagesc(PCA2(inx,:));
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
end
if g.savefig == 1
    saveas(gca,[g.savepath 'auROCmap.fig'])
end

% clustered Z-score maps

figure;
suptitle(sprintf('Z-score maps, clusters based on %s',g.norm_method));
for k = 1:size(auROC,2) % rows -- corresponding to different partitions
    subplot(size(auROC,2),8,[(k-1)*8+1:(k-1)*8+6])
    N = imagesc(time,cellnum,squeeze(Zpsth(inx,k,:)));
    c_limit = quantile(N.CData(:),[0.9995,0.0005]);
    colormapdefiner(c_limit(1),c_limit(2),0,100,length(g.PlotPartitions),8,[(k-1)*8+1:(k-1)*8+6],g.Trigger1Name,g.Trigger2Name,-1)
    hold on
    additionalplots(time,time_plot,clusterborders,cellnum(end),g,inx,taggednum,IDS);
    hold off
    subplot(size(auROC,2),8,(k-1)*8+7)
    imagesc(PCA2(inx,:));
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
end
if g.savefig == 1
    saveas(gca,[g.savepath 'Zscoremap.fig'])
end
sortedZpsth = Zpsth(inx,:,:);
sortedCellIDS = IDS(inx);
%--------------------------------
%mean Zscores
figure;
for j = 1:g.ClusterNum
    subplot(1,5,j)
    hold on
    for k = 1:length(g.PlotPartitions)
        borders = [0;clusterborders;cellnum(end)];
        avgZscore(j,k,:) = mean(sortedZpsth(borders(j)+1:borders(j+1),k,:),1);
        stdZscore(j,k,:) = std(sortedZpsth(borders(j)+1:borders(j+1),k,:),1)./sqrt(borders(j+1)-borders(j));
        set(gcf, 'Renderer', 'painters');
        errorshade(time, squeeze(avgZscore(j,k,:)),squeeze(stdZscore(j,k,:)), 'LineColor',g.avgplotcolors(k), 'ShadeColor',g. avgplotcolors(k), 'LineWidth', 3);
        hold on
        min(P.YTick)
        line([0,0],[min(P.YTick),max(P.YTick)],'Color','magenta','LineWidth',1.5);
        line([g.avgtimediff,g.avgtimediff],[min(P.YTick),max(P.YTick)],'Color','cyan','LineWidth',1.5);
        xlim([time_plot(1),time_plot(end)]);
        ylim([-2,4])
        setmyplot_balazs
        if g.savefig == 1
            saveas(gca,[g.savepath 'meanZscore_' int2str(j) '.fig'])
        end
    end
end

% Bar plot depicting explained variance by the PCA components
figure
hold on
B_VAR = zeros(length(g.ClusterPartitions),max(g.PCA_dim));
for i_cp = 1:length(g.ClusterPartitions)
    B_VAR(i_cp,1:g.PCA_dim(i_cp)) = VAR((1:g.PCA_dim(i_cp))+sum(g.PCA_dim(1:i_cp-1)));
end
bar(1:length(g.ClusterPartitions),B_VAR,'stacked')
set(gca,'xtick',1:length(g.ClusterPartitions));
set(gca,'xticklabel',{'Hit_4_H_z','FalseAlarm'});
ylabel('Explained variance (%)')
saveas(gca,[g.savepath 'var.fig'])

% -------------------------------
% sav data for futher processing
sortedClusters = ones(borders(2),1);
for cb = 2:length(borders)-1
    sortedClusters = [sortedClusters;ones(borders(cb+1)-borders(cb),1)*cb];
end


if g.savematrix == 1
    save([g.savepath 'DATA.mat'],'sortedCellIDS','sortedClusters','sortedZpsth','g','IDS_orig','clusterborders','inx','IDS','skipinx','VAR')
    filename = [g.savepath 'Clusters.xlsx'];
    xlswrite(filename,sortedCellIDS')
    xlswrite(filename,sortedClusters,'B1:B10000')
    %writetable('sortedClusters',filename,'Sheet',1,'Range','B1')
    
end

% -------------------------------------------------------------------------
function colormapdefiner(largest,smallest,indexValue,binnum,s1,s2,f,trig1,trig2,scalcefactor)
% logarithmic scale for the auROC map (red - green)
% Calculate where proportionally indexValue lies between minimum and maximum values
index = abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [logspace(0,scalcefactor,round(binnum*index))',...
    zeros(round(binnum*index),1),...
    zeros(round(binnum*index),1)];

customCMap2 = [zeros(round(binnum*(1-index)),1),...
    logspace(scalcefactor,0,round(binnum*(1-index)))',...
    zeros(round(binnum*(1-index)),1)];
customCMap = [customCMap1;customCMap2];  % Combine colormaps
colormap(subplot(s1,s2,f),customCMap);
caxis ([smallest,largest]);
xlabel(sprintf('time(s) - %s - %s',trig1,trig2));
ylabel('cell#');

% -------------------------------------------------------------------------
function additionalplots(time,time_plot,clusternum,cellnum,g,inx,taggednum,IDS)

line([time(1),time(end)],[clusternum+0.5,clusternum+0.5],'Color','yellow');
line([0,0],[0,cellnum+0.5]','Color','magenta','LineWidth',1.5);
line([g.avgtimediff,g.avgtimediff],[0,cellnum]','Color','cyan','LineWidth',1.5);
[t1,t2] = ismember(squeeze(inx),taggednum);
t1 = find(t1==1);
t3 = IDS(inx(t1));
%t2 = t2(t2~=0);
%[t2,t4] = sort(t2);
if isequal(g.showtagged,'*')
    scatter(ones(1,length(taggednum))*g.pwindow(1)+g.dt,t1,[],[0.5843 0.8157 0.9882],'*');
elseif isequal(g.showtagged,'IDS')
    ylabel('tagged cell IDs');
    set(gca,'ytick',[t1],'yticklabel',t3);
elseif isequal(g.showtagged,'showall')
    ylabel('cell IDs');
    set(gca,'ytick',[1:cellnum],'yticklabel',IDS(inx));
end
xlim([time_plot(1),time_plot(end)])