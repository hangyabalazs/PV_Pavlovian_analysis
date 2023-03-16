function PV_pavlovian_analysis_main(choosecells, resdir, preprocess, behavior, anatomy, recording, optogenetics)
%PV_PAVLOVIAN_ANALYSIS_MAIN main data analysis funtion of PV-pavlovian project. 
%   PV_PAVLOVIAN_ANALYSIS_MAIN(CHOOSECELLS, RESDIR, PREPROCESS, BEHAVIOR,
%   ANATOMY, RECORDING, OPTOGENETICS) main data analysis function of in
%   vivo recordings of PV+ basal forebrain neurons. The function calls
%   preprocessing functions, analyzing behavior during a probabilistic
%   pavlovian conditioning task, recording, optogenetics manipulation and
%   anatomy. These features are controlled by PREPROCESS(1x5 logical)
%   BEHAVIOR(1x2 logical), ANATOMY (1x1 logical), RECORDING(1x9 logical)
%   and OPTOGENETICS(1x2 logical) inputs, respectively.
%
%   PREPROCESS(1x5 logical) controls response sorter, autocorrelograms,
%   croscorrelograms, tagging statistics and PETH for individual neurons.
%   BEHAVIOR(1x2 logical) controls example and average PETH for behavior.
%   ANATOMY (1x1 logical) adds cell location data. RECORDING(1x9 logical)
%   controls cluster quality statistics, average PETH, latency of
%   punishment response, PETH partitioned by median splint of trials, PETH
%   of response based on previous trial outcome and PETH showing burst and
%   single spikes separately. OPTOGENETICS(1x2 logical) controls statistics
%   for optogenetic inhibition.

% Balazs Hangya, Panna Hegedus, 12/03/2023
% Institut of Experimental Medicine, Budapest
% hangya.balazs@koki.hu


% Input argument check
narginchk(0,7)
if nargin < 1
    choosecells = 1;
end

if choosecells % If you want to perform analysis of in vivo recordings, you need to choose cells
    choosecb('PV_pavlovian')   % choose CellBase
    
    % Getting cellids for VP neurons
    loadcb % load chosen CellBase
    cellids_all = select_pv_cells(cellids); % select PV+ cells from CELLIDLIST of chosen CellBase
    cellids = setdiff(cellids_all, {'HDB23_180221a_5.2' 'HDB17_170810a_4.1'}); % remove cells with poor tagging
    
    cellids = [cellids 'HDB30_181002a_2.1' 'HDB30_181002a_2.2']; % add HDB30 neurons
    NumCells = length(cellids); % use the subset with good cluster quality
end

if nargin < 2 % Create direcotry where results are saved
    resdir = fullfile(getpref('cellbase','datapath'), '_data_analysis');
    if ~isfolder(resdir) % make results directory if not exist
        mkdir(resdir)
    end
end

% Control which analysis modules to perform
% Preprocess
if nargin < 3
    preprocess = [1 1 1 1 1];
end

response_analysis = preprocess(1);   % Response sorter
perform_acg = preprocess(2);   % Autocorrelograms
perform_ccg = preprocess(3);   % Cross-correlograms
tagging_stat = preprocess(4);  % Tagging statistics
perform_quickanalysis = preprocess(5);   % PETH for individual neurons

% Behavior analysis
if nargin < 4
    behavior = [1 1];
end

lick_example = behavior(1);   % example lick PETH and raster
lick_average = behavior(2);   % average lick PETH and raster

% Anatomical location
if nargin < 5
    anatomy = [1];
end

loc = anatomy(1);

% In vivo recording analysis
if nargin < 6
    recording = [1 1 1 1 1 1 1 1];
end

cluster = recording(1);             % Cluster quality statistics
avg_PSTH = recording(2);            % Average PETH
expectation_PSTH = recording(3);    % Average PETH partitioned by trial type
latency = recording(4);             % Latency of punishment response
preprocess_latency = recording(5);  % Data preprocessing
msplit = recording(6);              % PETH of partitioned by median splint of trials
ms_preproc = recording(7);          % Data preprocessing
preproc_prevtrials = recording(8);   % Response based on previous trial outcome
burst_spikes = recording(9);        % PETH showing burst and single spikes separately

% Optogenetic experiment analysis
if nargin < 6
    optogenetics = [1 1];
end

optoinh = optogenetics(1);      % Optogenetic inhibition
singleplot = optogenetics(2);   % PETH for individual sessions

%%% DATA PREPROCESSING %%%
if response_analysis 
    
    % Grouping neurons based on their response during behavior - properties added to TheMatrix
    response_resdir = fullfile(resdir,'responsesorter');   % results directory
    
    pv_responsesorter(cellids,1,response_resdir,'cue');
    pv_responsesorter(cellids,1,response_resdir,'rew');
    pv_responsesorter(cellids,1,response_resdir,'pun');
    pv_responsesorter(cellids,1,response_resdir,'om');
end

% Auto-correlation
acg_resdir = fullfile(resdir,'ACG');   % results directory
if perform_acg
    pv_acg(cellids,acg_resdir,true);
end

% Cross-correlation
ccg_resdir = fullfile(resdir,'CCG');   % results directory
if perform_ccg
    pv_ccg(cellids,ccg_resdir,true);
end

% Tagging efficiency
impdir_LS = fullfile(resdir,'taggedprop');
if tagging_stat
    taggedprop(cellids, true);
end

% Individual neuronal responses
resdir_quickanalysis = fullfile(resdir, 'quickanalysis');
if perform_quickanalysis
    quickanalysis_pavlovian_PV(cellids, resdir_quickanalysis)
end

%%% BEHAVIOR ANALYSIS %%%
if lick_example % Plot behavior in individual sessions
    lick_raster_resdir = fullfile(resdir, 'lick_rasters'); % Results directory
    
    if ~isfolder(lick_raster_resdir) % Make folder if not exist
        mkdir(lick_raster_resdir)
    end

    % Example lick PETH + raster
    for i =1:NumCells
        current_cell = cellids{i};
        [animal, session, ~, ~] = cellid2tags(current_cell);
        
        figure;  
        viewlick({animal session},'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav',...
            'ShowEvents',{{'StimulusOn' 'StimulusOff'}}, 'sigma', 0.08,...
            'Partitions','#TrialType','window',[-5 5]   )% lick rasters aligned to stimulus onset
        fnm = fullfile(lick_raster_resdir, [animal '_' session '.jpg']);
        fnm2 = fullfile(lick_raster_resdir, [animal '_' session '.fig']);
        set(gcf, 'renderer', 'painters');
        fnm3 = fullfile(lick_raster_resdir, [animal '_' session '.eps']);
        saveas(gcf, fnm)
        saveas(gcf, fnm2)
        saveas(gcf, fnm3)
    end
end

if lick_average % Plot average behavior PETH
    resdir_lick_psth = fullfile(resdir,'lick_psth');   % Results directory
    resdir_lick_psth_allcond = fullfile(resdir,'lick_psth_allcond');   % Results directory
    lick_psth_summary(cellids,resdir_lick_psth,true)  % Average lick PETH
    lick_psth_summary_allcond(cellids,resdir_lick_psth_allcond,true) % Average lick PETH parititoned to all possible outcomes
end

%%% ANATOMICAL LOCATION %%%
if loc % Add tetrode position log (cell location) to cellbase
    insertdata('H:\_personal\auditory_pavlovian_PV_cellbase\_data_analysis\anatomy\cell_location\PV_neuron_location_AP.xlsx','type','prop','name','AP')
    insertdata('H:\_personal\auditory_pavlovian_PV_cellbase\_data_analysis\anatomy\cell_location\PV_neuron_location_DV.xlsx','type','prop','name','DV')
    insertdata('H:\_personal\auditory_pavlovian_PV_cellbase\_data_analysis\anatomy\cell_location\PV_neuron_location_L.xlsx','type','prop','name','L')
    insertdata('H:\_personal\auditory_pavlovian_PV_cellbase\_data_analysis\anatomy\cell_location\PV_neuron_location_Area.xlsx','type','prop','name','Area1')
end

%%% RECORDING ANALYSIS %%%
% Cluster quality statistics
if cluster
    resdir_C = fullfile(resdir, 'clusters');
    cluster_quality_stats(cellids, resdir_C, 'PV'); % Cluster quality statistics
    
    sc = selectcell('"ID_PC">20&"Lr_PC"<0.15');  % cells with criterion quality
    ntcells = setdiff(sc,cellids);  % untagged population
    cluster_quality_stats(ntcells, resdir_C, 'untagged')  % untagged
    
    % H index distribution
    tagging(CELLIDLIST,true)
    H_indices = nan(1,length(CELLIDLIST));
    for n=1:length(CELLIDLIST)
        cellid = CELLIDLIST{n};
        [H_indices(n) ~] = nbisstim(cellid);
        disp(H_indices(n))
    end
    
    cells = CELLIDLIST; % Set H-indices for all neurons
    for n=1:length(cells)
        st = setvalue(cells{n},'Hindex',H_indices(n));   
    end
    
    H = Hindex_distribution(ntcells,cellids);
    fnm = fullfile(resdir_C, 'tagging_histogram.fig');
    saveas(H,fnm);
    set(H, 'renderer', 'painters')
    fnm = fullfile(resdir_C, 'tagging_histogram.eps');
    saveas(H,fnm);
    close(H)
    
    % Latency and jitter CDF for light activation
    PV_R_L_J(cellids, resdir_C, 'tag', true)
    
    % Color coded PSTH for light stimulation
    PV_lightpsth_colored(cellids,fullfile(resdir,'taggingsummary'), resdir_C);
end

% Example tagged cluster
    [selcluster, neighborclusters, selLS] = example_cluster('HDB36_190512a_3.2');
    fnm = fullfile(resdir_C, 'example_cluster.mat');
    save(fnm,'selcluster','neighborclusters','selLS');

% Average PETH of cue, reward and punishment response
if avg_PSTH
    population = getcellresp(cellids); % group neurons based on their response
    
    resdir_avg_psth = fullfile(resdir, 'avg_PSTH_nonpart');
    avg_psth_PV2(cellids,'cueresponse',resdir_avg_psth); 
    avg_psth_PV2(cellids,'rewardresponse',resdir_avg_psth);
    avg_psth_PV2(cellids,'punishresponse',resdir_avg_psth);
    avg_psth_PV2(cellids,'omissionresponse',resdir_avg_psth);
    
    % Neurons partitioned based on activation and inhibition to cue and reinforcement
    resdir_avg_psth2 = fullfile(resdir, 'avg_PSTH_nonpart2');
    avg_psth_PV(population.cue.excitation,population.cue.inhibition,'cue',resdir_avg_psth2)
    avg_psth_PV(population.reward.excitation,population.reward.inhibition,'reward',resdir_avg_psth2)
    avg_psth_PV(population.punishment.excitation,population.punishment.inhibition,'punishment',resdir_avg_psth2)
    
    % Color coded matrix of cell responses
    [~, G] = response_matrix(cellids);
    saveas(G,fullfile(resdir_avg_psth,'imagesc_cellreponses.fig'))
    saveas(G,fullfile(resdir_avg_psth,'imagesc_cellreponses.eps'))
end

% avg PETH of response to expected and unexpected punishment
if expectation_PSTH
    
    % Only neurons activated by punishment
    resdir_avg_psth_expectation = fullfile(resdir, 'avg_PSTH_expectation_punishcells');
    avg_psth_expectation_PV(population.punishment.excitation, 'cue', 'none', resdir_avg_psth_expectation, 'real');
    avg_psth_expectation_PV(population.punishment.excitation, 'reward', 'none', resdir_avg_psth_expectation, 'real');
    avg_psth_expectation_PV(population.punishment.excitation, 'punish', 'none', resdir_avg_psth_expectation, 'real');
    
    % All PVBFNs
    resdir_avg_psth_expectation = fullfile(resdir, 'avg_PSTH_expectation_allcellid');
    avg_psth_expectation_PV(cellids, 'cue', 'none', resdir_avg_psth_expectation, 'real');
    avg_psth_expectation_PV(cellids, 'reward', 'none', resdir_avg_psth_expectation, 'real');
    avg_psth_expectation_PV(cellids, 'punish', 'none', resdir_avg_psth_expectation, 'real');
end

% Latency of response to expected and unexpected punishment
if latency
    responsesorter_expectation = fullfile(resdir, 'responsesorter_expectation');
    if preprocess_latency
        pv_responsesorter(cellids,1,responsesorter_expectation,'pun_expectation');
    end
    latencies = latency_dist_PV(population.punishment.excitation,responsesorter_expectation, 'punishment_expectation');  % latency properties of cue, reward and punishment response
    
    % Boxplots
    maxvalue1 = latencies.punishment_expectation.maxvalue{1};
    maxvalue2 = latencies.punishment_expectation.maxvalue{2};
    
    latency_peak1 = latencies.punishment_expectation.excited.latency_peak{1};
    latency_peak2 = latencies.punishment_expectation.excited.latency_peak{2};
    
    label1 = latencies.punishment_expectation.tags{1};
    label2 = latencies.punishment_expectation.tags{2};
    boxstat(maxvalue1, maxvalue2, label1, label2, 0.005, 'paired');
    fnm = fullfile(responsesorter_expectation, 'boxstat_maxvalue_punishment_expectation.jpg');
    fnm2 = fullfile(responsesorter_expectation, 'boxstat_maxvalue_punishment_expectation.eps');
    fnm3 = fullfile(responsesorter_expectation, 'boxstat_maxvalue_punishment_expectation.fig');
    set(gcf, 'renderer', 'painters')
    saveas(gcf, fnm)
    saveas(gcf, fnm2)
    saveas(gcf, fnm2)
    close(gcf)
    
    boxstat(latency_peak1, latency_peak2, label1, label2, 0.005, 'paired');
    fnm = fullfile(responsesorter_expectation, 'boxstat_latency_punishment_expectation.jpg');
    fnm2 = fullfile(responsesorter_expectation, 'boxstat_latency_punishment_expectation.eps');
    fnm3 = fullfile(responsesorter_expectation, 'boxstat_latency_punishment_expectation.fig');
    set(gcf, 'renderer', 'painters')
    saveas(gcf, fnm)
    saveas(gcf, fnm2)
    saveas(gcf, fnm2)
    close(gcf)
end

% Median split
if msplit
    if ms_preproc
        add_TE_event(cellids, 'MedialSplit'); % add TE and TrialEvents a MedianSplit variable - make if condition!
        quickanalysis_pavlovian_PV(cellids, resdir_quickanalysis)% individual PETHs
    end
    
    resdir_avg_psth_mediansplit = fullfile(resdir, 'avg_mediansplit'); % Results directory
    avg_psth_expectation_PV(cellids, 'mediansplit', 'none', resdir_avg_psth_mediansplit, 'real'); % avg PETH
end

% Plot burst and single spikes separately
if burst_spikes
    resdir_avg_psth_burst = fullfile(resdir, 'burst_spikes'); % Results directory
    burst_detection(cellids, 'burst', resdir_avg_psth_burst) % Detecting bursts
    burst_detection(cellids, 'single', resdir_avg_psth_burst) % Detecting single spikes
    avg_psth_PV(cellids, cellids, 'punishresponse',resdir_avg_psth_burst, 'burst'); % Average PETH
    acg_matrix = load(fullfile(resdir, 'ACG', 'ACG_matrices.mat'));
   
    % Plot refractory period against burst index
    figure;
    scatter(acg_matrix.Refractory, acg_matrix.BurstIndex, 'o')
    xlabel('Refractory period')
    ylabel('Burst Index')
    
    % Plot refractory period and burst index histograms
    figure;
    histogram(acg_matrix.Refractory, 20)
    title('Refractory period')
    figure;
    histogram(acg_matrix.BurstIndex, 20)
    title('Burst Index')
    
    thr = 0.3; % treshold for selecting bursting neurons
    resdir_sorting_burst = fullfile(resdir, 'ACG', 'bursting');
    resdir_sorting_nonburst = fullfile(resdir, 'ACG', 'nonbursting');
    acg_sorting(cellids, acg_resdir, resdir_sorting_burst, resdir_sorting_nonburst, thr) % sorting bursting and non-bursting neurons
    
    bursting_cells = cellids((acg_matrix.BurstIndex>thr)&acg_matrix.Refractory<=2); % Bursting neurons
    nonbursting_cells = setdiff(cellids, bursting_cells); % Non-bursting neurons
    
    bursting = getcellresp(bursting_cells); % Cell responses
    nonbursting = getcellresp(nonbursting_cells);
    
    % Avg PETH
    resdir_avg_psth_VP1 = fullfile(resdir,'avg_PSTH_bursting');
    avg_psth_VP(bursting.cue.excitation,bursting.cue.inhibition,'cueresponse',resdir_avg_psth_VP1);
    avg_psth_VP(bursting.reward.excitation,bursting.reward.inhibition,'rewardresponse',resdir_avg_psth_VP1);
    avg_psth_VP(bursting.punishment.excitation,bursting.punishment.inhibition,'punishresponse',resdir_avg_psth_VP1);
    
    resdir_avg_psth_VP2 = fullfile(resdir, 'avg_PSTH_nonbursting');
    avg_psth_VP(nonbursting.cue.excitation,nonbursting.cue.inhibition,'cueresponse',resdir_avg_psth_VP2);
    avg_psth_VP(nonbursting.reward.excitation,nonbursting.reward.inhibition,'rewardresponse',resdir_avg_psth_VP2);
    avg_psth_VP(nonbursting.punishment.excitation,nonbursting.punishment.inhibition,'punishresponse',resdir_avg_psth_VP2);
    
    %      Compare responsive neuronal proportions - chi square-test
    [tbl1, chi2stat1, pval1] = chi2_vp(length([bursting.cue.excitation' bursting.cue.inhibition']),length(bursting_cells),length([nonbursting.cue.excitation' nonbursting.cue.inhibition']),length(nonbursting_cells));
    [tbl2, chi2stat2, pval2] = chi2_vp(length([bursting.reward.excitation' bursting.reward.inhibition']),length(bursting_cells),length([nonbursting.reward.excitation' nonbursting.reward.inhibition']),length(nonbursting_cells));
    [tbl3, chi2stat3, pval3] = chi2_vp(length([bursting.punishment.excitation' bursting.punishment.inhibition']),length(bursting_cells),length([nonbursting.punishment.excitation' nonbursting.punishment.inhibition']),length(nonbursting_cells));
    
    % Maxvalue stats
    resdir_maxval_comp_bursting = fullfile(resdir,'compare_maxval_bursting');
    avg_psth_compare_maxval_VP(bursting.cue.excitation,'cueresponse',nonbursting.cue.excitation,'cueresponse','excitation',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.reward.excitation,'rewardresponse',nonbursting.reward.excitation,'rewardresponse','excitation',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.punishment.excitation,'punishresponse',nonbursting.punishment.excitation, 'punishresponse','excitation',resdir_maxval_comp_bursting);
    
    avg_psth_compare_maxval_VP(bursting.cue.inhibition,'cueresponse',nonbursting.cue.inhibition,'cueresponse','inhibition',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.reward.inhibition,'rewardresponse',nonbursting.reward.inhibition,'rewardresponse','inhibition',resdir_maxval_comp_bursting);
    avg_psth_compare_maxval_VP(bursting.punishment.inhibition,'punishresponse',nonbursting.punishment.inhibition,'punishresponse','inhibition',resdir_maxval_comp_bursting);
    
    % Average ACGs
    segfilter = 'stim_excl_vp'; % exclude recording segments with stimulation
    filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
    resdir_acg_average = fullfile(resdir,'ACG','acg_average');
    if ~isfolder(resdir_acg_average)
        mkdir(resdir_acg_average)
    end
    avg_vpacg(bursting_cells,0.5,'resdir',resdir_acg_average,...
        'segfilter',segfilter,'filterinput',filterinput,...
        'minspikeno',100,'maxspikeno',10000,'issave',true);  % Average ACG for bursting neurons
    avg_vpacg(nonbursting_cells,0.5,'resdir',resdir_acg_average,...
        'segfilter',segfilter,'filterinput',filterinput,...
        'minspikeno',100,'maxspikeno',10000,'issave',true);  % Average ACG for non-bursting neurons
    
    % Spike shape analysis
    if ~ismember('SpikeShape',listtag('prop')) % perform spike shape analysis if properties are not added to CellBase yet
        resdir_spikeshape = fullfile(resdir,'spikeshapeanalysis'); % saveresults here
        vp_spikeshape_analysis(resdir_spikeshape);
    end
    
    % Average spike shape
    [~, b_inx] = intersect(cellids, bursting_cells);
    [~, nb_inx] = intersect(cellids, nonbursting_cells);
    SpikeShapes_cellids = getvalue('SpikeShape_test', cellids);
    spike_waveforms = nan(NumCells,27);
    [valley_cellid, PeakToPostValleyTime, integral_cellid] = deal(nan(1,NumCells));
    
    figure;
    for s = 1:NumCells
        spike_waveforms(s,:) = SpikeShapes_cellids{s,1}{1,1}.Spike;
        valley_cellid(s) = SpikeShapes_cellids{s,1}{1,1}.Valley;
        PeakToPostValleyTime(s) = SpikeShapes_cellids{s,1}{1,1}.PeakToPostValleyTime;
        integral_cellid(s) = SpikeShapes_cellids{s,1}{1,1}.Integral;
        
        spike_waveforms(s,:) = spike_waveforms(s,:)./max(spike_waveforms(s,:));
        plot(spike_waveforms(s,:))
        hold on
    end
    
    % Spike shape scatter
    ValleyToIntegral = valley_cellid./integral_cellid;
    figure;
    scatter(ValleyToIntegral, PeakToPostValleyTime, 'o')
    xlabel('ValleyToIntegral')
    ylabel('PeakToPostValley')
    hold on
    scatter(ValleyToIntegral(b_inx), PeakToPostValleyTime(b_inx), 'ro')
    
        inx1 = find((ValleyToIntegral >= 0.06) & (PeakToPostValleyTime <= 240));
    inx2 = find((ValleyToIntegral < 0.06) & (PeakToPostValleyTime > 240));
    inx3 = find(PeakToPostValleyTime <= 240);
    
    cellids_inx1 = cellids(inx1);
    cellids_inx2 = cellids(inx2);
    cellids_inx3 = cellids(inx3);
    
    % Calculate regression
    polypredcicall(ValleyToIntegral,acg_matrix.BurstIndex,0.95,'robust')
    xlabel('ValleyToIntegral');
    ylabel('Burst index');
    
    % Calculate regression
    polypredcicall(PeakToPostValleyTime,acg_matrix.BurstIndex,0.95,'robust')
    xlabel('PeakToPostValleyTime');
    ylabel('acg_matrix.BurstIndex');
    
    % Average normalized spike shape
    figure;
    errorshade([1:27],nanmean(spike_waveforms(b_inx,:)),nanstd(spike_waveforms(b_inx,:))/sqrt(size(spike_waveforms(b_inx,:),1)),...
        'LineColor',[0 0 1],'ShadeColor',[1 0 0]) % excitation
    hold on
    errorshade([1:27],nanmean(spike_waveforms(nb_inx,:)),nanstd(spike_waveforms(nb_inx,:))/sqrt(size(spike_waveforms(nb_inx,:),1)),...
        'LineColor',[0 0 1],'ShadeColor',[0 0 1]) % excitation
    
    % CCG pairs: tagged- non tagged
    ccg_tagged_vs_nontagged = fullfile(resdir, 'CCG', 'tagged_vs_nontagged');
    pv_nonpv_ccg(cellids,ccg_tagged_vs_nontagged,true);
end


% Response based on previous trials
if preproc_prevtrials
    add_TE_event(cellids, 'PreviousTrial1'); % add TE and TrialEvents a MedianSplit variable
    add_TE_event(cellids, 'PreviousTrial2'); % add TE and TrialEvents a MedianSplit variable
    add_TE_event(cellids, 'AllPreviousTrial'); % add TE and TrialEvents a MedianSplit variable
    quickanalysis_pavlovian_PV(cellids, resdir_quickanalysis); % individual PSTHs
    
    resdir_regression = fullfile(resdir, 'regression2previoustrial'); % Results directory
    regression_prevtrial(cellids, resdir_regression);
    keyboard
end


%%% OPTOGENETIC INHIBITION DURING PAVLOVIAN CONDIITONING %%%
if optoinh
    resdir_optoinh = fullfile(resdir, 'Opto_inhibition', 'Single_plots'); % Results directory
    subresdir_Control = fullfile(resdir_optoinh, 'Control'); % Control subfolder in single plots
    subresdir_ArchT = fullfile(resdir_optoinh, 'ArchT'); % ArchT subfolder in single plots
    root = 'E:\HDB_PV\auditory_pavlovian_PV_cellbase\_data_analysis\_ArchT_pavlovian\Analysis\Data2Use'; % Root directory where data is saved
    
    % Make results directory if not exist
    if ~isfolder(resdir_optoinh)
        mkdir(resdir_optoinh);
    end
    
    if ~isfolder(subresdir_Control)
        mkdir(subresdir_Control);
    end
    
    if ~isfolder(subresdir_ArchT)
        mkdir(subresdir_ArchT);
    end
    
    % List ArchT and Control animals
    ArchTdir = fullfile(root, 'ArchT');
    Controldir = fullfile(root, 'Control');
    
    ArchT_folders = dir(ArchTdir);
    ArchT_animals = {ArchT_folders(3:end).name};
    
    Control_folders = dir(Controldir);
    Control_animals = {Control_folders(3:end).name};
    
    % Plot individual behavior sessions
    if singleplot
        all_animals = [ArchT_animals Control_animals];
        for i = 1:length(all_animals)
            c_animal = all_animals{i};
            if ismember(c_animal, ArchT_animals)
                animal_folders = dir(fullfile(ArchTdir, c_animal));
                fullpth = ArchTdir;
                ControlOrNot = 'ArchT';
            else
                animal_folders = dir(fullfile(Controldir, c_animal));
                fullpth = Controldir;
                ControlOrNot = 'Control';
            end
            
            c_sessions = {animal_folders(3:end).name};
            NumSession = length(c_sessions);
            
            for j = 1:NumSession
                TE = solo2trialevents_pavlovian_optoinh(fullfile(fullpth,c_animal, c_sessions{j},[c_animal '_' c_sessions{j} '.mat']),1);
                
                %Lick raster
                H = figure;
                viewlick_opto(fullpth, {c_animal c_sessions{j}},'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#TrialType','window',[-5 5]);
                maximize_figure(H)
                
                resdir_singleplot = fullfile(horzcat(resdir,'\',ControlOrNot), horzcat(c_animal, horzcat('_',ControlOrNot)));
                if ~isfolder(resdir_singleplot)
                    mkdir(resdir_singleplot)
                end
                
                fnm = fullfile(resdir_singleplot,[c_animal '_' c_sessions{j} '_LR' horzcat('_',ControlOrNot) '.jpg']);% save
                saveas(H,fnm)
                close(H)
            end
        end
    end
    
    % Average lick PSTH plot for different training stages in control and ArchT animals
    resdir_box = fullfile(resdir, 'Opto_inhibition_1200ms_win'); % Results directory
    [SummaryDataOutputArchT, Hit_psth_ArchT, FA_psth_ArchT] = lick_psth_summary_PV([], root, resdir_box, 'ArchT');
    [SummaryDataOutputControl, Hit_psth_Control, FA_psth_Control] = lick_psth_summary_PV([], root, resdir_box, 'Control');
    
    % Difference of reaction time to cue1 and cue2 in control and ArchT animals
    ArchT_RT1 = [SummaryDataOutputArchT.AnimalAverageRT1]; % ArchT and Control reaction times
    ArchT_RT2 = [SummaryDataOutputArchT.AnimalAverageRT2];
    Control_RT1 = [SummaryDataOutputControl.AnimalAverageRT1];
    Control_RT2 = [SummaryDataOutputControl.AnimalAverageRT2];
    diff_RT_ArchT = ArchT_RT2 - ArchT_RT1; % reaction time difference per session
    diff_RT_Ctrl = Control_RT2 - Control_RT1;
    
    % Plot reaction time to cue1 and cue2 (control and ArchT)
    figure;
    bar([1 2 ], [mean(Control_RT1) mean(Control_RT2)]);
    hold on
    plot([1 2], [Control_RT1' Control_RT2']);
    
    figure;
    bar([1 2 ], [median(ArchT_RT1) median(ArchT_RT2)]);
    hold on
    plot([1 2], [ArchT_RT1' ArchT_RT2']);
    
    % Statistical comparison of reaction times
    boxstat(diff_RT_Ctrl, diff_RT_ArchT, 'Control', 'ArchT', 0.05) % Reaction time difference
    boxstat(Control_RT1, Control_RT2, 'ControlT1', 'ControlT2', 0.05, 'paired') % Reactin time in control animals
    boxstat(ArchT_RT1, ArchT_RT2, 'AchTT1', 'ArchTT2', 0.05, 'paired') % Reaciton time in ArchT animals
    
    % Plot lick rates to cue1 and cue2 (control and ArchT)
    Control_LR1 = [SummaryDataOutputControl.AnimalAverageT1]; % Control animals' lick rate for cue1
    Control_LR2 = [SummaryDataOutputControl.AnimalAverageT2]; % Control animals' lick rate for cue2
    ArchT_LR1 = [SummaryDataOutputArchT.AnimalAverageT1]; % ArchT animals' lick rate for cue1
    ArchT_LR2 = [SummaryDataOutputArchT.AnimalAverageT2]; % ArchT animals' lick rate for cue2
    
    % Plot control lick rate bar graph
    figure;
    bar([1 2 ], [median(Control_LR1) median(Control_LR2)]);
    hold on
    plot([1 2], [Control_LR1' Control_LR2']);
    
    % Plot ArchT lick rate bar graph
    figure;
    bar([1 2 ], [median(ArchT_LR1) median(ArchT_LR2)]);
    hold on
    plot([1 2], [ArchT_LR1' ArchT_LR2'])
    
    % Statistical comparison of lick rates
    boxstat(Control_LR1, Control_LR2, 'ControlT1', 'ControlT2', 0.05, 'paired')
    boxstat(ArchT_LR1, ArchT_LR2, 'AchTT1', 'ArchTT2', 0.05, 'paired')
    
    % Statistical comparison of lick rate differences
    Control_diff = Control_LR1-Control_LR2;
    ArchT_diff = ArchT_LR1-ArchT_LR2;
    boxstat(Control_diff, ArchT_diff, 'Control_diff', 'ArchT_diff', 0.05)
    
    % Bar plot of lick rate differences
    figure;
    bar([1 2 ], [mean(Control_diff) mean(ArchT_diff)]);
    hold on
    err1 = std(Control_diff)/sqrt(length(Control_diff));
    err2 = std(ArchT_diff)/sqrt(length(ArchT_diff));
    errorbar([1 2], [mean(Control_diff) mean(ArchT_diff)], [err1 err2],[err1 err2])

    %CDF of reaction times
    % Control animals
    figure; cdfplot(Control_RT1);
    hold on
    cdfplot(Control_RT2);
    
    % ArchT animals
    figure; cdfplot(ArchT_RT1);
    hold on
    cdfplot(ArchT_RT2);
    
    % CDF of reaction time differences    
    figure;
    cdfplot(diff_RT_Ctrl);
    hold on
    cdfplot(diff_RT_ArchT);
    
    % ROC test for the lick rate differences to cue1 and cue2 in ArchT and Control animals
    pval_archt = roc_PSTH(Hit_psth_ArchT, FA_psth_ArchT, [-5 5], [0 1.2]);
    pval_ctrl = roc_PSTH(Hit_psth_Control, FA_psth_Control, [-5 5], [0 1.2]);
    
    % Plot p-values
    figure;
    plot(1:121,pval_ctrl)
    hold on
    plot(1:121,pval_archt)
    hold on
    plot(1:121,repmat(0.05, 101))
    
    % ROC analysis of lick rate differences between control and ArchT animals
    pval_all = roc_PSTH(Hit_psth_Control-FA_psth_Control, Hit_psth_ArchT-FA_psth_ArchT, [-5 5], [0 1.2]);
    figure;
    plot(1:121,pval_all)
    hold on
    plot(1:121,repmat(0.05, 101))
    
    resdir_box_avg = fullfile(resdir_box, 'BoxPlots_06_1');
    lick_boxplot_average(SummaryDataOutputControl, SummaryDataOutputArchT, resdir_box_avg)
end

% -------------------------------------------------------------------------
function responses = getcellresp(cellids)

% VP neuron responses
cueresp = getvalue('cueresponse',cellids);
rewardresp = getvalue('rewardresponse2',cellids);
punishmentresp = getvalue('punishresponse',cellids);
omissionresp = getvalue('omissionresponse',cellids);

% Cue response
cue_e = cellids(cueresp == 1);  % activated
cue_i = cellids(cueresp == -1); % inhibited
cue_n = cellids(cueresp == 0);  % non-responsive

%reward response
rew_e = cellids(rewardresp == 1);  % activated
rew_i = cellids(rewardresp == -1); % inhibited
rew_n = cellids(rewardresp == 0);  % non-responsive

%punishment response
pun_e = cellids(punishmentresp == 1);  % activated
pun_i = cellids(punishmentresp == -1); % inhibited
pun_n = cellids(punishmentresp == 0);  % non-responsive

om_e = cellids(omissionresp == 1);
om_i = cellids(omissionresp == -1);
om_n = cellids(omissionresp == 0);

% Output
responses = struct;
responses.cue.excitation = cue_e';
responses.cue.inhibition = cue_i';
responses.cue.none = cue_n';
responses.reward.excitation = rew_e';
responses.reward.inhibition = rew_i';
responses.reward.none = rew_n';
responses.punishment.excitation = pun_e';
responses.punishment.inhibition = pun_i';
responses.punishment.none = pun_n';
responses.omission.excitation = om_e';
responses.omission.inhibition = om_i';
responses.omission.none = om_n';

% -------------------------------------------------------------------------
function add_TE_event(cellids, fieldname)
for t = 1:length(cellids)
    cellid = cellids{t};
    TE = loadcb(cellid,'TrialEvents');   % load trial events
    NumTrials = length(TE.NTrials);
    switch fieldname % we can add more cases, like 'PunishedTrials' or 'RewardedTrials'
        case 'MedialSplit'
            if ~isfield(TE, fieldname)
                TE.(fieldname) = nan(1,NumTrials);
                median_trialcount = median([1:NumTrials]);
                TE.(fieldname)(1:median_trialcount) = 1; % first median
                TE.(fieldname)(median_trialcount+1:end) = 2; % second median
                TE.(fieldname)(isnan(TE.Punishment)) = nan;
                [mouse, session, ~, ~] = cellid2tags(cellid);
                save(fullfile(getpref('cellbase', 'datapath'), mouse, session, 'TrialEvents.mat'), '-struct', 'TE');
            end
        case 'PreviousTrial1'
            if ~isfield(TE, fieldname)
                TE.(fieldname) = nan(1,NumTrials);
                for c = 2:NumTrials
                    if ~isnan(TE.AllReward(c-1)) && TE.Punishment(c) == 1% reward
                        TE.(fieldname)(c) = 1;
                    elseif ~isnan(TE.Punishment(c-1)) && TE.Punishment(c) == 1 % punishment
                        TE.(fieldname)(c) = 2;
                    elseif ~isnan(TE.Omission(c-1)) && TE.Punishment(c) == 1 % omission
                        TE.(fieldname)(c) = 3;
                    end
                end
                TE.(fieldname)(TE.Punishment~=1) = nan;
                [mouse, session, ~, ~] = cellid2tags(cellid);
                save(fullfile(getpref('cellbase', 'datapath'), mouse, session, 'TrialEvents.mat'), '-struct', 'TE');
            end
        case 'PreviousTrial2'
            if ~isfield(TE, fieldname)
                TE.(fieldname) = nan(1,NumTrials);
                for c = 2:NumTrials
                    if ~isnan(TE.AllReward(c-1)) && TE.Punishment(c) == 2% reward
                        TE.(fieldname)(c) = 1;
                    elseif ~isnan(TE.Punishment(c-1)) && TE.Punishment(c) == 2 % punishment
                        TE.(fieldname)(c) = 2;
                    elseif ~isnan(TE.Omission(c-1)) && TE.Punishment(c) == 2% omission
                        TE.(fieldname)(c) = 3;
                    end
                end
                TE.(fieldname)(TE.Punishment~=2) = nan;
                [mouse, session, ~, ~] = cellid2tags(cellid);
                save(fullfile(getpref('cellbase', 'datapath'), mouse, session, 'TrialEvents.mat'), '-struct', 'TE');
            end
        case 'AllPreviousTrial'
            if ~isfield(TE, fieldname)
                TE.(fieldname) = nan(1,NumTrials);
                for c = 2:NumTrials
                    if ~isnan(TE.AllReward(c-1)) && ~isnan(TE.Punishment(c))% reward
                        TE.(fieldname)(c) = 1;
                    elseif ~isnan(TE.Punishment(c-1)) && ~isnan(TE.Punishment(c)) % punishment
                        TE.(fieldname)(c) = 2;
                    elseif ~isnan(TE.Omission(c-1)) && ~isnan(TE.Punishment(c))% omission
                        TE.(fieldname)(c) = 3;
                    end
                end
                TE.(fieldname)(TE.Punishment~=2) = nan;
                [mouse, session, ~, ~] = cellid2tags(cellid);
                save(fullfile(getpref('cellbase', 'datapath'), mouse, session, 'TrialEvents.mat'), '-struct', 'TE');
            end
    end
    
end
