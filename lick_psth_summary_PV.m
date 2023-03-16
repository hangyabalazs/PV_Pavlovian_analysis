function [SummaryDataOutput, Hit_allpsth, FA_allpsth] = lick_psth_summary_PV(animals, root, resdir, tag)
%LICK_PSTH_SUMMARY_PV  Average lick PETH.
%   LICK_PSTH_SUMMARY_PV(ANIMALS, ROOT, RESDIR, TAG) plots average lick
%   PETH (beam break time stamps aligned to an event) for each animal
%   defined in ANIMALS. Last 10 sessions are used in case of each animal to
%   calculate average PETH. Behavioral data is loaded from ROOT, PETHs are
%   saved to RESDIR. TAG indicates experimental group(control or ArchT).
%
%   [SummaryDataOutput, Hit_allpsth, FA_allpsth] =
%   lick_psth_summary_PV(animals, root, resdir, tag) generates output data
%   structure (SummaryDataOutput) containing the ID of the animal, the
%   training stage, tag, lick rate and reaction time to cue1 and cue2.
%   Hit_allpsth and FA_allpsth contain lick peths for cue1 and cue2 for
%   each animal.

%   See also ULTIMATE_PSTH.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu
%   20-Dec-2019

if isempty(animals) || ~exist('animals','var')
    if strcmp(tag, 'ArchT')
        animals = {'PVINH12','PVINH21','PVINH22','PVINH31','PVINH41','PVINH52'};
    elseif strcmp(tag, 'Control')
        animals = {'PVINH11','PVINH14','PVINH33','PVINH34','PVINH42','PVINH44','PVINH51','PVINH54'};
    elseif isempty(tag)
        animals = {'PVINH12','PVINH21','PVINH22','PVINH31','PVINH41','PVINH52','PVINH11','PVINH14','PVINH33','PVINH34','PVINH42','PVINH44','PVINH51','PVINH54'};
    end
end

    choosecb('PV_opto_pavlovian');
    fullpth = getpref('cellbase', 'datapath');


if isempty(root) || ~exist('root','var')
    root = 'C:\Users\Victoria\Desktop\KOKI\Data';
end

if ~isfolder(resdir)
    mkdir(resdir)
end

prompt = {'0 - all stages,  1 - stage 1,  2 - stage 2,  3 - stage 3,  4 - stage 4,  5 - full task only'};
dlgtitle = 'Which stage would you like to analyse? ';
dims = [3 135];
definput = {'', 'hsv'};
WhichStageUserInput = str2double(inputdlg(prompt,dlgtitle,dims,definput));
if isempty(WhichStageUserInput) || ~exist('WhichStageUserInput', 'var')
    error 'More information is needed. Please specify which stage you would like to analyse!'
end
WhichStageUserInputString = num2str(WhichStageUserInput); % for naming saved plots later


BigData ={};
p = 0; % for populating full_sessions array

for i = 1:length(animals)  % for each animal
    c_animal = animals{i};

    subfolders_gross= dir(fullfile(fullpth,c_animal));
    subfolders = subfolders_gross(~ismember({subfolders_gross(:).name},{'.','..'})); % to filter out unwanted elements
    
    % Print folder names to command window.
    sessionIDs = {}; % subfolders
    for k = 1 : length(subfolders)
        fprintf('Sub folder #%d = %s\n', k, subfolders(k).name);
        sessionIDs{k}=subfolders(k).name;        
    end
    
        % all existing session folders within 1 animal's folder:
    for m= 1:length(sessionIDs)       
        T=load(fullfile(fullpth,c_animal,sessionIDs{m}, 'TE.mat'));     % The TE file containing the session data is uploaded    

        p=p+1;
        newstruct1 = 'NameOfAnimal';
        newstruct2 = 'SessionName';
        newstruct4 = 'ArchTORControl';
        BigData{p} = T; %information from all experiments is stored
        BigData{p}.(newstruct1) = c_animal;
        BigData{p}.(newstruct2) = sessionIDs{m};
        BigData{p}.(newstruct4) = tag;

    end
end

% Time window
wn = [-5 5];   % in seconds
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector

% PETH
StageData = session2stage(BigData, WhichStageUserInput);

if WhichStageUserInput ~= 1
    [SummaryDataOutput] = lick_boxplot(StageData, root, resdir, WhichStageUserInput, animals);
end

NumSessions = size(StageData,2);   % number of sessions, columns
[Hit_allpsth, FA_allpsth, Hit_allbinrast, FA_allbinrast, Hit_meanFR, FA_meanFR, Hit_maxFR, FA_maxFR] = deal([]);

for iS = 1:NumSessions
    SN = StageData{iS}.SessionName;
    cellid = {StageData{1,iS}.NameOfAnimal, SN}
    [spsth_hit, spsth_fa, binrast_hit, binrast_fa] = main(cellid,wn,dt,tag, WhichStageUserInputString, resdir);
    Hit_allpsth = [Hit_allpsth; spsth_hit]; 
    FA_allpsth = [FA_allpsth; spsth_fa];
    Hit_allbinrast = [Hit_allbinrast; nanmean(binrast_hit, 1)];
    FA_allbinrast = [FA_allbinrast; nanmean(binrast_fa,1)];
    Hit_meanFR = [Hit_meanFR, mean(spsth_hit(:,5001:5201))];
    FA_meanFR = [FA_meanFR, mean(spsth_fa(:,5001:5201))];
    Hit_maxFR = [Hit_maxFR, max(spsth_hit(:,5001:5201))];
    FA_maxFR = [FA_maxFR, max(spsth_fa(:,5001:5201))];
end

% Plot & save
figure;
G = gcf; % current figure handle
green = [51 204 51] / 255;   % colors for plotting
red = [216 41 0] / 255;
errorshade(time,nanmean(Hit_allpsth),nanstd(Hit_allpsth)/sqrt(size(Hit_allpsth,1)),...
    'LineColor',green,'ShadeColor',green)
hold on
errorshade(time,nanmean(FA_allpsth),nanstd(FA_allpsth)/sqrt(size(FA_allpsth,1)),...
    'LineColor',red,'ShadeColor',red)

fnm = fullfile(resdir,strcat('average_lickPSTH_', tag, '_Stage_', WhichStageUserInputString,'.fig'));
set(G, 'renderer', 'painters')
fnm2 = fullfile(resdir,strcat('average_lickPSTH_', tag, '_Stage_', WhichStageUserInputString,'.eps'));
fnm3 = fullfile(resdir,strcat('average_lickPSTH_', tag, '_Stage_', WhichStageUserInputString,'.jpg'));
saveas(G,fnm)
saveas(G,fnm2)
saveas(G,fnm3)
close(G)


% -------------------------------------------------------------------------
function [spsth_hit, spsth_fa, binrast_hit, binrast_fa] = main(cellid, wn, dt, tag, WhichStageUserInputString, resdir)

% Filter input
filterinput_hit = 'TrialType==1';
filterinput_fa = 'TrialType==2';

% Calcualte lick PSTH
[~, spsth_hit, ~, ~, binrast_hit, ~] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.08,'event_filter','custom','filterinput',filterinput_hit,...
    'isadaptive',0,'maxtrialno',Inf);
[~, spsth_fa, ~, ~, binrast_fa, ~] = ...
    ultimate_psth(cellid, 'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.08,'event_filter','custom','filterinput',filterinput_fa,...
    'isadaptive',0,'maxtrialno',Inf);

%--------------------
% % Lick raster , singe plots for double-checking --> average plot
H = figure;
viewlick(cellid,'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav',...
    'ShowEvents',{{'StimulusOn' 'StimulusOff'}},...
    'Partitions','#TrialType','window',wn)
%maximize_figure(H)
cellidt = [cellid{1} '_' cellid{2}];
fnm = fullfile(resdir(1:end-13), 'Single plots', [cellidt '_' tag '_Stage' WhichStageUserInputString  '_LICK.fig']);
saveas(H,fnm)
fnm = fullfile(resdir(1:end-13), 'Single plots', [cellidt '_' tag '_Stage' WhichStageUserInputString  '_LICK.jpg']);
saveas(H,fnm)
close(H)