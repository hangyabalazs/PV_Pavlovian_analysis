function [AllAverage] = lick_boxplot(StageData, resdir, WhichStageUserInput,TrialType)
%LICK_BOXPLOT Compare the number of licks for the two cues of the Pavlovian task.
%   LICK_BOXPLOT(STAGEDATA, RESDIR, WHICHSTAGEUSERINPUT, TRIALTYPE)Creates
%   and saves boxplots to RESDIR, comparing the number of licks for two
%   cues defined by TRIALTYPE from the stage of the task defined in
%   WHICHSTAGEUSERINPUT based on the extracted data STAGEDATA. Average
%   lickrates, licknumbers and reaction times are returned in the struct
%   ALLAVERAGE.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu

% Initialize
AllAverage = {};

resdir_boxplot = fullfile( resdir, 'BoxPlots_06_1'); % subfolder for saving results into
if ~isfolder(resdir_boxplot)
    mkdir(resdir_boxplot);
end
T = StageData;

for i = 1:length(T)
    
    % Extracting Licks
    T1_unfiltered = [T{i}.StimulusOn];
    T1_inx = T{i}.(TrialType)==1;
    T1_filtered = T1_unfiltered(T1_inx);
    BeamBreak = [T{i}.LickIn];
    BeamBreak_T1 = BeamBreak(T1_inx);
    
    NumTrials1 = length(BeamBreak_T1);
    [T1_anticipatory, T1_RT] = deal(nan(1,NumTrials1));
    
    Lick_T1 = ((nansum(T{i}.Hit(T{i}.(TrialType)==1))) + (nansum(T{i}.LickOmission(T{i}.(TrialType)==1))) + (nansum(T{i}.FalseAlarm(T{i}.(TrialType)==1)))) / sum(T{i}.(TrialType) == 1);
    Lick_T2 = ((nansum(T{i}.FalseAlarm(T{i}.(TrialType)==2))) + (nansum(T{i}.LickOmission(T{i}.(TrialType)==2))) + (nansum(T{i}.Hit(T{i}.(TrialType)==2)))) / sum(T{i}.(TrialType) == 2);
    
    % Calculating reaction time to cue1
    for n = 1:NumTrials1
        T1_bb = BeamBreak_T1{n};
        inx = (T1_filtered(n)+0.6 <= T1_bb) & (T1_bb < T1_filtered(n)+1.1);
        T1_anticipatory(n) = sum(double(inx))/1.1; %how many ones
        if ~isempty(T1_bb)
            T1_RT(n) = min(T1_bb(T1_filtered(n) <= T1_bb))- T1_filtered(n);
        end
    end
    
    T2_unfiltered = [T{i}.StimulusOn];
    T2_inx = T{i}.(TrialType)==2;
    T2_filtered=T2_unfiltered(T2_inx);
    
    if ~isempty(T2_filtered)
        BeamBreak = [T{i}.LickIn];
        BeamBreak_T2 = BeamBreak(T2_inx);
        
        NumTrials2 = length(BeamBreak_T2);
        [T2_anticipatory, T2_RT] = deal(nan(1,NumTrials2));
        
        %Calculating reaction time to cue2
        for m = 1:NumTrials2
            T2_bb = BeamBreak_T2{m};
            inx = (T2_filtered(m)+0.6 <= T2_bb) & (T2_bb < T2_filtered(m)+1.1);
            T2_anticipatory(m) = sum(double(inx))/1.1; %how many ones
            if ~isempty(T2_bb)
                T2_RT(m) = min(T2_bb(T2_filtered(m) <= T2_bb))- T2_filtered(m);
            end
        end
        NameOfAnimal = T{i}.NameOfAnimal;
        tag = T{i}.ArchTORControl;
        sessionID = T{i}.SessionName;
        StageString = num2str(WhichStageUserInput);
        
        % Compater anticipatory licks between the 2 cues
        boxstat(T1_anticipatory,T2_anticipatory,'Reward cue','Punishment cue', 0.05)
        B = gcf;
        
        % Save figs
        fnm = fullfile(resdir_boxplot,strcat('boxplot_lickPSTH_', NameOfAnimal, '_', tag, '_', sessionID, '_Stage', StageString, '.fig')); %%%
        set(B, 'renderer', 'painters')
        fnm2 = fullfile(resdir_boxplot,strcat('boxplot_lickPSTH_', NameOfAnimal, '_', tag, '_', sessionID, '_Stage ', StageString, '.eps'));
        fnm3 = fullfile(resdir_boxplot,strcat('boxplot_lickPSTH_', NameOfAnimal, '_',tag, '_', sessionID, '_Stage ', StageString, '.jpg'));
        
        saveas(B,fnm)
        saveas(B,fnm2)
        saveas(B,fnm3)
        close(B)
        %-------------------------
        
        % Calculate average likc rates and reaction times
        AnimalAverageT1 = (sum(T1_anticipatory))/NumTrials1;
        AnimalAverageT2 = (sum(T2_anticipatory))/NumTrials2;
        AnimalAverageRT1 = nanmean(T1_RT);
        AnimalAverageRT2 = nanmean(T2_RT);
         
        newstruct1 = 'NameOfAnimal';
        newstruct2 = 'Stage';
        newstruct3 = 'Tag';
        newstruct4 = 'AnimalAverageT1';
        newstruct5 = 'AnimalAverageT2';
        newstruct6 = 'AnimalAverageRT1';
        newstruct7 = 'AnimalAverageRT2';
        newstruct8 = 'LickT1';
        newstruct9 = 'LickT2';
        
        AllAverage(i).(newstruct1) = NameOfAnimal;
        AllAverage(i).(newstruct2) = StageString;
        AllAverage(i).(newstruct3) = tag;
        AllAverage(i).(newstruct4) = [AnimalAverageT1];
        AllAverage(i).(newstruct5) = [AnimalAverageT2];
        AllAverage(i).(newstruct6) = [AnimalAverageRT1];
        AllAverage(i).(newstruct7) = [AnimalAverageRT2];
        AllAverage(i).(newstruct8) = [Lick_T1];
        AllAverage(i).(newstruct9) = [Lick_T2];
    end
    
end