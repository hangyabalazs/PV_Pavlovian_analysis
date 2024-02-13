function [StageData] = session2stage(T, WhichStageUserInput)
%SESSION2STAGE Sorting sessions into stages.
%   SESSION2STAGE(T,WHICHSATGEUSERINPUT)Identifying the
%   WHICHSTAGGEUSERINPUT stages of the sessions of data T based on the
%   ratio of tone type 1 and type 2, presense of punishment and freewater.

%   Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hegedus.panna@koki.mta.hu

StageData = {};
for p = 1:length(T)
    T1 = T{1,p}.ToneType1;
    T2 = T{1,p}.ToneType2;
    NoPunishment = T{1,p}.NoPunishment;
    if ((T1(1)==1) & (T2(1)==0) & (NoPunishment(1)==1))==1 %& (FreeWater(1)==20);
        T{p}.('Stage') = 1;  % 'Stage 1 + 2';
    elseif ((T1(1)==0.75) & (T2(1)==0.25) & (NoPunishment(1)==1))
        T{p}.('Stage') = 2;  % 'Stage 3';
    elseif ((T1(1)==0.75) & (T2(1)==0.25) & (NoPunishment(1)==0))
        T{p}.('Stage') = 3;  % 'Stage 4';
    elseif ((T1(1)==0.6) & (T2(1)==0.4) & (NoPunishment(1)==0))
        T{p}.('Stage') = 4;  % 'Stage 5';
    elseif ((T1(1)==0.5) & (T2(1)==0.5) & (NoPunishment(1)==0))
        T{p}.('Stage') = 5;  % 'Stage 6';
    else
    end
end

if WhichStageUserInput == 0
    StageData = T;
else
    v = 1; % for populating StageData
    for w = 1:length(T)
        if T{1,w}.Stage == (WhichStageUserInput)
            StageData{1,v} =  T{1,w};
            v = v+1;
        end
    end
end