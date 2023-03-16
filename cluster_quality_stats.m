function cluster_quality_stats(cellids, resdir, str)
%CLUSTER_QUALITY_STATS   ID and L-ratio statistics.
%   CLUSTER_QUALITY_STATS(CELLIDS,RESDIR,STR) calculates and saves median
%   and SE of median of Isolation Distance and L-ratio for CELLIDS. These
%   parameters are stored in CellBase by LRATIO2. STR string input is used
%   to generate output file name.
%
%   See also LRATIO2.

%   Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu
%   10-Dec-2020

% ID and L-ratio (saved in CellBase using Lratio2.m)
Lr = getvalue('Lr_PC',cellids);
ID = getvalue('ID_PC',cellids);

% Median, SE
median_Lr = median(Lr);
SE_Lr = se_of_median(Lr);

median_ID = median(ID);
SE_ID = se_of_median(ID);

% Save
if ~isfolder(resdir)
    mkdir(resdir)
end
save(fullfile(resdir,[str '_cluster_quality_variables.mat']),'Lr','ID',...
    'median_Lr','SE_Lr','median_ID','SE_ID');