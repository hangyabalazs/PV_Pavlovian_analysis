function H = Hindex_distribution(ntcells,taggedcells)
%HINDEX_DISTRIBUTION   Plot H-index distribution.
%   H = HINDEX_DISTRIBUTION(NTCELLS,TCELLS) plots tagging (H-) index
%   distribution for a population of tagged (TCELLS) and non-tagged cells
%   (NTCELLS).
%
%   See also SALT.

%   Balazs Hangya
%   hangya.balazs@koki.hu
%   12-14-2020

% H-index
ht = getvalue('Hindex',taggedcells);
hnt = getvalue('Hindex',ntcells);

% H-index histogram
H = figure;
histogram(hnt,70,'EdgeColor',[0.7 0.7 0.7],'FaceColor',[0.7 0.7 0.7])
hold on
histogram(ht,1,'EdgeColor',[0.1 0.1 0.75],'FaceColor',[0.1 0.1 0.75])