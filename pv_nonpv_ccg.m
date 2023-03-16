function pv_nonpv_ccg(cellids,resdir,issave)
%PV_NONPV_CCG(CELLIDS, RESDIR, ISSAVE)   Cross-correlation analysis.
%   PV_NONPV_CCG calculates cross-correlations for PV(+) and non-tagged
%   neurons in 50 ms windows. CCG results are saved for non-tetrode pairs
%   and tetrode pairs separately.
%
%   See also CCG.

%   Balazs Hangya, Panna Hegedus
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu

% Directories
if ~isfolder(resdir)
    mkdir(resdir)
end

% Input argument check
narginchk(0,3);
if nargin < 2
    issave = true;   % default saving behavior
end
if nargin < 1
    pv_cells = select_pv_cells([]);
else
    pv_cells = cellids;
end
pv_cells = pv_cells';

loadcb

% CCG
for p = 1:length(pv_cells)
    c_pv_cellid = pv_cells{p};
    [animal, session, ~, ~]=cellid2tags(c_pv_cellid);
    session_cells = findcell('rat', animal, 'session', session);
    non_pv_cells = setdiff(session_cells, c_pv_cellid);
    for n = 1:length(non_pv_cells)
        segfilter = 'stim_excl_vp';
        filterinput = {'light_activation_duration',[-5 5],'margins',[0 0]};
        ccg({c_pv_cellid non_pv_cells{n}},0.05,'whichcells','allpairs','resdir',resdir,...
            'segfilter',segfilter,'filterinput',filterinput,...
            'minspikeno',100,'maxspikeno',10000,'issave',issave);
        ccg_grouping_VP({c_pv_cellid non_pv_cells{n}},'whichcells','allpairs','resdir',resdir,...
            'segfilter',segfilter,'filterinput',filterinput,...
            'minspikeno',100,'maxspikeno',10000,'issave',issave);
    end
end
