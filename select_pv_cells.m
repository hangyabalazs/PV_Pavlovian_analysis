function cellids = select_pv_cells(cellids)
%SELECT_PV_CELLS   Select VP neurons from CellBase.
%   SELECT_PV_CELLS(CELLIDS) selects parvalbumin(+) (PV-)cells that are from 
%   CELLIDS (default, all cell IDs in CellBase). Well-isolated units 
%   (ID > 20, L_ratio < 0.15) are returned.
%
%   See also PV_PAVLOVIAN_ANALYSIS_MAIN.

%   Balazs Hangya and Panna Hegedus
%   Laboratory of Systems Neuroscience
%   Institute of Experimental Medicine, Budapest, Hungary


% Input arguments
narginchk(0,1)
if nargin < 1
    cellids = [];
end

% Load CellBase
load(getpref('cellbase','fname'));

% List of cellIDs
if isempty(cellids)
    Lratio = getvalue('Lr_PC');
    ID = getvalue('ID_PC');
    ispv = getvalue('ispv');   % select PV cells
    ptinx = ispv & ID > 20 & Lratio < 0.15;   % good clusters
    cellids = CELLIDLIST(ptinx);
else
    if isnumeric(cellids)  % 'cellids' can be index set or list of cell IDs
        cellids = CELLIDLIST(cellids);
    else
        if ischar(cellids)
            cellids = {cellids};   % only one cell ID
        end
    end
end
cellids = cellids(:)';   % convert to row vector