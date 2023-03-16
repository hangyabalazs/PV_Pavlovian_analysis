function pvals = roc_PSTH(hit_psth, fa_psth, win, testwin)
% ROC_PSTH performs ROC analysis.
% ROC_PSTH(HIT_PSTH, FA_PSTH, TAG, WIN, TESTWIN, STAGE) performs ROC
% analysis of input data (HIT_PSTH and FA_PSTH) of length WIN in a test window
% TESTWIN.

% See also ROCAREA

% Balazs Hangya, Panna Hegedüs
% hangya.balazs@koki.hu
% 13/03/2023

dt = 0.001;   % resolution, in seconds
time = win(1):dt:win(2);   % time vector
t1_inx = find(time==testwin(1));
t2_inx = 6201;

x = hit_psth(:,t1_inx:10:t2_inx); % take every 10th timepoint
y = fa_psth(:,t1_inx:10:t2_inx);

numSample = size(x,2);
[SE,pvals] = deal(nan(1,numSample));

% Perform ROC analysis
for n = 1:numSample
    [D, pvals(n), SE(n)] = rocarea(x(:,n),y(:,n), 'bootstrap', 200, 'display', false);
end