function [H Wp] = boxstat(d1,d2,l1,l2,alpha,str)
%BOXSTAT   Box plot with statistics.
%   [H WP] = BOXSTAT(D1,D2,L1,L2,ALPHA) plots box plot (H, figure handle)
%   of data sets D1 and D2 using x axis labels L1 and L2. It performs
%   Mann-Whitney U-test with significance level ALPHA (default is 0.05) and
%   returns p value (WP).
%
%   [H WP] = BOXSTAT(D1,D2,L1,L2,ALPHA,'PAIRED') performes a paired
%   non-parametric test (Wilcoxon's signed rank test) instead of the
%   Mann-Whitney U-test.
%
%   See also BOXPLOT.

% Input argument check
error(nargchk(4,6,nargin))
if nargin < 6
    str = 'nonpaired';
end
if nargin < 5 || isempty(alpha)
    alpha = 0.05;   % default significance level
end
d1 = d1(:)';   % row vector format
d2 = d2(:)';

% Box plot
switch str
    case 'nonpaired'
        H = figure;
        boxplot([d1 d2],[zeros(size(d1)) ones(size(d2))],'labels',[{l1} {l2}]);       
        [Wp, Wh] = ranksum(d1,d2,'alpha',alpha);   % Mann-Whitney U-test
        setmyplot_balazs
    case 'paired'
        [Wp, Wh] = signrank(d1,d2,'alpha',alpha);   % Wilcoxon signed rank test
        H = figure;
        line([ones(size(d1)); ones(size(d2))*2],[d1; d2],'Color',[0.8,0.8,0.8,0.8])
        hold on
        bar([1; 2],[nanmean(d1(~isnan(d1)&~isnan(d2))); nanmean(d2(~isnan(d1)&~isnan(d2)))])
        xlim([0.5,2.5])
        xticks([1,2])
        xticklabels({l1,l2})
        setmyplot_balazs
    otherwise
        error('boxstat:inputArg','Unsupported input argument.')
end
if Wh
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp),'Color',clr,'Horizontalalignment','center')