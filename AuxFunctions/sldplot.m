function sldplot(xcell,ycell,opts,labels,legendopt)
%SLDPLOT Control plot with a slider.
%
%   SLDPLOT(XCELL,YCELL) plots the arrays listed in XCELL and YCELL. XCELL
%   and YCELL must be of same length. Each element of XCELL must be a
%   N-by-1 vector. Each element of YCELL must be a M-by-N matrix. 
%
%   SLDPLOT(XCELL,YCELL,OPTS) passes the optional argument OPTS for
%   customize the plot. OPTS must be a length(XCELL)-by-1 cell array of
%   cell arrays.
%
%   SLDPLOT(XCELL,YCELL,OPTS,LABELS) passes the labels of each slide step
%   to title the interactive figure. LABELS must be a M-by-1 cell array of
%   strings.
%
%   SLDPLOT(XCELL,YCELL,OPTS,LABELS,LEGENDOPT) plot the legend. LEGENDOPT
%   is a cell array with the argument of the legend() command.
%
%   Examples:See the Furnace_wall.m file.
% 
%   ============================================================
%   Author: ataide@peq.coppe.ufrj.br
%   homepage: github.com/asanet
%   Contact me for help/personal classes!


if nargin <= 2
    opts = {};
    legendopt = {};
    labels = {};
end

nol = length(xcell);
if nol ~= length(ycell)
    error('Incomplete data set. X and Y must be in pairs.')
end

f = gcf;                       
ax = gca;

sz = length(ycell{1});

ax.Units = 'normalized';
ax.Position([2 4]) = [ax.Position(2) + 0.075, ax.Position(4) - 0.075];
uicontrol('Parent',f,'style','slide','units','normalized','position',[0.1 0.025 0.8 0.05 ], ...
                     'min',1,'max',sz,'Value',1,'sliderstep',[1/sz 1/sz],'Tag','sld'); 

h = gobjects(nol,1);                 
while ishandle(f)
    i = fix(get(f.Children(1),'Value'));
    h(1) = plot(xcell{1},ycell{1}(i,:));
    set(h(1),opts{1}{:})
    hold(ax,'on')
    for k = 2:nol
        h(k) = plot(xcell{k},ycell{k}(i,:));
        set(h(k),opts{k}{:})
    end
    hold(ax,'off')
    legend(legendopt{:})
    if ~isempty(labels)
        title(labels{i})
    end
    ax.NextPlot = 'replaceChildren';
    waitfor(f.Children(1),'Value')
end

                 
end

