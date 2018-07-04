function sldplot(xcell,ycell,opts,labels,legendopt)

if nargin <=2
    opts = {};
    legendopt = {};
    labels = {};
end

f = gcf;                       
ax = gca;

nol = length(xcell);
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

