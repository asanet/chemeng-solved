function f = figured(int)

if nargin == 1
    f = figure(int);
else
    f = figure;
end

set(f,'Units','normalized','PaperOrientation','landscape','PaperPositionMode','auto',...
                                        'Tag','main','Position',[0.3 0.25 0.4 0.5]); 
axes('YGrid','on','Box','on','FontSize',16,'NextPlot','replacechildren');

mh1 = uimenu('Parent',f,'Label','Additions');
uimenu(mh1,'Label','Save2PDF','Enable','on','Tag','saveinpdf','Accelerator','B');
h = guihandles(f);
set(h.saveinpdf,'Callback',@save2pdfemb);


guidata(f,h)
end

function save2pdfemb(ho,~)
h = guidata(ho);

[FileName,PathName] = uiputfile({'*.pdf','Portable Document File (.pdf)'}, 'Save Figure as...');
if FileName~=0
    save2pdflocal([PathName FileName],800,h.main)
end
end

function save2pdflocal(pdfFileName,dpi,handle)

% Backup previous settings
prePaperType = get(handle,'PaperType');
prePaperUnits = get(handle,'PaperUnits');
preUnits = get(handle,'Units');
prePaperPosition = get(handle,'PaperPosition');
prePaperSize = get(handle,'PaperSize');

% Make changing paper type possible
set(handle,'PaperType','<custom>');

% Set units to all be the same
set(handle,'PaperUnits','inches');
set(handle,'Units','inches');

% Set the page size and position to match the figure's dimensions
% paperPosition = get(handle,'PaperPosition');
position = get(handle,'Position');
set(handle,'PaperPosition',[0,0,position(3:4)]);
set(handle,'PaperSize',position(3:4));

% Save the pdf (this is the same method used by "saveas")
print(handle,'-dpdf',pdfFileName,sprintf('-r%d',dpi))

% Restore the previous settings
set(handle,'PaperType',prePaperType);
set(handle,'PaperUnits',prePaperUnits);
set(handle,'Units',preUnits);
set(handle,'PaperPosition',prePaperPosition);
set(handle,'PaperSize',prePaperSize);

end
