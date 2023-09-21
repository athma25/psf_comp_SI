function printPdf(fig,fName)

% Saves fig as fName.pdf in working directory

set(fig,'Units','Inches');
fPos=get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[fPos(3), fPos(4)]);
print(fig,sprintf('%s/%s',pwd,fName),'-dpdf','-r400');

end
