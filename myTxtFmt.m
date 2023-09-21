function myTxtFmt(t,s,tex)

% t: text object
% s: font size
% tex<-(0/1): 1 for latex interpretation

t.FontSize=s;
if tex==1
	t.Interpreter='Latex';
end
