% Plot ImageId structures currently in workspace

S = whos;
A = regexp(cellstr([S.name]),'ImageId_\d*','match');
CW = colormap(hsv(length(A{1})));

for ii = 1:length(A{1})
    thisvar = eval(char(A{1}(ii)));
    plotFibers(thisvar,'x','MarkerEdgeColor',CW(ii,:))
    hold on
end