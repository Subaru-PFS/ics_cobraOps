function output = plotFibers(data,varargin)
hold on
% plot(  [1.5082e+03 + 5.2461e+02i;  1.5228e+03 + 1.5450e+03i], 'rx')
% hold on


keyboard;
pos = [];
posid = [];
for ii=data.pId
    fldID = ['pId' num2str(ii)];
    if isfield(data.(fldID),'CLpos')
        posID = 'CLpos';
    elseif isfield(data.(fldID),'CCDpos')
        posID = 'CCDpos';
    end
    
    pos = [pos data.(fldID).(posID)];   
    text(real(data.(fldID).(posID)),imag(data.(fldID).(posID)),num2str(ii));
end

output = plot(pos,varargin{:});

% title(data.Name,'Interpreter','none')
axis equal

return