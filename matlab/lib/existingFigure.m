function h = existingFigure(name)

if isempty(findobj('type','figure','name',name))
    h = figure('name',name);
else
    h = findobj('type','figure','name',name);
    figure(h)
end

return;