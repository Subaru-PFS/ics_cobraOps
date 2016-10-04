function output = getLegendHandle(figureHandle)

h = get(figureHandle,'children');
hLeg = [];
for k = 1:length(h)
    if strcmpi(get(h(k),'Tag'),'legend')
        hLeg = h(k);
        break;
    end
end

output = hLeg;
end