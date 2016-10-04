function [Text] = Array2Text(Array)

size_array = length(Array);

Text = [num2str(size_array) ',100'];

for ii=1:size_array
    Text = [Text ',' num2str(Array(ii))];
end

filler = 100-size_array;

for jj = 1:filler
    Text = [Text ',0.0'];
end

Text = [Text ','];

