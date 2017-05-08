function output = NumArr2mMapStr(NumArr,mapSize)
% converts numerical array into a string for input into xml
% database file
%
% usage: output = NumArr2mMapStr(vector, mapSize=100)
    
if ~exist('mapSize','var'),mapSize=100;,end;

NumArr = [NumArr(:); zeros(mapSize+2-length(NumArr),1)];

mMapStr = '';

for ii=1:length(NumArr)
    mMapStr = [mMapStr num2str(NumArr(ii)) ','];
end


output = mMapStr;

end
