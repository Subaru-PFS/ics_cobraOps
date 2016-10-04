function output = NumArr2mMapStr(NumArr,mapSize)

if ~exist('mapSize','var'),mapSize=100;,end;

NumArr = [NumArr; zeros(mapSize+2-length(NumArr),1)];

mMapStr = '';

for ii=1:length(NumArr)
%     if ii==length(NumArr)
%         mMapStr = [mMapStr num2str(NumArr(ii))];
%     else
        mMapStr = [mMapStr num2str(NumArr(ii)) ','];
%     end
end


output = mMapStr;

end
