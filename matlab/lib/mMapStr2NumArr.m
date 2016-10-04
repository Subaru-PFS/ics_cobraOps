function output = mMapStr2NumArr(mMapStr)

mMapCellArr = regexp(mMapStr,',','split');

mMapNumArr = str2num(char(mMapCellArr));

output = mMapNumArr;

end
