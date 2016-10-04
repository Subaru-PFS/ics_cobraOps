function output=loadmats(file_expr,mypath)
% usage OUTPUT = loadmats(FILE_EXPR, MYPATH)
% example output = loadmats('*.mat','.');
% load every matfile that matches FILE_EXPR in MYPATH into OUTPUT
%
% if you want OUTPUT unpacked, run "eval(unpackstruct(OUTPUT));"

if ~exist('file_expr','var'), file_expr = '*.mat'; end;
if ~exist('mypath','var')   , mypath    = './'   ; end;

files = dir2cell(fullfile(mypath,file_expr));

output = struct();

for jj = 1:length(files)
  contents = load(files{jj});
  FieldName = fields(contents);
  for kk = 1:length(FieldName)
    output = setfield(output, FieldName{kk}, contents.(FieldName{kk}));
  end
end
