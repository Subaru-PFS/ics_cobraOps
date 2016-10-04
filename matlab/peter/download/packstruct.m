function output=packstruct(varargin)
%PACKSTRUCT - packs variables into structure
%
%Examples: Given structure myStruct, with fields a,b,c, & d
%
% (1) packstruct(myStruct)   %assign fields to variables
% 
%         ans =
%         'a = myStruct.a;b = myStruct.b;c = myStruct.c;d = myStruct.d;     
% 
%Usage:
%
%     assigns=packstruct(InputStructure)
% 
%     in:
%         field variables
% 
%     out:
%         InputStructure: A structure
%
%Typical usage:
% 
%     Structure=packstruct(a,b,c,d,...);
% 
% props to Matt Jacobson

for jj=1:length(varargin);
  fname = inputname(jj);
%  verbosity = ['output.' fname ' = ' fname ';'];
%  disp(verbosity);
  cmd = ['output.' fname ' = varargin{' num2str(jj) '};'];
  eval(cmd);
end
