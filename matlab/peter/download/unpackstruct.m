function assigns=unpackstruct(varargin)
%UNPACKSTRUCT - returns a string that, if executed, would 
%assign fields of a structure to individual variables of the same name 
%
%Examples: Given structure myStruct, with fields a,b,c, & d
%
% (1) unpackstruct(myStruct)   %assign fields to variables
% 
%         ans =
%         'a = myStruct.a;b = myStruct.b;c = myStruct.c;d = myStruct.d;     
% 
%Usage:
%
%     assigns=unpackstruct(InputStructure)
% 
%     in:
%         InputStructure: A structure
% 
%     out:
%         assigns: a text string containing the commands (see Examples above) 
%
%Typical usage:
% 
%     eval(unpackstruct(InputStructure));
% 
% props to Matt Jacobson

nn=length(varargin);
S=varargin{1};
    
fields=fieldnames(S);
sname=inputname(1); 

cmds = cellfun(@(f) [f ' = ' sname '.' f ';'],fields,'uniformoutput',0);
assigns = '';
for jj=1:length(cmds)
  assigns = strcat(assigns,cmds{jj});
%  disp(cmds{jj});
end
