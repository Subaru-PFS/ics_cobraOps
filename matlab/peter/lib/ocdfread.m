% ocdfread.m: a hack to read cdf's from octave
% usage: [data info] = ocdfread('file');

function [data info] = ocdfread(filename)

  stdout = system(["cdfread.pl " filename]);
  eval(stdout);
  endfunction