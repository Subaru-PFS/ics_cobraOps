README file for allan.m / allan_overlap.m / allan_modified.m
MH Nov2010 v2.x

Allan deviation analysis often involves the use of large data sets. allan.m will plot the entire data set by default, and this may take some time (minutes) for large data sets. To improve the performance of MATLAB when plotting large data sets, I recommend installing the m-file "dsplot.m" by Jiro Doke. You can download this file from the File Exchange:

	http://www.mathworks.com/matlabcentral/fileexchange/15850

If this file is present on your MATLAB path then allan.m will use it. Note that you can suppress the plotting of the data set by using VERBOSE < 2. For more help using allan.m, type "help allan" at the MATLAB prompt.

I welcome your feedback about this code and your experience with it- user feedback has been invaluable in achieving the current level of performance.

M.A. Hopcroft
mhopeng@ml1.net
18Oct2010
