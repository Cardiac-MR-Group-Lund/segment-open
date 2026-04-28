function clockstr = segtime2clock()
%Returns acquisition time as a string

%Sebastian Bidhult, updated 2023-09-19 by Jelena

global SET NO %#ok<*GVMIS> 

sectotal = SET(NO).AcquisitionTime;

% Calculate hours, minutes, and seconds from sectotal
hour = floor(sectotal / 3600);
minutes = floor(mod(sectotal, 3600) / 60);
sec = floor(mod(sectotal, 60));

% Create the clock string
clockstr = sprintf('%02d:%02d:%02d', hour, minutes, sec);
