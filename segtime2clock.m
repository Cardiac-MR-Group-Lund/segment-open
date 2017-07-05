function clockstr = segtime2clock()
%Returns acquisition time as a string

%Likely this could rewritten to a single line code...

%Sebastian Bidhult

global SET NO

sectotal = SET(NO).AcquisitionTime;

hour = sectotal/3600;
minutes = 60*(hour-floor(hour));
sec = 60*(minutes-floor(minutes));

hour = floor(hour);
minutes = floor(minutes);
sec = floor(sec);

clockstr = [num2str(hour) ':' num2str(minutes) ':' num2str(sec)];