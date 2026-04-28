function paramstr = makeunitstring(paramstr, unitstr) %skiptranslation
% function to put string together with parameter sting and unit string
% paramstr = makeunitstring(paramstr, unitstr)

%Jelena Bock

paramstr = sprintf('%s [%s]',paramstr,unitstr);
