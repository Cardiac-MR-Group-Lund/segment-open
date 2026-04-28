function outstr = mysplitstring(str) %skiptranslation
% split string into 2 parts

numparts = 2;
% calculate approximate length of new string
strlen = ceil(length(str)/numparts);
str = strsplit(str);
% length of each substring
substrlen = cellfun(@length,str);
% add 1 to account for spaces
substrlen = substrlen +1;
%cumulative length
cumlen = cumsum(substrlen);
% find the length closest to 0.5 of original string length
lendiff = abs(cumlen - strlen);
minind = find(lendiff == min(lendiff));
% create new array with splited string through strjoin 
outstr = repmat({''},1,numparts);
% first part
outstr{1} = strjoin(str(1:minind));
% second part
outstr{2} = strjoin(str(minind+1:end));

