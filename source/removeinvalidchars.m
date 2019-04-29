function outstri = removeinvalidchars(stri)
%Removes invalid characters for the license code

%Einar Heiberg

%This code can be written more elegant, but who cares, it works and is fast
%enough.

valid = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!#$@£{}[]()+-:,;<>|';

logvec = false(size(stri));
for loop = 1:length(stri)
  if sum(stri(loop)==valid)>0
    logvec(loop) = true;
  end;
end;
outstri = stri(logvec);
