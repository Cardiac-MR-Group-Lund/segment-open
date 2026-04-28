function stri = removeforbiddenchars(stri)
%Remove forbidden characters from string, suitable for creating valid filenames.
%
%See also REMOVEINVALIDCHARS

%Einar Heiberg

%stri = upper(stri);
stri = deblank(stri);
stri(stri==' ')='_';
stri(stri=='Å') = 'A';
stri(stri=='Ä') = 'A';
stri(stri=='Ö') = 'O';
stri(stri=='Ü') = 'U';
stri(stri=='É') = 'E';
stri(stri=='È') = 'E';
stri(stri=='å') = 'a';
stri(stri=='ä') = 'a';
stri(stri=='ö') = 'o';
stri(stri=='é') = 'e';
stri(stri=='è') = 'e';
stri(stri=='ü') = 'u';

logind = ...
  (stri==':')|...
  (stri=='\')|...
  (stri=='/')|...
  (stri=='?')|...
  (stri=='*')|...
  (stri=='"')|...
  (stri=='<')|...
  (stri=='>')|...
  (stri=='|');

stri(logind)='-';
takeind = (stri>64 & stri<91) | (stri=='-') | (stri=='_') | (stri=='^') | (stri>47 & stri<58) | (stri>96 & stri<123);
stri = stri(takeind);