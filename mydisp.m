function mydisp(stri)
%Internal display function
global DATA

stri = translation.dictionary(stri);

stri = char(textwrap({stri},60));

try
  if ~DATA.Silent
    disp(stri);
  end;
catch %#ok<CTCH>
 disp(stri);
end;