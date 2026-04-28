function mydisp(stri, logfid, dowrap) %skiptranslation
%Internal display function, masks patient info

arguments
  stri
  logfid = '';
  dowrap = true
end

stri = tools('maskpatientstrings',stri);
if dowrap
  stri = char(textwrap({stri},60));
end

if isempty(logfid)
  stri = reshape(stri',1,numel(stri));
  if strlength(stri)>3 && strncmp(stri,num2str(year(datetime('now'))),4)
    % this string already has the time stamp
    fprintf(stri)
  else
    logdisp(stri)
  end
else
  fprintf(logfid, stri);
end
