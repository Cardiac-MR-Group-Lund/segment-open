function addcomment
%Add a comment to current image stack

%Nils Lundahl
global DATA SET NO

if ~DATA.DataLoaded
  myfailed('No data loaded.');
  return;
end;

comment = struct('Username','Unknown','Time','','Text','');

%Extract environemnt variable containing username of windows machines
if ispc
  comment.Username = getenv('USERNAME');
else
  comment.Username = 'Unkown'; %Added Einar
end

txt = inputdlg(comment.Username,'Enter comment',[10 80]);

%Add it
if ~isempty(txt)
  stri = '';
  for li = 1:size(txt{1},1)
    stri = [stri strtrim(txt{1}(li,:)) '\n'];
  end
  comment.Text = sprintf(stri);
  comment.Time = datestr(now,31);
  if isfield(SET,'Comment') && ~isempty(SET(NO).Comment)
    SET(NO).Comment = [SET(NO).Comment comment];
  else
    SET(NO).Comment = comment;
  end
else
  myfailed('Empty comment. Ignored.');
end

