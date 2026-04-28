function addcommenttostack

%Add a comment to current image stack

%Nils Lundahl

%Renamed and move to main segment code by Einar
%Updated to myinputstruct by Fanny Månefjord 2024

global DATA SET NO %#ok<GVMIS> 

if ~DATA.DataLoaded
  myfailed('No data loaded.');
  return;
end

comment = struct('Username','Unknown','Time','','Text','');

%Extract environemnt variable containing username of windows machines
if ispc
  comment.Username = getenv('USERNAME');
else
  comment.Username = 'Unkown'; %Added Einar
end

%txt1 = myinputdlg(comment.Username,'Enter comment',[10 80]);
s = [];
s(1).Field = 'username';
s(1).Label = dprintf('Username');
s(1).Default = comment.Username;
s(2).Field = 'text';
s(2).Label = dprintf('Comment');
s(2).Default = "";

[s,ok] = myinputstruct(s, 'Add comment', 30);
if ok
  txt = {s.text};
  comment.Username = s.username;
else
  return; %when user presses cancel
end

%Add comment
if ~isempty(txt{1})
  stri = '';
  for li = 1:size(txt{1},1)
    stri = [stri strtrim(txt{1}(li,:))]; %#ok<AGROW> 
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

