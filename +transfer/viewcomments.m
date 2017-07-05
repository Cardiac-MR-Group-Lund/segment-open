function viewcomments
% View all comments on active image stacks

%Nils Lundahl

%Modified with fixcommentstri by Einar Heiberg

global SET

if ~isfield(SET,'Comment') || isempty([SET.Comment])
  myfailed('No comments available')
  return
end

stri = '';
%Loop over all image stacks
for no = 1:numel(SET)
  nostri = '';
  if ~isempty(SET(no).Comment)
    comments = SET(no).Comment;
    nostri = dprintf('Comments on image stack %d:\n\n',no);
    for cmt = comments
      txt = sprintf('%s @ %s\n%s\n\n',cmt.Username,cmt.Time,fixcommentstri(cmt.Text));
      nostri = [nostri txt];
    end
  end
  stri = [stri nostri];
end

msgbox(stri,'Comments')

%------------------------------------------
function commentstri = fixcommentstri(stri)
%------------------------------------------
%Ensures that commentstri is a vector of chars.

if size(stri,1)>1
  commentstri = '';
  for loop = 1:size(stri,1)
    commentstri = [commentstri sprintf('\n') stri(loop,:)]; %#ok<AGROW>
  end;
else
  commentstri = stri;
end;