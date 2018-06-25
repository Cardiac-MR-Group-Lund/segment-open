function [x1,y1]=mypoly2cw(x2,y2)
%Finds clockwise orientation of segmentation 
%Written by Klas.
try
  k=convhull(x2,y2);
%Convhull returns counterclockwise sorted segmentation do flip to get
%clockwise orientation.
catch
  %convhull fails for some small accidental drawings return nan in this
  %case
  x1=nan(size(x2));
  y1=nan(size(y2));
  return
end
%Check if index order is mostly flipped by convhull that means curve has
%mostly clockwise orientation
if sum(diff(k(2:end))<0)>length(k)/2
  x1=x2;
  y1=y2;
else
  x1=fliplr(x2);
  y1=fliplr(y2);
end

% %first direction is 
% dir(1)=det([1,1,1;x2(k(end-1)),x2(k(1)),x2(k(2));y2(k(end-1)),y2(k(1)),y2(k(2))]');
% for i=2:length(k)-1
%   dir(i)=det([1,1,1;x2(k(i-1)),x2(k(i)),x2(k(i+1));y2(k(i-1)),y2(k(i)),y2(k(i+1))]');
% end

%For direction at ind x is mapped to first ind in k smaller than or equal
%to x
%closed contour so last index equal to first
%x1=cat(2,x2(k(1):k(end-1)-1),x2(k(end-1):end),x2(2:k(1)));
%y1=cat(2,y2(k(1):k(end-1)-1),y2(k(end-1):end),y2(2:k(1)));


%Next level scheme which unties all knots aswell. However might not be
%desirable

% %pad with end and first
% x_tmp=[x2(end-1);x2(k);x2(2)];
% y_tmp=[y2(end-1);y2(k);y2(2)];
% 
% calculate determinant of orientation matrix
% O_tmp=ones(3);
% for i=2:N
%   O(2,:)=x_tmp(i-1:i+1);
%   O(3,:)=y_tmp(i-1:i+1);
%   direction(i)=det(O);
% end
% 
% get leftmost point
% [~,indleft]=min(x2);
% [~,indright]=max(x2);
% 
% these must be located on the convexhull
% xleft=x2(k==indleft);
% xright=x2(k==indleft);
% [~,minind]=[indleft,indright]; 
% 
% 
% right clockwise direction and clockwise contour i.e left to right must be
% top ergo bottom right to left
% if minind==1
%   topx=x2(indleft:indright-1);
%   topy=y2(indleft:indright-1);
%   
%   bottomx=x2([indright:end,2:indleft-1]);
%   bottomy=y2([indright:end,2:indleft-1]);
% else
%   topx=x2([indleft:end,2:indright-1]);
%   topy=y2([indleft:end,2:indright-1]);
%   
%   topx=x2(indright:indleft-1);
%   topy=y2(indright:indleft-1);
% end



% %pad with end and first
% x_tmp=[x2(end-1);x2(k);x2(2)];
% y_tmp=[y2(end-1);y2(k);y2(2)];
% 
% calculate determinant of orientation matrix
% O_tmp=ones(3);
% for i=2:N
%   O(2,:)=x_tmp(i-1:i+1);
%   O(3,:)=y_tmp(i-1:i+1);
%   direction(i)=det(O);
% end


