function myrequestfocus(fig)
%Take away focus from current object (typically slider) and set to figure.
%This is to avoid the problem of sliders that "steal" arrow keys after one
%have used them. 
%
%We need to check if this problem exists in newer Matlabs.

%Einar Heiberg

%For now it does not do anything in newer Matlab. 
if strcmp(version('-release'),'2022a')
  warning('off','MATLAB:ui:javaframe:PropertyToBeRemoved')
  jFig = get(fig,'JavaFrame'); %#ok<JAVFM>
  jFig.requestFocus;
else
  logdisp('Ignored myrequest focus call.')
end
