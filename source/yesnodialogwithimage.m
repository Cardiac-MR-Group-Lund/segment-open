function usersanswer = yesnodialogwithimage(imagedata,txtmessage)
% function to show a user dialog that shows an image and two buttons with 
% options "YES,continue" and "NO, return"
global DATA
% get current screensize in pixels
set(0,'units','pixels')
screensizeinpix = get(0,'screensize');
usersanswer = false;
% calculate figure position
figwidth = 300;
figheight = 500;
figstartx = floor(screensizeinpix(3)-figwidth)/2;
figstarty = floor(screensizeinpix(4)-figheight)/2;
% create an invisible figure
fig = figure(...
      ...'Visible', 'off',...
       'Name',dprintf('Warning'), ...
       'MenuBar','none',...
       'Resize', 'off',...
       'keypressfcn',@keypressed,...
       'Position',[figstartx figstarty figwidth figheight],...
       'Color',DATA.GUISettings.BackgroundColor,...
       'NumberTitle','off'...
       );
     
% load in image ad crop
offset = 20;
if ~exist('imagedata','var')|| isempty(imagedata)
  if isdeployed
    tmpimg = imread('drawing_guidance.png');
  else
    tmpimg = imread(['+straintagging' filesep 'drawing_guidance.png']);
  end
  imwidth = floor(size(tmpimg,2)/4);
  imheight = size(tmpimg,1)- 50;
  img = imcrop(tmpimg,[2*imwidth 0 imwidth imheight]);
else
  imwidth = size(imagedata,2);
  imheight = size(imagedata,1);
  img = imagedata;
end

% add axes to the figure and show image
ax = axes(fig);
ax.Units = 'pixels';
axeswidth = figwidth-2*offset;
axesheight = floor(axeswidth*imheight/imwidth);
ax.Position = [offset 135 axeswidth axesheight];
imagesc(ax,img)
axis off;

% add text message
if ~exist('txtmessage','var')|| isempty(txtmessage)
  txtmessage = dprintf('Warning');
end
txtposition = [offset 70 axeswidth 60];
msgtextcontrol = uicontrol(fig,...
                          'Style','text',...
                          'Units','pixels',...
                          'Position',txtposition,...
                          'keypressfcn',@yeskeypressed,... 
                          'HorizontalAlignment', 'left',...
                          'BackgroundColor',DATA.GUISettings.BackgroundColor,...
                          'ForegroundColor',DATA.GUISettings.ForegroundColor,...
                          'String',txtmessage); %#ok<NASGU>
% add buttons
buttonheight = 30;
yesbuttonposition = [offset offset axeswidth/2-10 buttonheight];
yesbuttoncontrol = uicontrol(fig,...
                          'Style','pushbutton',...
                          'Units','pixels',...
                          'Position',yesbuttonposition,...
                          'Callback',@selection,...
                          'BackgroundColor',DATA.GUISettings.BackgroundColor,...
                          'ForegroundColor',DATA.GUISettings.ForegroundColor,...
                          'String',dprintf('YES, continue'),...
                          'UserData','yes'); %#ok<NASGU>
                        
nobuttonposition = [offset+axeswidth/2+10 offset axeswidth/2-10 buttonheight];
nobuttoncontrol = uicontrol(fig,...
                          'Style','pushbutton',...
                          'Units','pixels',...
                          'Position',nobuttonposition,...
                          'Callback',@selection,...                          
                          'BackgroundColor',DATA.GUISettings.BackgroundColor,...
                          'ForegroundColor',DATA.GUISettings.ForegroundColor,...
                          'String',dprintf('NO, return'),...
                          'UserData','no'); %#ok<NASGU>
                        
uiwait(fig)  
    function selection(src,~)            
      str = src.UserData;
      if contains(lower(str),'yes')
        usersanswer = true;
      else
        usersanswer = false;
      end
      uiresume(fig)
      close(fig)            
    end 
  
    %----------------------------------
    function keypressed(fignum,evnt) %#ok<INUSL>
    %----------------------------------
      key = getkey(evnt);
      switch key
        case {'y','Y'}
          usersanswer = true;
          uiresume(fig)
          close(fig)
        case {'n','N'}
          usersanswer = true;
          uiresume(fig)
          close(fig)
        otherwise
          
      end
      
    end
end
