function reportradvel
%--------------------
%GUI for radial velocity.
global DATA SET NO

if SET(NO).TSize < 2
  myfailed('Data needs to be time resolved.',DATA.GUI.Segment);
  return;
end

if isempty(SET(NO).EndoX)
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end

tempnos=NO;
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
  return;
end

ind = findfunctions('findslicewithendo',NO);
pos = find(ind);
nslices = length(pos);
radvel = calcfunctions('calcradialvelocity',NO);

if nslices<1
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end

fig = mygui('reportradvel.fig');

% Generate a structure of handles to pass to callbacks, and store it.
handles = fig.handles;
handles.fig = handles.figure1;

%Calculate over sectors
[meanradvel,meanx,meany,sectors] = calcfunctions('findmeaninsector',...
  'endo',radvel,find(ind),6);

%temp = meanradvel(:,:,2:end);
maxv = max(abs(meanradvel(:)));
konst = 32/maxv;

axes(handles.coloraxes); 
h = imagesc([-maxv maxv]);
colorbar;
set(h,'visible','off');
axis off;

axes(handles.imageaxes); 
tf = SET(NO).CurrentTimeFrame;

%Calculate "current slice"
cs = pos(round(length(pos)/2));
image(65+64*SET(NO).IM(:,:,SET(NO).CurrentTimeFrame,cs));
colormap([jet;gray]);
axis image
axis off;

hold on
h = line('XData',SET(NO).EndoY(:,tf,cs),'YData',SET(NO).EndoX(:,tf,cs),'Color','r','LineStyle','-');
set(h,'linewidth',3);
for loop=1:6
  line('XData',[meany(tf,cs) SET(NO).EndoY(sectors(loop,1,tf),tf,cs)],...
    'YData',[meanx(tf,cs) SET(NO).EndoX(sectors(loop,1,tf),tf,cs)],...
    'Color','r','LineStyle','-');
  text(...
    0.5*(SET(NO).EndoY(sectors(loop,1,tf),tf,cs)+SET(NO).EndoY(sectors(loop+1,1,tf),tf,cs)),...
    0.5*(SET(NO).EndoX(sectors(loop,1,tf),tf,cs)+SET(NO).EndoX(sectors(loop+1,1,tf),tf,cs)),...
    sprintf('%d',loop));
end
hold off;
axis off;
axis image;

axes(handles.axes1); 
helpimageposneg(meanradvel,1,konst);

axes(handles.axes2); 
helpimageposneg(meanradvel,2,konst);

axes(handles.axes3); 
helpimageposneg(meanradvel,3,konst);

axes(handles.axes4);
helpimageposneg(meanradvel,4,konst);

axes(handles.axes5); 
helpimageposneg(meanradvel,5,konst);

axes(handles.axes6); 
helpimageposneg(meanradvel,6,konst);

%---------------------------------------
function helpimageposneg(im,sector,konst)
%---------------------------------------
%Helper function to display image.
global SET NO DATA

im = im(sector,:,:);
im = squeeze(im);
im = im(:,2:end);
h = image(32+konst*im);
set(h,'xdata',SET(NO).TIncr*[0.5 SET(NO).TSize+0.5]);
set(gca,'yticklabel','','xlim',[0 (SET(NO).TSize-1)*SET(NO).TIncr]);
set(gca,'XColor',DATA.GUISettings.ForegroundColor,'YColor',DATA.GUISettings.ForegroundColor);
xlabel(dprintf('Time'),'Color',DATA.GUISettings.ForegroundColor);
title(dprintf('Sector %d',sector),'Color',DATA.GUISettings.ForegroundColor);
