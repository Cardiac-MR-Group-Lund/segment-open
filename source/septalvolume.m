function [varargout] = septalvolume(varargin)
%Septal Volume Gui for Katarina and Pia illustrating the effect of translation.

if nargin==0
varargin={'init'};
end

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%---------------------------------------------------------
function init
%---------------------------------------------------------
global DATA SET NO

if isempty(SET(NO).EpiX)
  myfailed('No epicardium');
  return
end;

if ~isempty(find(strcmp(SET(NO).Point.Label,'RV insertion'))) && ~isempty(find(strcmp(SET(NO).Point.Label,'RV insertion Inferior')))
   myfailed('Found inconsistence in insertion point name. There is a mix of RV insertion points with RV insertion Inferior/Anterior points.');
     return
end
if ~isempty(find(strcmp(SET(NO).Point.Label,'RV insertion'))) && ~isempty(find(strcmp(SET(NO).Point.Label,'RV insertion Anterior')))
   myfailed('Found inconsistence in insertion point name. There is a mix of RV insertion points with RV insertion Inferior/Anterior points.');
  return
end

if isempty(find(strcmp(SET(NO).Point.Label,'RV insertion'))) && isempty(find(strcmp(SET(NO).Point.Label,'RV insertion Inferior')))
   myfailed('No RV insertion point');
     return
end
if isempty(find(strcmp(SET(NO).Point.Label,'RV insertion'))) && isempty(find(strcmp(SET(NO).Point.Label,'RV insertion Anterior')))
   myfailed('No RV insertion point');
  return
end

% if find(strcmp(SET(NO).Point.Label,'RV insertion'))
%   mymsgbox('RV insertion points found. You are using older method to calculate septum motion');
% end
% 
% if isempty(find(strcmp(SET(NO).Point.Label,'RV insertion'))) 
%    myfailed('No RV insertion point');
%   return
% end

%reboot if started
if isopengui(fullfile('septalvolume.fig'))
  close_Callback
end

DATA.GUI.septalvolume = mygui(fullfile('septalvolume.fig'));
gui = DATA.GUI.septalvolume;
gui.no = NO;


[septalvolume,rvlateralvolume,lvlateralvolume,numused,numignored,slicestodo,graphics] = utilityprivate('calcseptumvolume',false);%,true);

[septalvolumecorr,rvlateralvolumecorr,lvlateralvolumecorr,numused,numignored,slicestodo,graphicscorr]=utilityprivate('calcseptumvolumegui',false);

%store the graphical values to be used for each plot.
gui.graphics = graphics;
gui.slicestodo = slicestodo;
gui.graphicscorr=graphicscorr;
gui.volumedata=[septalvolume, rvlateralvolume, lvlateralvolume; septalvolumecorr, rvlateralvolumecorr, lvlateralvolumecorr];
gui.noslice=[numused, numignored];


set(gui.handles.septalvolumewithoutcorrectiontext,'String', sprintf('%s',num2str(septalvolume)));
set(gui.handles.septalvolumewithcorrectiontext,'String', sprintf('%s',num2str(septalvolumecorr)));
set(gui.handles.lateralvolumewithoutcorrectiontext,'String', sprintf('%s',num2str(lvlateralvolume)));
set(gui.handles.lateralvolumewithcorrectiontext,'String', sprintf('%s',num2str(lvlateralvolumecorr)));
set(gui.handles.slicetext1,'String', sprintf('%s',num2str(lvlateralvolumecorr)));

if ~isnan(rvlateralvolume)
set(gui.handles.rvvolumewithoutcorrectiontext,'String', sprintf('%s',num2str(rvlateralvolume)));
set(gui.handles.rvvolumewithcorrectiontext,'String', sprintf('%s',num2str(rvlateralvolumecorr)));
else
  set(gui.handles.rvvolumewithoutcorrectiontext,'String', sprintf('-'));
set(gui.handles.rvvolumewithcorrectiontext,'String', sprintf('-'));
end
%Set the slider so that it matches the middle slice in the slider callback
%all images are updated aswell
gui.sliceind = 1;%slicestodo(1);%SET(gui.no).CurrentSlice; 

if length(slicestodo)==1%SET(gui.no).ZSize==1
  set(gui.handles.sliceslider,'visible','off')
  set(gui.handles.slicetext1,'String', sprintf('%s',num2str(slicetodo)));
else
  set(gui.handles.sliceslider,'Max',length(slicestodo))
  set(gui.handles.sliceslider,'Value',1+length(slicestodo)-gui.sliceind)
  set(gui.handles.sliceslider,'Min',1)  
  set(gui.handles.sliceslider,'Sliderstep',[1/(length(slicestodo)-1),0.1])
  set(gui.handles.slicetext1,'String', sprintf('%s',num2str(slicestodo(gui.sliceind))));
end

%plot1
slider_Callback
%update is called in slider callback therefore we are done after this
exportseptumvolumegui

%-----------------------
function update%#ok<DEFNU>
%-----------------------
global DATA
gui = DATA.GUI.septalvolume;
%This function updates images and results

im = getanatomicimage;
drawanatomicimage(im)

set(gui.handles.slicetext1,'String', sprintf('%s',num2str(gui.slicestodo(gui.sliceind))));
%Show plots
plot1
plot2
plot3
plot4
plot5
plot6
if isempty(gui.graphics.rvepixed)
  set(gui.handles.rvlateralaxes,'visible','off')
  set(gui.handles.rvlateralcorraxes,'visible','off')
else
  plot7
  plot8
end

%---------------------
function plot1
%---------------------
global DATA SET
gui = DATA.GUI.septalvolume;
%if ~isempty(gui.graphics(gui.sliceind).rvepiyed)

plot(gui.handles.plotaxes,gui.graphics(gui.sliceind).lvepiyedcen,gui.graphics(gui.sliceind).lvepixedcen, '+b',gui.graphics(gui.sliceind).lvepiyescen,...
  gui.graphics(gui.sliceind).lvepixescen,'+r', gui.graphics(gui.sliceind).lvepiyed, gui.graphics(gui.sliceind).lvepixed,'b',...
  gui.graphics(gui.sliceind).lvepiyes,gui.graphics(gui.sliceind).lvepixes,'r',...
  gui.graphics(gui.sliceind).rvepiyed,gui.graphics(gui.sliceind).rvepixed,'g',...
  gui.graphics(gui.sliceind).rvepiyes,gui.graphics(gui.sliceind).rvepixes,'m',...
  gui.graphics(gui.sliceind).antedy,gui.graphics(gui.sliceind).antedx,'*b',...
  gui.graphics(gui.sliceind).infedy, gui.graphics(gui.sliceind).infedx,'*b', gui.graphics(gui.sliceind).antesy,gui.graphics(gui.sliceind).antesx, 'or',...
  gui.graphics(gui.sliceind).infesy,gui.graphics(gui.sliceind).infesx,'or')

set(gui.handles.plotaxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [],'Ydir','reverse')
%'xlim',[0 140],'ylim',[0 140]
axis normal

legendaxesplot=plot(gui.handles.legendaxes,gui.graphics(gui.sliceind).lvepiyedcen,gui.graphics(gui.sliceind).lvepixedcen, '+b',gui.graphics(gui.sliceind).lvepiyescen,...
  gui.graphics(gui.sliceind).lvepixescen,'+r', gui.graphics(gui.sliceind).lvepiyed, gui.graphics(gui.sliceind).lvepixed,'b',...
  gui.graphics(gui.sliceind).lvepiyes,gui.graphics(gui.sliceind).lvepixes,'r',...
  gui.graphics(gui.sliceind).rvepiyed,gui.graphics(gui.sliceind).rvepixed,'g',...
  gui.graphics(gui.sliceind).rvepiyes,gui.graphics(gui.sliceind).rvepixes,'m',...
  gui.graphics(gui.sliceind).antedy,gui.graphics(gui.sliceind).antedx,'*b',...
  gui.graphics(gui.sliceind).infesy,gui.graphics(gui.sliceind).infesx,'or');

set(gui.handles.legendaxes, 'Visible', 'off');
set(legendaxesplot,'Visible','off');

legend(gui.handles.legendaxes,{'LV ED cent point','LV ES cent point','LV ED','LV ES','RV ED','RV ES', 'RV insertion ES', 'RV insertion ED'},'Fontsize',8,'Location','northwest','Orientation','vertical');%,'NumColumns',3)%'Fontsize',7,


%---------------------
function plot2
%---------------------
global DATA SET
gui = DATA.GUI.septalvolume;

%plot(gui.handles.septalaxes,gui.graphics(gui.sliceind).areaxsep,gui.graphics(gui.sliceind).areaysep)
plot(gui.handles.septalaxes,gui.graphics(gui.sliceind).areaysep,gui.graphics(gui.sliceind).areaxsep)
set(gui.handles.septalaxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [],'Ydir','reverse')
xb = get(gui.handles.plotaxes,'xlim');
yb = get(gui.handles.plotaxes,'ylim');
set(gui.handles.septalaxes,'xlim',xb,'ylim',yb)

%---------------------
function plot3
%---------------------
global DATA SET
gui = DATA.GUI.septalvolume;
% plot(gui.handles.lateralaxes,gui.graphics(gui.sliceind).areaxlat,gui.graphics(gui.sliceind).areaylat)
plot(gui.handles.lateralaxes,gui.graphics(gui.sliceind).areaylat,gui.graphics(gui.sliceind).areaxlat)
set(gui.handles.lateralaxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [],'Ydir','reverse')
xb = get(gui.handles.plotaxes,'xlim');
yb = get(gui.handles.plotaxes,'ylim');
set(gui.handles.lateralaxes,'xlim',xb,'ylim',yb)


%---------------------
function plot4
%---------------------
global DATA SET
gui = DATA.GUI.septalvolume;
%   gui.graphics(gui.sliceind).rvepiyed,gui.graphics(gui.sliceind).rvepixed,'g',...
%    gui.graphics(gui.sliceind).rvepiyes,gui.graphics(gui.sliceind).rvepixes,'m',...

plot(gui.handles.plotcorraxes,gui.graphicscorr(gui.sliceind).lvepiyedcen,gui.graphicscorr(gui.sliceind).lvepixedcen, '+b',gui.graphicscorr(gui.sliceind).lvepiyescen,...
  gui.graphicscorr(gui.sliceind).lvepixescen,'+r', gui.graphicscorr(gui.sliceind).lvepiyed, gui.graphicscorr(gui.sliceind).lvepixed,'b',...
  gui.graphicscorr(gui.sliceind).rvepiyed,gui.graphicscorr(gui.sliceind).rvepixed,'g',...
   gui.graphicscorr(gui.sliceind).rvepiyes,gui.graphicscorr(gui.sliceind).rvepixes,'m',...
  gui.graphicscorr(gui.sliceind).lvepiyes,gui.graphicscorr(gui.sliceind).lvepixes,'r',gui.graphicscorr(gui.sliceind).antedy,gui.graphicscorr(gui.sliceind).antedx,'*b',...
  gui.graphicscorr(gui.sliceind).infedy, gui.graphicscorr(gui.sliceind).infedx,'*b', gui.graphicscorr(gui.sliceind).antesy,gui.graphicscorr(gui.sliceind).antesx, 'or',...
  gui.graphicscorr(gui.sliceind).infesy,gui.graphicscorr(gui.sliceind).infesx,'or')
set(gui.handles.plotcorraxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [],'Ydir','reverse')
xb = get(gui.handles.plotaxes,'xlim');
yb = get(gui.handles.plotaxes,'ylim');
set(gui.handles.plotcorraxes,'xlim',xb,'ylim',yb)


%---------------------
function plot5
%---------------------
global DATA SET 
gui = DATA.GUI.septalvolume;

%plot(gui.handles.septalcorraxes,gui.graphicscorr(gui.sliceind).areaxsep,gui.graphicscorr(gui.sliceind).areaysep)
plot(gui.handles.septalcorraxes,gui.graphicscorr(gui.sliceind).areaysep, gui.graphicscorr(gui.sliceind).areaxsep)
set(gui.handles.septalcorraxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [],'Ydir','reverse')
xb = get(gui.handles.plotaxes,'xlim');
yb = get(gui.handles.plotaxes,'ylim');
set(gui.handles.septalcorraxes,'xlim',xb,'ylim',yb)


function plot6
%---------------------
global DATA SET 
gui = DATA.GUI.septalvolume;

plot(gui.handles.lateralcorraxes,gui.graphicscorr(gui.sliceind).areaxlat,gui.graphicscorr(gui.sliceind).areaylat)
plot(gui.handles.lateralcorraxes,gui.graphicscorr(gui.sliceind).areaylat,gui.graphicscorr(gui.sliceind).areaxlat)
set(gui.handles.lateralcorraxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [],'Ydir','reverse')
xb = get(gui.handles.plotaxes,'xlim');
yb = get(gui.handles.plotaxes,'ylim');
set(gui.handles.lateralcorraxes,'xlim',xb,'ylim',yb)

%---------------------
function plot7
%---------------------
global DATA SET
gui = DATA.GUI.septalvolume;

% plot(gui.handles.lateralaxes,gui.graphics(gui.sliceind).areaxlat,gui.graphics(gui.sliceind).areaylat)
plot(gui.handles.rvlateralaxes,gui.graphics(gui.sliceind).areayrvlat,gui.graphics(gui.sliceind).areaxrvlat)
set(gui.handles.rvlateralaxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [],'Ydir','reverse')
xb = get(gui.handles.plotaxes,'xlim');
yb = get(gui.handles.plotaxes,'ylim');
set(gui.handles.rvlateralaxes,'xlim',xb,'ylim',yb)

%---------------------
function plot8
%---------------------
global DATA SET
gui = DATA.GUI.septalvolume;

% plot(gui.handles.lateralaxes,gui.graphics(gui.sliceind).areaxlat,gui.graphics(gui.sliceind).areaylat)
plot(gui.handles.rvlateralcorraxes,gui.graphicscorr(gui.sliceind).areayrvlat,gui.graphicscorr(gui.sliceind).areaxrvlat)
set(gui.handles.rvlateralcorraxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [],'Ydir','reverse')
xb = get(gui.handles.plotaxes,'xlim');
yb = get(gui.handles.plotaxes,'ylim');
set(gui.handles.rvlateralcorraxes,'xlim',xb,'ylim',yb)

%-----------------------
function slider_Callback %#ok<DEFNU>
%-----------------------
%callback from slice slider to change image slice in display

global DATA %SET
gui = DATA.GUI.septalvolume;

gui.sliceind = (1+length(gui.slicestodo)-max(1,round(get(gui.handles.sliceslider,'Value'))));
update

%------------------------------------------------
function im = getanatomicimage
%------------------------------------------------
global DATA SET
gui = DATA.GUI.septalvolume;
imed = SET(gui.no).IM(:,:,SET(gui.no).EDT,gui.slicestodo(gui.sliceind));
imes = SET(gui.no).IM(:,:,SET(gui.no).EST,gui.slicestodo(gui.sliceind));
im = [imed,imes];

%------------------------------------------------
function drawanatomicimage(im)
%------------------------------------------------
global DATA SET
gui = DATA.GUI.septalvolume;
cmap = gray(256);
im = imresize(im,2,'bicubic');
c =SET(gui.no).IntensityMapping.Contrast;
b = SET(gui.no).IntensityMapping.Brightness;
rim = calcfunctions('remapuint8viewim',im,gui.no,cmap(:,1),c,b);
gim = calcfunctions('remapuint8viewim',im,gui.no,cmap(:,2),c,b);
bim = calcfunctions('remapuint8viewim',im,gui.no,cmap(:,3),c,b);
im = cat(3,rim, gim, bim);
gui.handles.imagehandle = image(im,'parent',gui.handles.imageaxes);
set(gui.handles.imageaxes,'dataaspectratio',...
  [1/SET(gui.no).ResolutionY ...
  1/SET(gui.no).ResolutionX 1],'xtick', [], 'ytick', [])


%--------------------------
function exportseptumvolumegui %#ok<DEFNU>
%--------------------------
%Export septum volumes for current file 
global DATA SET NO
gui = DATA.GUI.septalvolume;

%Create output matrix
%Filename, ED,ES,SV,SEPVOL,SEP%
outdata = cell(3,16); 
outdata{1, 1} = 'Filename';
outdata{1, 2} = 'Status';
outdata{1, 3} = 'LVED [ml]';
outdata{1, 4} = 'LVES [ml]';
outdata{1, 5} = 'LVSV [ml]';
outdata{1, 6} = 'RVED [ml]';
outdata{1, 7} = 'RVES [ml]';
outdata{1, 8} = 'RVSV [ml]';
outdata{1, 9} = 'SepVol [ml]';
outdata{1,10} = 'SepVol%LVSV [%]';
outdata{1,11} = 'RVlatVol [ml]';
outdata{1,12} = 'RVlatVol%RVSV [%]';
outdata{1,13} = 'LVlatVol [ml]';
outdata{1,14} = 'LVlatVol%LVSV [%]';
outdata{1,15} = 'Used slices';
outdata{1,16} = 'Ignored slices';

%find cine short-axis
[normalno] = findfunctions('findcineshortaxisno');
NO = normalno;
if NO==0
  NO=1;
end;
if isnan(NO)
  NO = 1;
end;
if isempty(NO)
  NO = 1;
end;

%Set filename
outdata{2,1} = SET(NO).PatientInfo.Name;
doplot = false;

%[sepvolume,rvlateralvolume,lvlateralvolume,numused,numignored] = calcseptumvolume(doplot);

if isnan(gui.volumedata(1,2)) %(rvlateralvolume)
  dorv = false;
else
  dorv = true;
end;

outdata{2, 2} = 'Uncorrected';
outdata{2, 3} = SET(NO).EDV;
outdata{2, 4} = SET(NO).ESV;
outdata{2, 5} = SET(NO).SV;

outdata{2, 6} = SET(NO).RVEDV;
outdata{2, 7} = SET(NO).RVESV;
outdata{2, 8} = SET(NO).RVSV;

outdata{2, 9} = gui.volumedata(1,1);        %sepvolume;
outdata{2, 10} = 100*gui.volumedata(1,1)/SET(NO).SV; %Bugfixed: was ESV!!!

if dorv
  outdata{2,11} = gui.volumedata(1,2);  %rvlateralvolume;
  outdata{2,12} = 100*gui.volumedata(1,2)/SET(NO).RVSV;
else
  outdata{2,11} = NaN;
  outdata{2,12} = NaN;
end;

outdata{2,13} = gui.volumedata(1,3);     %lvlateralvolume;
outdata{2,14} = 100*gui.volumedata(1,3)/SET(NO).SV;

outdata{2,15} = gui.noslice(1); %numused;
outdata{2,16} = gui.noslice(2); %numignored;

%Set filename
outdata{2,1} = SET(NO).PatientInfo.Name;
doplot = false;

% [sepvolume,rvlateralvolume,lvlateralvolume,numused,numignored] = calcseptumvolumegui(doplot);
% 
% if isnan(rvlateralvolume)
%   dorv = false;
% else
%   dorv = true;
% end;
outdata{3, 1} = outdata{2,1};
outdata{3, 2} = 'Corrected';
outdata{3, 3} = SET(NO).EDV;
outdata{3, 4} = SET(NO).ESV;
outdata{3, 5} = SET(NO).SV;

outdata{3, 6} = SET(NO).RVEDV;
outdata{3, 7} = SET(NO).RVESV;
outdata{3, 8} = SET(NO).RVSV;

outdata{3, 9} = gui.volumedata(2,1);  %corrected septalvolume
outdata{3, 10} = 100*gui.volumedata(2,1)/SET(NO).SV; %Bugfixed: was ESV!!!

if dorv
  outdata{3,11} = gui.volumedata(2,2);  % corrected rvlateralvolume
  outdata{3,12} = 100*gui.volumedata(2,2)/SET(NO).RVSV;
else
  outdata{3,11} = NaN;
  outdata{3,12} = NaN;
end;

outdata{3,13} = gui.volumedata(2,3);  %corrected lvlateralvolume;
outdata{3,14} = 100*gui.volumedata(2,3)/SET(NO).SV;

outdata{3,15} = gui.noslice(1);   %numused;
outdata{3,16} = gui.noslice(2);  %numignored;

%--- Output to a string
segment('cell2clipboard',outdata);



%------------------------------------------------
function close_Callback
%------------------------------------------------
global DATA 
gui = DATA.GUI.septalvolume;
DATA.GUI.septalvolume = close(gui);
