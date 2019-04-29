function varargout = calcfunctions(varargin)
% CALCFUNCTIONS
% Functions for doing calculations

% Moved out from segment_main by Nisse Lundahl

%Invoke subfunction
macro_helper(varargin{:}); %future macro recording use
if (nargout)
  [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
else
  feval(varargin{:}); % FEVAL switchyard
end

%------------------------------------------
function [edmax, MaxdiameterPoint, Zslice_no]  = maxsaxdiameter(no,tf,type)
%------------------------------------------
%Get maximum endocardial diameter in short-axis cine stack
% Input:    no - no of stack, 
%            tf - time frame, 
%            type - 'RV', 'LV'
% Output:   edmax - max LV/RV diameter
%           MaxdiameterPoint - points with the largest distance in form [x,y]
%           Zslice_no - slice number with largest LV/RV diameter
% 
global SET
xres = SET(no).ResolutionX;
yres = SET(no).ResolutionY;

switch type
    case 'LV'
        edxall = squeeze(SET(no).EndoX(:,tf,:));
        edyall = squeeze(SET(no).EndoY(:,tf,:));
        
        edmax = 0;
        edmax_temp = 0;
        MaxdiameterPoint = [];
        Zslice_no=[];
        
        for z = 1:SET(no).ZSize
            edx = edxall(:,z);
            if ~isnan(edx(1))
                edy = edyall(:,z);
                edxmat = repmat(edx,1,length(edx));
                edymat = repmat(edy,1,length(edy));
                eddist = sqrt((xres*(edxmat'-edxmat)).^2+ ...
                    (yres*(edymat'-edymat)).^2);
                
                % edmax = max(edmax,max(eddist(:)));
                edmax_temp=max(eddist(:));
                if edmax_temp>edmax
                    edmax=edmax_temp;
                    [row, col] = find(ismember(eddist, edmax));
                    MaxdiameterPoint(1,:)=[edx(row(1)),edy(row(1))];
                    MaxdiameterPoint(2,:)=[edx(col(1)),edy(col(1))];
                    Zslice_no=z;  % #slice number
                end
            end
        end
        %% In case of RV max diameter, first we find center point of LV based on EPI contour
        % Then we finiding septum, and calculate the middle point of septum
        % The last part is to find the intersection between line (center LV and center septum) with RV endo contur.
        % This intersecton defines our maximal RV dimater.
    case 'RV'
        edmax = 0;
        edmax_temp = 0;
        MaxdiameterPoint=[];
        Zslice_no=[];
        
        edxall = squeeze(SET(no).RVEndoX(:,tf,:));
        edyall = squeeze(SET(no).RVEndoY(:,tf,:));
        epxLVall = squeeze(SET(no).EpiX(:,tf,:));
        epyLVall = squeeze(SET(no).EpiY(:,tf,:));
        
      
        %Center of Epi LV ROI for each slice
        xLVcen=mean(epxLVall);
        yLVcen=mean(epyLVall);
        
        % Convhull RV roi
        for z = 1:SET(no).ZSize
            edx = edxall(:,z);
            edxLV=epxLVall(:,z);
            if isnan(edxLV(1))
                continue
            end
            if ~isnan(edx(1))
                edy = edyall(:,z);
                k=convhull(edx,edy);
                
                xhull=edx(k);
                yhull=edy(k);
                
                len=diff(xhull).^2+diff(yhull).^2;
                
                %%maximum length line segment is made out of the septum points. An addition
                %to make this more robust is to pick out the n largest line segments and check the area of the concavity the maximum area concavity should be the LV.
                [~,ind]=max(len);
                a=[xhull(ind),yhull(ind)];
                b=[xhull(ind+1),yhull(ind+1)];
                
                %find closest points on RV contour
                [~,ind_a]=min((edx-a(1)).^2+(edy-a(2)).^2);
                [~,ind_b]=min((edx-b(1)).^2+(edy-b(2)).^2);
               
                %%number of poinst between A and B
                if (edx(ind_a) < edx(ind_a+1))&&(edx(ind_b) < edx(ind_b+1))
                    noABx=[edx(1:ind_a-1); edx(ind_b+1:end)];
                else
                    noABx=edx(ind_a+1:ind_b-1);
                end
               
                %% point which is in the middle between A and B
                ind_midAB=ceil(length(noABx)/2);
                % middle point
                if ind_a <= ind_midAB
                    midRVab=[edx(end-(ind_midAB-ind_a)), edy(end-(ind_midAB-ind_a))];
                else
                    midRVab=[edx(ind_a-ind_midAB), edy(ind_a-ind_midAB)];
                    
                end
            
                %%line between center LV and centre of septal RV
                x=[xLVcen(z), midRVab(1)];
                y=[yLVcen(z), midRVab(2)];
                p = polyfit(x,y,1);
          
                %%calculate line
                x_space=linspace(min(edx), max(edx),100);
                y_space=p(1)*x_space+p(2);
               
                %% Intersection between the line  LV - septal centre and RV segmentation
                [xi,yi] = polyxpoly(edx,edy,x_space,y_space);
                
                %% Max RV diameter
                edmax_temp=sqrt((xres*diff(xi)).^2+(yres*diff(yi)).^2);
                if edmax_temp>edmax
                    edmax=edmax_temp;
                    MaxdiameterPoint(1,:)=[xi(1),yi(1)];
                    MaxdiameterPoint(2,:)=[xi(2),yi(2)];
                    Zslice_no=z;  % #slice number
                end
            end
        end
end


%----------------------------------------
function [N,xc,yc,zc] = lsplanefit(x,y,z) %#ok<DEFNU>, used by makecut
%---------------------------------------
%This function calculates the plane least squares fit to the points given by the
%x,y,z coordinates. It uses the null space of a matrix formulation of the .

%An idea is to subtract the first point this will force the plane to lie on
%the first point.
%x=rand(1,10),y=rand(1,10),z=rand(1,10)

xc =  mean(x);
yc = mean(y);
zc = mean(z);

x= x - xc;
y= y - yc;
z= z - zc;

% %using leastsquares
% Sxx = x*x';
% Sxy = x*y';
% Syy = y*y';
% Sxz = x*z';
% Syz = y*z';
% 
% %using that there is no problem with fixating one param in the equation and
% %that the point cloud is zero centered i.e sum(x)=0 etc we get that it is
% %sufficient to solve the below system
% A = [Sxx,Sxy;Sxy,Syy];
% B = -[Sxz;Syz];
% 
% params = A\B;
% N = [params;1];

%using svd
A=[x;y;z];
[U,~,~] = svd(A);
N = cross(U(:,1),U(:,2));

% point = [0,0,0];
% normal = N';

%# a plane is a*x+b*y+c*z+d=0
%# [a,b,c] is the normal. Thus, we have to calculate
%# d and we're set
% d = -point*normal'; %'# dot product for less typing

% %# create x,y
% [xx,yy]=ndgrid(linspace(-1,1,10),linspace(-1,1,10));
% 
% %# calculate corresponding z
% zz = (-normal(1)*xx - normal(2)*yy - d)/normal(3);

%# plot the surface
% figure
% surf(xx,yy,zz)
% hold on 
% plot3(x',y',z','k*')

%----------------------
function calcvolume(no) %#ok<DEFNU>
%----------------------
%Calculate volume of segmentation and updates. Updates both
%lv and rv segmentation. Calls subfunctions to do the work.

calclvvolume(no);
calcrvvolume(no);

volume_helper(no); %Find peak ejection rate, and empty volumes etc.

%-------------------------------------------
function varargout = calclvvolume(no,docomp)
%-------------------------------------------
%Calculate LV volume. Docomp if to use longaxis motion, see below.
%Uses area*(thickness+slicedist)
%NOTE: the exported LVM do NOT include Papillary volume (PV)

global DATA SET

if nargin<2
  docomp=true;
elseif SET(no).Longaxis > 1%SET(NO).Longaxis > 1
  docomp = true;
end;

if SET(no).Rotated
  calclvvolumepolar(no);
  return;
end;

ind = (findfunctions('findslicewithendo',no))|(findfunctions('findslicewithepi',no)); %Accepts empty endo and epi if exist

if ~any(ind)
  SET(no).LVV = nan(1,SET(no).TSize);
  SET(no).EPV = SET(no).LVV;
  SET(no).LVM = SET(no).LVV;
  SET(no).PV = zeros(size(SET(no).LVV));
  SET(no).PFR = 0;
  SET(no).PER = 0;
  SET(no).PFRT = 1;
  SET(no).PERT = 1;
  SET(no).ESV = 0;
  SET(no).EDV = 0;
  SET(no).EF = 0;
  SET(no).SV = 0;
  if nargout>0
    varargout = cell(1,1);
    varargout{1} = [];
  end;
  if nargout>1
    varargout{2} = [];
  end;  
  return;
end;

%Find what slices to do
if SET(no).ZSize>1
  pos = find(ind);
else
  pos = 1;
end;
LVVall = zeros(length(pos),SET(no).TSize);
EPVall = zeros(length(pos),SET(no).TSize);

% if isempty(SET(no).EndoX)
%   %This far then create
%   SET(no).EndoX = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
%   SET(no).EndoY = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);  
% end;
% 
% if isempty(SET(no).EpiX)
%   %This far then create
%   SET(no).EpiX = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);
%   SET(no).EpiY = nan(DATA.NumPoints,SET(no).TSize,SET(no).ZSize);  
% end;

%Loop over all segmented slices
for sloop=1:length(pos)
  for tloop=1:SET(no).TSize
    if ~isempty(SET(no).EndoX)&&~isnan(SET(no).EndoX(1,tloop,pos(sloop)))
      A = stablepolyarea(...
        SET(no).ResolutionY*SET(no).EndoY(1:end-1,tloop,pos(sloop)),...
        SET(no).ResolutionX*SET(no).EndoX(1:end-1,tloop,pos(sloop)));
      LVVall(sloop,tloop)=A*(SET(no).SliceThickness+SET(no).SliceGap)/1000; %to cm3
    end;
    if ~isempty(SET(no).EpiX)&&~isnan(SET(no).EpiX(1,tloop,pos(sloop)))
      A = stablepolyarea(...
        SET(no).ResolutionY*SET(no).EpiY(1:end-1,tloop,pos(sloop)),...
        SET(no).ResolutionX*SET(no).EpiX(1:end-1,tloop,pos(sloop)));
      EPVall(sloop,tloop)=A*(SET(no).SliceThickness+SET(no).SliceGap)/1000; %to cm3
    end;
  end;
end;

%Sum to get total volume
if length(pos)>1
  SET(no).LVV = sum(LVVall);
  SET(no).EPV = sum(EPVall);
else
  SET(no).LVV = LVVall;
  SET(no).EPV = EPVall;
end;
SET(no).LVM = SET(no).EPV-SET(no).LVV+SET(no).PV;
  
if nargout>0
  varargout = cell(1,1);
  varargout{1} = LVVall;
end;

if nargout>1
  varargout{2} = EPVall;
end;

if nargout>2
  varargout{3} = 1.05*(EPVall-LVVall); %LVM in g, NOTE not include Papillary volume (PV)
end;

%Update EDV,ESV
SET(no).EDV = SET(no).LVV(SET(no).EDT);
if ~isequal(SET(no).EDT,SET(no).EST)
  SET(no).ESV = SET(no).LVV(SET(no).EST);
else
  SET(no).ESV = NaN;
end

if (SET(no).ZSize>2)&&docomp
  LVVnocomp = SET(no).LVV;
  EPVnocomp = SET(no).EPV;

  for rloop=1:1 %Converges very quickly

    %Calculate compensation mechanism
    if (SET(no).LVV(SET(no).EDT)-SET(no).LVV(SET(no).EST)~=0)
      pp = SET(no).LVV;
      pp = pp-SET(no).LVV(SET(no).EST);
      pp = pp./(SET(no).LVV(SET(no).EDT)-SET(no).LVV(SET(no).EST));
      pp = 1-pp;
      pp = min(max(pp,0),1);
    else
      pp = ones(size(SET(no).LVV));
    end;
    err = zeros(1,20);

    %Compensate
    %pp = pp.*pp;

    %Check if autodetect
    if SET(no).AutoLongaxis
      for loop=1:20

        %Calculate slices
        slices = (loop-1); %Convert to mm
        slices = slices/(SET(no).SliceThickness+SET(no).SliceGap); %Convert to slices

        sloop=1;
        SET(no).LVV = LVVnocomp; %Restore
        SET(no).EPV = EPVnocomp;
        while (slices>0)&&(sloop<size(LVVall,1))
          %Remove whole slice or fraction of it.
          SET(no).LVV = SET(no).LVV-min(slices,1)*LVVall(sloop,:).*pp;
          SET(no).EPV = SET(no).EPV-min(slices,1)*EPVall(sloop,:).*pp;
          slices = slices-1;
          sloop = sloop+1;
        end;

        %Calculate error
        err(loop) = max(SET(no).EPV-SET(no).LVV+SET(no).PV)-min(SET(no).EPV-SET(no).LVV+SET(no).PV);
      end; %loop
      [~,inde] = min(err);
      SET(no).Longaxis = inde;
    end; %autodetect

    %Convert to slices => 1mm between, first is zero
    slices = (SET(no).Longaxis-1);
    if isempty(slices)
      slices = 0;
    end
    slices = slices/(SET(no).SliceThickness+SET(no).SliceGap);

    %Compensate the last time, now store it!
    sloop=1;
    SET(no).LVV = LVVnocomp; %Restore
    SET(no).EPV = EPVnocomp;

    %Check if need to fix with outline
    if ~isempty(SET(no).EndoX)
      zloop = find(not(isnan(SET(no).EndoX(1,SET(no).CurrentTimeFrame,:))));
    elseif ~isempty(SET(no).EpiX)
      zloop = find(not(isnan(SET(no).EpiX(1,SET(no).CurrentTimeFrame,:))));
    else
      zloop = [];
    end
    if isempty(zloop)
      slices = 0;
      zloop = SET(no).CurrentSlice;
    else
      zloop = zloop(1);
    end;

    %Make sure that it is updated
    %if docomp
    %  updatemodeldisplay;
    %end;
    
    while (slices>0)&&(sloop<=size(LVVall,1))
      %Remove whole slice or fraction of it.
      SET(no).LVV = SET(no).LVV-min(slices,1)*LVVall(sloop,:).*pp;
      SET(no).EPV = SET(no).EPV-min(slices,1)*EPVall(sloop,:).*pp;
     
      [xofs,yofs] = calcoffset(zloop,'montage');

      %Need to fix with outline?
      if 0%get(DATA.Handles.volumeoutlinecheckbox,'value')
        for tloop=1:SET(no).TSize
          if ~isempty(SET(no).EndoX)
            SET(no).EndoXView((1+(zloop-1)*(DATA.NumPoints+1)):(zloop*(DATA.NumPoints+1)-1),tloop) = ...
              (SET(no).EndoX(:,tloop,zloop)-SET(no).CenterX)*(1-pp(tloop)*min(1,slices))+SET(no).CenterX+xofs;
            SET(no).EndoYView((1+(zloop-1)*(DATA.NumPoints+1)):(zloop*(DATA.NumPoints+1)-1),tloop) = ...
              (SET(no).EndoY(:,tloop,zloop)-SET(no).CenterY)*(1-pp(tloop)*min(1,slices))+SET(no).CenterY+yofs;
          end
          if ~isempty(SET(no).EpiX)
            SET(no).EpiXView((1+(zloop-1)*(DATA.NumPoints+1)):(zloop*(DATA.NumPoints+1)-1),tloop) = ...
              (SET(no).EpiX(:,tloop,zloop)-SET(no).CenterX)*(1-pp(tloop)*min(1,slices))+SET(no).CenterX+xofs;
            SET(no).EpiYView((1+(zloop-1)*(DATA.NumPoints+1)):(zloop*(DATA.NumPoints+1)-1),tloop) = ...
              (SET(no).EpiY(:,tloop,zloop)-SET(no).CenterY)*(1-pp(tloop)*min(1,slices))+SET(no).CenterY+yofs;
          end
        end;
      else
        for tloop=1:SET(no).TSize
          if ~isempty(SET(no).EndoX)
            SET(no).EndoXView((1+(zloop-1)*(DATA.NumPoints+1)):(zloop*(DATA.NumPoints+1)-1),tloop) = ...
              SET(no).EndoX(:,tloop,zloop)+xofs;
            SET(no).EndoYView((1+(zloop-1)*(DATA.NumPoints+1)):(zloop*(DATA.NumPoints+1)-1),tloop) = ...
              SET(no).EndoY(:,tloop,zloop)+yofs;
          end
          if ~isempty(SET(no).EpiX)
            SET(no).EpiXView((1+(zloop-1)*(DATA.NumPoints+1)):(zloop*(DATA.NumPoints+1)-1),tloop) = ...
              SET(no).EpiX(:,tloop,zloop)+xofs;
            SET(no).EpiYView((1+(zloop-1)*(DATA.NumPoints+1)):(zloop*(DATA.NumPoints+1)-1),tloop) = ...
              SET(no).EpiY(:,tloop,zloop)+yofs;
          end
        end;
      end;
      zloop = zloop+1;
      slices = slices-1;
      sloop = sloop+1;
    end;
  end; %rloop
  
  if DATA.Silent
    return;
  end;
  
  for loop=1:length(DATA.ViewPanels)
    if isequal(no,DATA.ViewPanels(loop))
      if 0%isequal(get(DATA.Handles.volumeoutlinecheckbox,'value'),1) &&...
          ismember(DATA.ViewPanelsType{DATA.CurrentPanel},{'montage','montagerow','montagefit','sax3'})
        set(DATA.Handles.endocontour(loop),'LineStyle',':');
        set(DATA.Handles.epicontour(loop),'LineStyle',':');
      else
        set(DATA.Handles.endocontour(loop),'LineStyle','-');
        set(DATA.Handles.epicontour(loop),'LineStyle','-');
      end;
    end;
  end;

end; %If docompensation


%-----------------------------
function calclvvolumepolar(no)
%-----------------------------
%Calculate volume of lv when rotated image stacks.
%no is the imsage stack. The LV volume calculation includes
%longaxis compensation taken from the setting in the GUI and
%stored in the SET structure.

global SET

ind = (findfunctions('findslicewithendo',no))|(findfunctions('findslicewithepi',no));

if ~any(ind)
  SET(no).LVV = nan(1,SET(no).TSize);
  SET(no).EPV = SET(no).LVV;
  SET(no).LVM = SET(no).LVV;
  SET(no).PV = zeros(size(SET(no).LVV));
  SET(no).PFR = 0;
  SET(no).PER = 0;
  SET(no).PFRT = 1;
  SET(no).PERT = 1;
  SET(no).ESV = 0;
  SET(no).EDV = 0;
  SET(no).EF = 0;
  SET(no).SV = 0;
  return;
end;

%Find rotation axis, rotation around x-axis
my = SET(no).RotationCenter;

%'Angle' increment
alphapart = 1/(2*SET(no).ZSize);

%--- Calc endo volume
if ~isempty(SET(no).EndoX) && ~all(isnan(SET(no).EndoX(:)))
  %Calc dx/ds
  temp = double(SET(no).ResolutionX)*double(cat(1,SET(no).EndoX(:,:,ind),SET(no).EndoX(1,:,ind)))/10; %cm
  if ndims(temp)==2
    dxds = conv2(temp,[1;-1],'same');
    dxds = dxds(1:(end-1),:,:); %Remove "outside" data
  else
    dxds = econv3(temp,[-1;1]);
    dxds = dxds(2:end,:,:); %Different convolves have different "outside"
  end;
  y = double(SET(no).ResolutionY)*double((SET(no).EndoY(:,:,ind)-my))/10; %cm
  vol = pi*alphapart*y.^2.*sign(y).*dxds;
  vol = nansum(nansum(vol,1),3);
  SET(no).LVV = vol; %Store
else
  SET(no).LVV = zeros(1,SET(no).TSize);
end;

%--- Calc epi volume
if ~isempty(SET(no).EpiX) && ~all(isnan(SET(no).EpiX(:)))
  temp = double(SET(no).ResolutionX)*double(cat(1,SET(no).EpiX(:,:,ind),SET(no).EpiX(1,:,ind)))/10; %cm
  if ndims(temp)==2
    dxds = conv2(temp,[1;-1],'same');
    dxds = dxds(1:(end-1),:,:); %Remove "outside" data
  else
    dxds = econv3(temp,[-1;1]);
    dxds = dxds(2:end,:,:); %Different convolves have different "outside"
  end;
  y = double(SET(no).ResolutionY)*double(SET(no).EpiY(:,:,ind)-my)/10; %cm
  vol = pi*alphapart*y.^2.*sign(y).*dxds;
  vol = nansum(nansum(vol,1),3);
  SET(no).EPV = vol; %Store
else
  SET(no).EPV = zeros(1,SET(no).TSize);  
end;

%------------------------------------
function [varargout] = calcrvvolume(no)
%------------------------------------
%Calculate RV volume. no is the image stack.
%The RV volume calculation does not involve any longaxis 
%compensation.

global SET

needtodo = false;
if ~isempty(SET(no).RVEndoX)
  needtodo = true;
end;
if ~isempty(SET(no).RVEpiX)
  needtodo = true;
end;
if SET(no).Rotated
  calcrvvolumepolar(no);
  return;
end;
  
if ~needtodo
  SET(no).RVV = zeros(1,SET(no).TSize);
  SET(no).RVEPV = zeros(1,SET(no).TSize);
  SET(no).RVEDV = 0;
  SET(no).RVESV = 0;  
  SET(no).RVSV = 0;
  SET(no).RVEF = 0;
  SET(no).RVM = 0;
  varargout = cell(1,nargout);
  return;
end;

ind = (findfunctions('findslicewithrvendo',no))|(findfunctions('findslicewithrvepi',no)); %Accepts empty endo and epi if exist

%Find what slices to do
if SET(no).ZSize>1
  pos = find(ind);
else
  pos = 1;
end;
RVVall = zeros(length(pos),SET(no).TSize);
EPVall = zeros(length(pos),SET(no).TSize);

%Loop over all segmented slices
SET(no).RVV = zeros(1,SET(no).TSize);
if ~isempty(SET(no).RVEndoX)
  for sloop=1:length(pos)
    for tloop=1:SET(no).TSize
      if ~isnan(SET(no).RVEndoX(1,tloop,pos(sloop)))
        A = stablepolyarea(...
          SET(no).ResolutionY*SET(no).RVEndoY(1:end-1,tloop,pos(sloop)),...
          SET(no).ResolutionX*SET(no).RVEndoX(1:end-1,tloop,pos(sloop)));
        RVVall(sloop,tloop)=A*(SET(no).SliceThickness+SET(no).SliceGap)/1000;
      end;
    end;
  end;
end;
SET(no).RVV=sum(RVVall,1); %Buggfix for RV calculation. Addded ,1 EH:
SET(no).RVEDV = SET(no).RVV(SET(no).EDT);
SET(no).RVESV = SET(no).RVV(SET(no).EST);

SET(no).RVEPV = zeros(1,SET(no).TSize);
if ~isempty(SET(no).RVEpiX)
   for sloop=1:length(pos)
    for tloop=1:SET(no).TSize
      if ~isnan(SET(no).RVEpiX(1,tloop,pos(sloop)))
        A = stablepolyarea(...
          SET(no).ResolutionY*SET(no).RVEpiY(1:end-1,tloop,pos(sloop)),...
          SET(no).ResolutionX*SET(no).RVEpiX(1:end-1,tloop,pos(sloop)));
        EPVall(sloop,tloop)=A*(SET(no).SliceThickness+SET(no).SliceGap)/1000;
      end;
    end;
  end;
end;
SET(no).RVEPV=sum(EPVall,1); %Buggfix for RV calculation. Addded ,1 EH:

if nargout>0
  varargout = cell(1,1);
  varargout{1} = RVVall;
end;

if nargout>1
  varargout{2} = EPVall;
end;

%-----------------------------
function calcrvvolumepolar(no)
%-----------------------------
%Calculate RV volume when rotated image stacks.

global SET

ind = (findfunctions('findslicewithrvendo',no))|(findfunctions('findslicewithrvepi',no));

if ~any(ind)
  SET(no).RVV = nan(1,SET(no).TSize);
  SET(no).RVEPV = SET(no).RVV;
  SET(no).RVM = SET(no).RVV;
  SET(no).RVESV = 0;
  SET(no).RVEDV = 0;
  SET(no).RVEF = 0;
  SET(no).RVSV = 0;
  return;
end;

%Find rotation axis, rotation around x-axis
my = SET(no).RotationCenter;

%'Angle' increment
alphapart = 1/(2*SET(no).ZSize);

%--- Calc endo volume
if ~isempty(SET(no).RVEndoX)
  %Calc dx/ds
  temp = SET(no).ResolutionX*cat(1,SET(no).RVEndoX(:,:,ind),SET(no).RVEndoX(1,:,ind))/10; %cm
  if ndims(temp)==2
    dxds = conv2(temp,[1;-1],'same');
    dxds = dxds(1:(end-1),:,:); %Remove "outside" data
  else
    dxds = econv3(temp,[-1;1]);
    dxds = dxds(2:end,:,:); %Different convolves have different "outside"
  end;
  y = SET(no).ResolutionY*(SET(no).RVEndoY(:,:,ind)-my)/10; %cm
  vol = pi*alphapart*y.^2.*sign(y).*dxds;
  vol = nansum(nansum(vol,1),3);
  SET(no).RVV = vol; %Store
else
  SET(no).RVV = zeros(1,SET(no).TSize);
end;

%--- Calc epi volume
if ~isempty(SET(no).RVEpiX)
  temp = SET(no).ResolutionX*cat(1,SET(no).RVEpiX(:,:,ind),SET(no).RVEpiX(1,:,ind))/10; %cm
  if ndims(temp)==2
    dxds = conv2(temp,[1;-1],'same');
    dxds = dxds(1:(end-1),:,:); %Remove "outside" data
  else
    dxds = econv3(temp,[-1;1]);
    dxds = dxds(2:end,:,:); %Different convolves have different "outside"
  end;
  y = SET(no).ResolutionY*(SET(no).RVEpiY(:,:,ind)-my)/10; %cm
  vol = pi*alphapart*y.^2.*sign(y).*dxds;
  vol = nansum(nansum(vol,1),3);
  SET(no).RVEPV = vol; %Store
else
  SET(no).RVEPV = zeros(1,SET(no).TSize);  
end;

%-------------------------
function volume_helper(no)
%-------------------------
%Helper to find peak ejection, and empty volumes etc.
%Lots of checks to prevent NaN or Inf to be presented
%when some of the data are missing.

global SET 

%Find zeros in volume
SET(no).LVV(SET(no).LVV==0) = NaN;
SET(no).EPV(SET(no).EPV==0) = NaN;
SET(no).RVV(SET(no).RVV==0) = NaN;
SET(no).RVEPV(SET(no).RVEPV==0) = NaN;
SET(no).LVM = SET(no).EPV-SET(no).LVV+SET(no).PV;
SET(no).RVM = SET(no).RVEPV-SET(no).RVV;

%Prevent that when endo is empty that epi may be viewed
if not(isnan(SET(no).EPV(1)))&&(isnan(SET(no).LVV(1)))
  SET(no).LVV = zeros(size(SET(no).LVV));
end;

%Find peak ejection rate, peak filling rate
if SET(no).TSize>2
  %dLVV = conv2(SET(no).LVV,[1 -1]/SET(no).TIncr,'valid');
  dLVV = diff(SET(no).LVV)./diff(SET(no).TimeVector);
  [SET(no).PFR,SET(no).PFRT] = max(dLVV); %Peak filling rate
  [SET(no).PER,SET(no).PERT] = min(dLVV); %Peak ejection rate
  SET(no).PER = -SET(no).PER;
  dRVV = diff(SET(no).RVV)./diff(SET(no).TimeVector);
  [SET(no).RVPFR,SET(no).RVPFRT] = max(dRVV); %Peak filling rate
  [SET(no).RVPER,SET(no).RVPERT] = min(dRVV); %Peak ejection rate
  SET(no).RVPER = -SET(no).RVPER;
else
  SET(no).PFR = 0;
  SET(no).PER = 0;
  SET(no).PFRT = 1;
  SET(no).PERT = 1;
  SET(no).RVPFR = 0;
  SET(no).RVPER = 0;
  SET(no).RVPFRT = 1;
  SET(no).RVPERT = 1;
end;

%Store LV-EDV,ESV,SV
SET(no).EDV = SET(no).LVV(SET(no).EDT)-SET(no).PV(SET(no).EDT);
if ~isequal(SET(no).EDT,SET(no).EST)
  SET(no).ESV = SET(no).LVV(SET(no).EST)-SET(no).PV(SET(no).EST);
  SET(no).SV = SET(no).EDV-SET(no).ESV; %Stroke volume
else
  SET(no).ESV=NaN;
  SET(no).SV=NaN;
end

%LV-EF
if SET(no).EDV>0
  SET(no).EF = SET(no).SV/SET(no).EDV; %Ejection fraction
else
  SET(no).EF = NaN;
end;

if ~isempty(SET(no).RVEndoX)
  %Store RV-EDV,ESV,SV
  SET(no).RVEDV = SET(no).RVV(SET(no).EDT);
  SET(no).RVESV = SET(no).RVV(SET(no).EST);
  SET(no).RVSV = SET(no).RVEDV-SET(no).RVESV; %Stroke volume
  
  %RV-EF
  if SET(no).RVEDV==0 || isnan(SET(no).RVEDV) || isempty(SET(no).RVEDV)
    SET(no).RVEF = 0;
  else
    SET(no).RVEF = SET(no).RVSV/SET(no).RVEDV; %Ejection fraction
  end;
end;

if isnan(SET(no).RVEDV)
  SET(no).RVEDV = 0;
end;
if isnan(SET(no).RVESV)
  SET(no).RVESV = 0;
end;
if isnan(SET(no).RVSV)
  SET(no).RVSV = 0;
end;
if isnan(SET(no).RVEF)
  SET(no).RVEF = 0;
end;
if isequal(SET(no).EDT,SET(no).EST)
  SET(no).RVESV=NaN;
  SET(no).RVSV=NaN;
  SET(no).RVEF=NaN;
end

%-------------------------------------
function bsa = calcbsa(weight,height) %#ok<DEFNU>
%-------------------------------------
%Calculates BSA. Formula based on Mosteller
%weight in kilo and height in cm.
bsa = sqrt(weight*height/3600);

%-----------------------------------
function ticks = calcticks(num,res) %#ok<DEFNU>
%-----------------------------------
%Helper fcn to calculate length of ticks in volumegraph

totlength = floor(num*res/10-1e-4); %length in cm
ticks = 10*(1:totlength); %ticks i mm
ticks = ticks/res; %ticks i pixels

%--------------------------
function im = truedata2im(z,no) %#ok<DEFNU>
%--------------------------
%Inverse of calctruedata

global SET NO

if nargin<2
  no = NO;
end;
if not(isa(z,'int16')) && ~isempty(SET(no).IntensityOffset)
  %z = im*SET(no).IntensityScaling+SET(no).IntensityOffset;
  im = (z-SET(no).IntensityOffset)/SET(no).IntensityScaling;
else
  im = z;
end

%-------------------------------
function z = calctruedata(im,no) 
%-------------------------------
%Calculate true image intensities (as before Segment internal
%normalization). Uses IntensityScaling and IntensityOffset stored
%in SET structure. im is input image, and no is image stack, 
%where to take the scaling from.

global SET NO

if nargin<2
  no = NO;
end;

if ~isempty(SET(no).IntensityScaling)
  if not(isa(im,'int16'))
    z = im*SET(no).IntensityScaling+SET(no).IntensityOffset;
  else
    z = single(im);
  end
elseif ~isempty(SET(no).VENC) && SET(no).VENC ~= 0
  z = (im-0.5)*2*SET(no).VENC;
else
  z = im;
end;

%----------------------------------------
function radvel = calcradialvelocity(no) %#ok<DEFNU>
%----------------------------------------
%Calculate radial velocity of the endocardium. 
%Forward difference is used.

global SET NO

if nargin<1
  no = NO;
end;

%Calc radius
rad = calcendoradius(no);

%Derivate
if SET(no).Cyclic
  temp = econvn(cat(2,rad,rad(:,1,:)),[1 -1]/SET(no).TIncr);
else
  if SET(no).TIncr~=0
    temp = econvn(rad,[1 -1]/SET(no).TIncr);
    temp = cat(2,temp,temp(:,end,:));
  else
    temp = zeros(size(rad));
  end;
end;
temp = temp(:,2:end,:); %Remove first timeframe, wrong derivate

radvel = temp/10; %=> cm/s

%--------------------------
function [meanarea,area]=calcroiarea(no,roino,thisframeonly) %#ok<DEFNU>
%--------------------------
%Calculates roi area (helper fcn)

global SET

if nargin < 3
  thisframeonly = false;
end
if thisframeonly
  tvec = SET(no).CurrentTimeFrame;
else
  tvec = SET(no).Roi(roino).T;
end

area=nan(1,SET(no).TSize);
for tloop=tvec
    if not(isnan(SET(no).Roi(roino).Y(1,tloop)))
      area(tloop)= (1/100)*stablepolyarea(...
        SET(no).ResolutionY*SET(no).Roi(roino).Y(:,tloop),...
        SET(no).ResolutionX*SET(no).Roi(roino).X(:,tloop));
    end
end
if thisframeonly
  area = area(tvec);
end
meanarea=mynanmean(area);

%-----------------------------------------------------------------------------
function [m,sd,rmin,rmax] = calcroiintensity(no,roino,normalize,thisframeonly) %#ok<DEFNU>
%-----------------------------------------------------------------------------
%Calculates intensity within a ROI (helper fcn)
global SET
if nargin < 3
  normalize = false;
end
if nargin < 4
  thisframeonly = false;
end

if thisframeonly
  tvec = SET(no).CurrentTimeFrame;
else
  tvec = SET(no).Roi(roino).T;
end

m = nan(1,SET(no).TSize);
sd = m;
rmin = m;
rmax = m;
z = SET(no).Roi(roino).Z;
for tloop=tvec
  if not(isnan(SET(no).Roi(roino).Y(1,tloop)))
    if normalize
      temp = SET(no).IM(:,:,tloop,z);
    else
      temp = calctruedata(SET(no).IM(:,:,tloop,z),no);
    end;
    roimask = segment('createmask',...
      [SET(no).XSize SET(no).YSize],...
      SET(no).Roi(roino).Y(:,tloop),...
      SET(no).Roi(roino).X(:,tloop));
    ind = find(roimask);
    if ~isempty(ind)
      m(tloop) = mean(temp(ind));
      sd(tloop) = std(temp(ind));
      rmin(tloop) = min(temp(ind));
      rmax(tloop) = max(temp(ind));
    end;
  end
end

if thisframeonly
  m = m(tvec);
  sd = sd(tvec);
  rmin = rmin(tvec);
  rmax = rmax(tvec);
end


%-----------------------------------------------
function res = calcmyocardvolume(numsectors,no,tf) %#ok<DEFNU>
%-----------------------------------------------
%Calculate myocardvolume. numsectors is the number
%of sectors to calculate in, and no is the image stack.
%uses findinmeansectorslice to do some of the calculation.

global DATA SET NO

if nargin<2
  no = NO;
end;
if nargin < 3
  tf = SET(no).CurrentTimeFrame;
end
  
res = nan([numsectors SET(no).ZSize]);

if isempty(SET(no).EndoX)
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end;

if isempty(SET(no).EpiX)
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end;

numpoints = DATA.Pref.RadialProfiles;

%Upsample model
if not(numpoints==DATA.NumPoints)
  [endox,endoy] = resamplemodel(SET(no).EndoX(:,tf,:),SET(no).EndoY(:,tf,:),numpoints);
  [epix,epiy] = resamplemodel(SET(no).EpiX(:,tf,:),SET(no).EpiY(:,tf,:),numpoints);
else
  endox = SET(no).EndoX(:,tf,:);
  endoy = SET(no).EndoY(:,tf,:);
  epix = SET(no).EpiX(:,tf,:);
  epiy = SET(no).EpiY(:,tf,:);
end;

for slice=1:SET(no).ZSize  
  if not(isnan(epix(1,slice)))
    if numsectors == 1
      res(1,slice) = 1e-3*(SET(no).SliceThickness+SET(no).SliceGap)*...
        (stablepolyarea(...
        SET(no).ResolutionY*epiy(:,:,slice),...
        SET(no).ResolutionX*epix(:,:,slice)) - ...
        stablepolyarea(...
        SET(no).ResolutionY*endoy(:,:,slice),...
        SET(no).ResolutionX*endox(:,:,slice)));
    else
      %Find sectors
      [meanxepi,meanyepi,episectors] = findmeaninsectorslice('epi',numpoints,tf,slice,numsectors,no);
      
      if isnan(endox(1,slice))
        meanxendo = meanxepi;
        meanyendo = meanyepi;
      else
        [meanxendo,meanyendo,endosectors] = findmeaninsectorslice('endo',numpoints,tf,slice,numsectors,no);
      end;
      
      for loop=1:numsectors
        if episectors(loop)<=episectors(loop+1)
          tempepix = epix(episectors(loop):episectors(loop+1),slice);
          tempepiy = epiy(episectors(loop):episectors(loop+1),slice);
        else
          tempepix = [...
            epix(episectors(loop):numpoints,slice) ; ...
            epix(1:episectors(loop+1),slice)];
          tempepiy = [...
            epiy(episectors(loop):numpoints,slice) ; ...
            epiy(1:episectors(loop+1),slice)];
        end;
        
        if isnan(endox(1,slice))
          tempendox = repmat(meanxendo,size(tempepix));
          tempendoy = repmat(meanyendo,size(tempepiy));
        else
          if endosectors(loop)<=endosectors(loop+1)
            tempendox = endox(endosectors(loop):endosectors(loop+1),slice);
            tempendoy = endoy(endosectors(loop):endosectors(loop+1),slice);
          else
            tempendox = [...
              endox(endosectors(loop):numpoints,slice) ; ...
              endox(1:endosectors(loop+1),slice)];
            tempendoy = [...
              endoy(endosectors(loop):numpoints,slice) ; ...
              endoy(1:endosectors(loop+1),slice)];
          end;
        end;
        
        %Calc volume in [ml]
        res(loop,slice) = 1e-3*(SET(no).SliceThickness+SET(no).SliceGap)*...
          stablepolyarea(...
          SET(no).ResolutionY*[tempendoy;flipud(tempepiy);tempendoy(1)],...
          SET(no).ResolutionX*[tempendox;flipud(tempepix);tempendox(1)]);
        
      end;
      %correct the volume so the total volume in the current slice is correct
      if ~isnan(endox(1,slice))
        slicevolume = 1e-3*(SET(no).SliceThickness+SET(no).SliceGap)*...
          stablepolyarea(...
          SET(no).ResolutionY*[endoy(:,slice);flipud(epiy(:,slice));endoy(1,slice)],...
          SET(no).ResolutionX*[endox(:,slice);flipud(epix(:,slice));endox(1,slice)]);
      else
        slicevolume = 1e-3*(SET(no).SliceThickness+SET(no).SliceGap)*...
          stablepolyarea(...
          SET(no).ResolutionY*[tempendoy;flipud(epiy(:,slice));tempendoy],...
          SET(no).ResolutionX*[tempendox;flipud(epix(:,slice));tempendox]);
      end
      correction = slicevolume/sum(res(:,slice));
      res(:,slice) = res(:,slice)*correction;
    end
  end;
end;

%---------------------------------
function rad = calcendoradius(no)
%---------------------------------
%Calculates endocardial radius. Respects setting in preferences
%if to use endocardial or epicardial center as reference. Loops
%over timeframes and slices.

global DATA SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).EndoX)
  rad = zeros(DATA.Pref.RadialProfiles,SET(no).TSize,SET(no).ZSize);
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end;

%Upsample model
if not(DATA.Pref.RadialProfiles==DATA.NumPoints)
  [endox,endoy] = resamplemodel(SET(no).EndoX,SET(no).EndoY,DATA.Pref.RadialProfiles);
else
  endox = SET(no).EndoX;
  endoy = SET(no).EndoY;
end;

rad = zeros(DATA.Pref.RadialProfiles,SET(no).TSize,SET(no).ZSize);
for zloop=1:SET(no).ZSize
  for tloop=1:SET(no).TSize    
    if SET(no).EndoCenter
      %Calc mean
      meanx = mean(SET(no).EndoX(:,tloop,zloop));
      meany = mean(SET(no).EndoY(:,tloop,zloop));
    else
      meanx = mean(SET(no).EpiX(:,tloop,zloop));
      meany = mean(SET(no).EpiY(:,tloop,zloop));
    end;
    
    %Calc radius
    rad(:,tloop,zloop) = sqrt(...
      ((endox(:,tloop,zloop)-meanx)*SET(no).ResolutionX).^2+...
      ((endoy(:,tloop,zloop)-meany)*SET(no).ResolutionY).^2);
  end;
end;

%--------------------------------
function rad = calcepiradius(no)
%--------------------------------
%Calculates epicardial radius. Respects setting if to use
%endocardial or epicardial center in calculation. Loops over
%timeframes and slices.

global DATA SET NO

if nargin<1
  no = NO;
end;

if isempty(SET(no).EpiX)
  rad = zeros(DATA.Pref.RadialProfiles,SET(no).TSize,SET(no).ZSize);
  myfailed('No LV epicardium available.',DATA.GUI.Segment);
  return;
end;

%Upsample model
if not(DATA.Pref.RadialProfiles == DATA.NumPoints)
  [epix,epiy] = resamplemodel(SET(no).EpiX,SET(no).EpiY,DATA.Pref.RadialProfiles);
else
  epix = SET(no).EpiX;
  epiy = SET(no).EpiY;
end;

rad = zeros(DATA.Pref.RadialProfiles,SET(no).TSize,SET(no).ZSize);
for zloop=1:SET(no).ZSize
  for tloop=1:SET(no).TSize
    %Calc mean
    if SET(no).EndoCenter && ~isnan(SET(no).EndoX(1,tloop,zloop))
      %Calc mean
      meanx = mean(SET(no).EndoX(:,tloop,zloop));
      meany = mean(SET(no).EndoY(:,tloop,zloop));    
    else
      meanx = mean(SET(no).EpiX(:,tloop,zloop));
      meany = mean(SET(no).EpiY(:,tloop,zloop));
    end;
  
    %Calc radius
    rad(:,tloop,zloop) = sqrt(...
      (SET(no).ResolutionX*(epix(:,tloop,zloop)-meanx)).^2+...
      (SET(no).ResolutionY*(epiy(:,tloop,zloop)-meany)).^2);
  end;
end;

%------------------------------------------------------
function [wallthickness,endox,endoy,epix,epiy] = calcwallthickness(sectors,no) %#ok<DEFNU>
%------------------------------------------------------
%Calculate wallthickness. Uses calcendoradius and calcepiradius
%to do the work.

global DATA SET NO

if nargin<2
  no = NO;
end

if isempty(SET(no).EndoX)
  myfailed('No LV endocardium available.',DATA.GUI.Segment);
  return;
end;

if isempty(SET(no).EpiX)
  myfailed('No LV epicardium available.',DATA.GUI.Segment);
  return;
end;

pos = find(findfunctions('findslicewithendo',no)&findfunctions('findslicewithepi',no));
nslices = length(pos);

if length(nslices)<1
  myfailed('Need more than one segmented slice.',DATA.GUI.Segment);
  return;
end;

%Calc radius
endorad = calcendoradius(no);
epirad = calcepiradius(no);
[endorad,~,~,sector_inds_endo,~] = findmeaninsector('endo',endorad,1:SET(no).ZSize,sectors,no);
[epirad,~,~,sector_inds_epi,~] = findmeaninsector('epi',epirad,1:SET(no).ZSize,sectors,no);

%If endo is missing in the most apical slices, fill in endo in the most apical slices
for t = 1:size(endorad,3)
  lastendo = find(~isnan(endorad(1,:,t)),1,'last');
  lastepi = find(~isnan(epirad(1,:,t)),1,'last');
  endorad(:,lastendo+1:lastepi,t) = 0;
end

%Wallthickness determed as mean epiradius minus mean endoradius
%within the sector.
wallthickness = epirad-endorad;


if nargout==5
  sector_inds_endo = sector_inds_endo(1:end-1,:,:);
  sector_inds_epi = sector_inds_epi(1:end-1,:,:);
  sector_inds_endo = sector_inds_endo.*(~isnan(endorad));
  sector_inds_epi = sector_inds_epi.*(~isnan(epirad));
  
  %sector_inds_endo(sector_inds_endo==0)=nan;
  %sector_inds_epi(sector_inds_epi==0)=nan;
  endox= nan(size(sector_inds_endo,1),size(sector_inds_endo,3),size(sector_inds_endo,2));
  endoy=endox;
  epix=endox;
  epiy=endox;
  
  for t = 1:SET(no).TSize
    for z = 1:SET(no).ZSize
      tmp = find(sector_inds_endo(:,z,t));
      if ~isempty(tmp)
        for p =sector_inds_endo(:,z,t)
          endox(:,t,z) = SET(no).EndoX(p,t,z);
          endoy(:,t,z) = SET(no).EndoY(p,t,z);
        end
        
        for p = sector_inds_epi(:,z,t)
          epix(:,t,z) = SET(no).EpiX(p,t,z);
          epiy(:,t,z) = SET(no).EpiY(p,t,z);
        end
      end
    end
  end
end

%---------------------------------------------------------------------
function [xout,yout] = calcsegmentationintersections(no,type,viewtype) %#ok<DEFNU>
%---------------------------------------------------------------------
%Calculates intersections between segmentation in image stacks.
%no is image stack and type is endo or epi, can also be rvendo etc.
%the type is a string that dynamically calls a field of the SET struct.

global DATA SET

if nargin < 3
  viewtype = 'one';
end

%Always return something
xout = nan;
yout = nan;
espot = regexp(type,'[e]');
xfield = [upper(type(1:espot)) type(espot+1:end) 'X'];
yfield = [upper(type(1:espot)) type(espot+1:end) 'Y'];

if DATA.Silent
  return;
end;

%Find existing segmentation
panelix = find(DATA.ViewPanels>0);
noseg = zeros(size(panelix));
typeseg = cell(size(panelix));

for loop=panelix
  noloop = DATA.ViewPanels(loop);
  if ~isempty(SET(noloop).(xfield))
    noseg(loop) = noloop;
    typeseg{loop} = DATA.ViewPanelsType{loop};
    %break %This causes segmentation from only one stack to be displayed
  end;
end;
typeseg = typeseg(noseg>0);
noseg = noseg(noseg>0);

if length(noseg)<1
  return;
end;
  
%Find equation for no plane
if ~ismember(viewtype,{'hla','vla','gla'})
  [pos,~,~,zdir] = getimageposandorientation(no,SET(no).CurrentSlice,viewtype);
else
  [pos,~,~,zdir] = getimageposandorientation(no,SET(no).(upper(viewtype)).slice,viewtype);
  if strcmp(viewtype,'gla')
    ang = SET(no).GLA.angle;
    res = SET(no).ResolutionY*cos(ang)+SET(no).ResolutionX*abs(sin(ang));
  end
end
d = sum(pos.*zdir);

xbuild = [];
ybuild = [];
for loop=1:length(noseg)
  if ~isequal(noseg(loop),no) || ~isequal(typeseg{loop},viewtype)
    x = SET(noseg(loop)).(xfield);
    y = SET(noseg(loop)).(yfield);

    if ~isempty(x)% || any(findfunctions(['findslicewith',type,'all'],noseg(loop)))
      %calculate corresponding timeframe
      if SET(no).TSize == 2
        if SET(no).CurrentTimeFrame == SET(noseg(loop)).CurrentTimeFrame
          tfs = SET(no).CurrentTimeFrame;
        else
          tfs = [];
          tfdiff = [];
        end
       elseif any(findfunctions(['findslicewith',type,'all'],noseg(loop)))
        alltf = (1+((1:40)-1)/(SET(noseg(loop)).TSize-1)*(SET(no).TSize-1));
        [~,closestsegind] = min(abs(SET(no).CurrentTimeFrame-round(alltf)));
        tfs = max(1,min(SET(noseg(loop)).TSize,closestsegind));
        tfdiff = abs(alltf(tfs)-SET(no).CurrentTimeFrame);  
      elseif SET(no).TSize>2
        %find a corresponding (closest) tf in no, for each tf in noseg
        alltf = (1+((1:40)-1)/(SET(noseg(loop)).TSize-1)*(SET(no).TSize-1));
        tfs = max(1,min(SET(noseg(loop)).TSize,find(SET(no).CurrentTimeFrame==round(alltf))));
        tfdiff = abs(alltf(tfs)-SET(no).CurrentTimeFrame);
%         tf = round(1+(SET(no).CurrentTimeFrame-1)/(SET(no).TSize-1)*(SET(noseg(loop)).TSize-1));
        %[~,tf] = min(abs(SET(noseg(loop)).TimeVector-SET(no).TimeVector(SET(no).CurrentTimeFrame)));
     
       % 
      else
        tfs = 1;
      end
      %calculate intersection points for all slices
      for zloop=1:SET(noseg(loop)).ZSize
        if ~isnan(x(1,SET(noseg(loop)).CurrentTimeFrame,zloop))
          if length(tfs) == 1
            pos = calcfunctions('xyz2rlapfh',...
              noseg(loop),...
              x(:,tfs,zloop),...
              y(:,tfs,zloop),...
              repmat(zloop,length(x(:,tfs,zloop)),1));%DATA.NumPoints,1));
          else
            for posloop = 1:length(tfs)
              temppos{posloop} = calcfunctions('xyz2rlapfh',...
                noseg(loop),...
                x(:,tfs(posloop),zloop),...
                y(:,tfs(posloop),zloop),...
                repmat(zloop,length(x(:,tfs(posloop),zloop)),1));%DATA.NumPoints,1));
            end
            [~,tfsorder] = sort(tfdiff);
            for posloop = tfsorder
              thispos = temppos{posloop};
              if ~isempty(thispos) && ~isnan(thispos(1,1))
                pos = temppos{posloop};
                break;
              end
            end
          end

          if ~isempty(pos) %pos may be empty from rotated stacks.
            %Find crossings
            disttoplane = pos*(zdir')-d;
            ind = find((disttoplane(1:(end-1)).*disttoplane(2:end))<0);

            if ~isempty(ind)
              crosspos = 0.5*(pos(ind,:)+pos(ind+1,:));

              xyzno = rlapfh2xyz(no,crosspos(:,1),crosspos(:,2),crosspos(:,3));
              
              switch viewtype
                case 'hla'
                  xbuild = [xbuild xyzno(3,:)]; %#ok<AGROW>
                  ybuild = [ybuild xyzno(2,:)]; %#ok<AGROW>
                case 'vla'
                  xbuild = [xbuild xyzno(3,:)]; %#ok<AGROW>
                  ybuild = [ybuild xyzno(1,:)]; %#ok<AGROW>
                case 'gla'
                  xbuild = [xbuild xyzno(3,:)]; %#ok<AGROW>
                  [~,glay] = sax2gla(xyzno(1,:),xyzno(2,:),0,noseg(loop));
                  ybuild = [ybuild glay]; %#ok<AGROW>
                otherwise
                  xbuild = [xbuild xyzno(1,:)]; %#ok<AGROW>
                  ybuild = [ybuild xyzno(2,:)]; %#ok<AGROW>
              end
            end;
          end;
          
        end; %empty slice
      end; %zloop
      if ~isempty(xbuild)
        xout = xbuild;
        yout = ybuild;
      end
    end; %empty segmentation
  end; %equal segmation
end;

%--------------------------------------------------------------------
function [x,y] = calcplaneintersections(NO,no,TYPE,type,slice,hslice) %#ok<DEFNU>
%--------------------------------------------------------------------
%Calculate intersections between image planes.
%  NO is the viewed plane
%  no is the (potentially) intersecting plane
%  TYPE is the viewtype of the viewed plane
%  type is the viewtype of the intersecting plane
%  slice is the slice of the intersecting plane
%  hslice is the horizontal slice of the viewed plane (if derived longaxis)

global SET

x = [];
y = [];

if isequal(no,0)
  return;
end;
if isequal(NO,0)
  return;
end;

if nargin < 5
  slice = SET(no).CurrentSlice;
end
if nargin < 4
  TYPE = 'one';
  type = 'one';
end

xszNO = SET(NO).XSize;
xresNO = SET(NO).ResolutionX;
yszNO = SET(NO).YSize;
yresNO = SET(NO).ResolutionY;
%--- Find equation of planes no and NO
if ~ismember(type(1:3),{'hla','vla','gla'})
  [posno,~,~,zdirno] = getimageposandorientation(no,slice,type);
else
  [posno,~,~,zdirno] = getimageposandorientation(no,SET(no).(upper(type(1:3))).slice,type(1:3));
end
%Plane equation is ax+by+cz=d, a,b,c is given direcly by zdir (Kossan).
d = zdirno*(posno'); %d=ax+by+cz

%---- Find equation of plane NO
if ~ismember(TYPE(1:3),{'hla','vla','gla'})
  [posNO,xdirNO,ydirNO] = getimageposandorientation(NO,SET(NO).CurrentSlice,TYPE);
else
  if nargin < 6
    hslice = SET(NO).(upper(TYPE(1:3))).slice;
  end
  xszNO = SET(NO).ZSize;
  xresNO = SET(NO).SliceThickness + SET(NO).SliceGap;
  if strcmp(TYPE(1:3),'vla')
    yszNO = SET(NO).XSize;
    yresNO = SET(NO).ResolutionX;
  elseif strcmp(TYPE(1:3),'gla')
    glaangle = SET(NO).GLA.angle;
    res = SET(NO).ResolutionY*cos(glaangle)+SET(NO).ResolutionX*abs(sin(glaangle));
    yszNO = floor(SET(NO).YSize*abs(cos(glaangle))+SET(NO).XSize*abs(sin(glaangle)));
    yresNO = res;    
  end
%   if strcmp(type,'vla')
%     xszNO = SET(NO).YSize;
%     xresNO = SET(NO).ResolutionY;
%     yszNO = SET(NO).ZSize;
%     yresNO = SET(NO).SliceThickness + SET(NO).SliceGap;
%   end
  [posNO,xdirNO,ydirNO] = getimageposandorientation(NO,hslice,TYPE(1:3));
end


%lines are given orientation + zero pos, i.e k*dir+pos. k is number of mm
%along the line.

% zeropos
% +--- la ---+    => ImageOrientation(1:3)
% |          |    ||
% lb         ld   v ImageOrientation(4:6)
% |          |
% +--- lc ---+

ladir = ydirNO;
lapos = posNO; 
maxka = (yszNO-1)*yresNO; %in mm

lbdir = xdirNO;
lbpos = posNO; 
maxkb = (xszNO-1)*xresNO;

lcdir = ydirNO;
lcpos = posNO+xdirNO*(xszNO-1)*xresNO; 
maxkc = (yszNO-1)*yresNO;

lddir = xdirNO;
ldpos = posNO+ydirNO*(yszNO-1)*yresNO;
maxkd = (xszNO-1)*xresNO;

%The following equation will find the points on the line that
%fulfills the plane equation.
%
%k*dir*(zdirno')+pos*(zdirno') = d
%
%simplifying to and solving gives:
%k*q+p=0 => k*q=-p => k=-p/q.
%
%q = dir*(zdirno');
%p = pos*(zdirno')-d;

%Line a
qa = ladir*(zdirno');
pa = lapos*(zdirno')-d;
if abs(qa)>1e-3
  ka = -pa/qa;
else
  ka = NaN;
end;
if ka<0
  ka = NaN;
end;
if ka>maxka
  ka = NaN;
end;

%Line b
qb = lbdir*(zdirno');
pb = lbpos*(zdirno')-d;
if abs(qb)>1e-3
  kb = -pb/qb;
else
  kb = NaN;
end;
if kb<0
  kb = NaN;
end;
if kb>maxkb
  kb = NaN;
end;

%Line c
qc = lcdir*(zdirno');
pc = lcpos*(zdirno')-d;
if abs(qc)>1e-3
  kc = -pc/qc;
else
  kc = NaN;
end;
if kc<0
  kc = NaN;
end;
if kc>maxkc
  kc = NaN;
end;

%Line d
qd = lddir*(zdirno');
pd = ldpos*(zdirno')-d;
if abs(qd)>1e-3
  kd = -pd/qd;
else
  kd = NaN;
end;
if kd<0
  kd = NaN;
end;
if kd>maxkd
  kd = NaN;
end;

% zeropos
% +- => ka --+    => ImageOrientation(1:3)
% |          |    ||
% kb         kd   v ImageOrientation(4:6)
% |          |
% +- => kc --+

%Convert to top row is x, bottom row is y (segment coordinates)
temp = [...
  [0; ka/yresNO] ...            %ka
  [kb/xresNO;0] ....            %kb
  [xszNO-1;kc/yresNO] ... %kc
  [kd/xresNO;yszNO-1]];   %kd

ind = ~isnan(sum(temp,1));
temp = temp(:,ind);
x = temp(1,:)+1; %segment has 1 as origo
y = temp(2,:)+1;

%---------------------------------------
function [rows,cols] = calcrowscols(no,z)
%---------------------------------------
%Calculate a good montage setup, so maxiumum number of slices
%can be displayed on the screen.

global SET

if isequal(no,0);
  no = 1; 
end;

if nargin<2
  z = SET(no).ZSize;
end;

switch z
  case 1
    rows = 1;
    cols = 1;
  case 3
    rows = 1; cols = 3;
  case 6
    rows = 2; cols = 3;
  case {10,11,12}
    rows = 3; cols = 4;
  case {13,14,15}
    rows = 4; cols = 4;
  case {17,18,19,20}
    rows = 4; cols = 5;
  otherwise
    rows = round(sqrt(z));
    cols = ceil(z/rows);
end;

%-------------------------------------------------
function [xofs,yofs] = calcoffset(z,type,no,panel)
%-------------------------------------------------
%Calculate offset required to plot coordinates in viewing 
%mode specified by type.

global DATA SET NO

if nargin < 4
  panel = DATA.CurrentPanel;
end
if nargin < 3
  no = NO;
end
if nargin<2 || isempty(type)
  type = DATA.ViewPanelsType{panel};
end;  

if panel > numel(DATA.ViewPanelsMatrix) || ...
    isempty(DATA.ViewPanelsMatrix{panel})
  [rows,cols] = calcrowscols(no);
  DATA.ViewPanelsMatrix{panel}(1) = rows;
  DATA.ViewPanelsMatrix{panel}(2) = cols;
end;

switch type
  case {'montage','montagerow','montagefit'}
    c = 1+mod(z-1,DATA.ViewPanelsMatrix{panel}(2));
    r = ceil(z/DATA.ViewPanelsMatrix{panel}(2));
    yofs = (c-1)*SET(no).YSize;
    xofs = (r-1)*SET(no).XSize;
  case 'sax3'
    xofs = zeros(1,SET(no).TSize);
    yofs = zeros(1,SET(no).TSize);
    for t = 1:SET(no).TSize
      zind = find(SET(no).SAX3.slices(:,t) == z);
      if isempty(zind)
        xofs(t) = nan;
        yofs(t) = nan;
      else
        c = 1+mod(zind-1,DATA.ViewPanelsMatrix{panel}(2));
        r = ceil(zind/DATA.ViewPanelsMatrix{panel}(2));
        yofs(t) = (c-1)*SET(no).YSize;
        xofs(t) = (r-1)*SET(no).XSize;
      end
    end
  case 'montagesegmented'
    slicestoinclude = segment_main('getmontagesegmentedslices',no);
    if ismember(z,slicestoinclude)
    c = 1+mod(z-slicestoinclude(1),DATA.ViewPanelsMatrix{panel}(2));
    r = ceil((z-slicestoinclude(1)+1)/DATA.ViewPanelsMatrix{panel}(2));
    yofs = (c-1)*SET(no).YSize;
    xofs = (r-1)*SET(no).XSize;
    else
      yofs = nan;
      xofs = nan;
    end
  otherwise
    xofs = 0;
    yofs = 0;
end;

%------------------------------------------------------
function [cellx,celly] = calcoffsetcells(no,panel,type) %#ok<DEFNU>
%------------------------------------------------------
%Same as calculateoffset, but returns cells, used in 
%updatemodeldisplay.

global DATA SET NO
if nargin < 2
  panel = DATA.CurrentPanel;
end
if nargin==0
  no = NO;
end;
if nargin < 3
  type = DATA.ViewPanelsType{panel};
end

switch type
  case {'montage','montagerow','montagefit'}
    c = 1+mod((1:SET(no).ZSize)-1,DATA.ViewPanelsMatrix{panel}(2));
    r = ceil((1:SET(no).ZSize)/DATA.ViewPanelsMatrix{panel}(2));
    yofs = (c-1)*SET(no).YSize;
    xofs = (r-1)*SET(no).XSize;
    cellx = num2cell(xofs');
    celly = num2cell(yofs');
  case 'sax3'
    cellx = cell(SET(no).ZSize,1);
    celly = cell(SET(no).ZSize,1);
    for z = 1:SET(no).ZSize
      xofs = zeros(1,SET(no).TSize);
      yofs = zeros(1,SET(no).TSize);
      for t = 1:SET(no).TSize
        zind = find(SET(no).SAX3.slices(:,t) == z);
        if isempty(zind)
          xofs(t) = nan;
          yofs(t) = nan;
        else
          c = 1+mod(zind-1,DATA.ViewPanelsMatrix{panel}(2));
          r = ceil(zind/DATA.ViewPanelsMatrix{panel}(2));
          yofs(t) = (c-1)*SET(no).YSize;
          xofs(t) = (r-1)*SET(no).XSize;
        end
      end
      cellx{z} = xofs;
      celly{z} = yofs;
    end
  case 'montagesegmented'
    slicestoinclude = segment_main('getmontagesegmentedslices',no);
    cellx = cell(SET(no).ZSize,1);
    celly = cell(SET(no).ZSize,1);
    for z = 1:SET(no).ZSize
      if ismember(z,slicestoinclude)
        c = 1+mod(z-slicestoinclude(1),DATA.ViewPanelsMatrix{panel}(2));
        r = ceil((z-slicestoinclude(1)+1)/DATA.ViewPanelsMatrix{panel}(2));
        yofs = (c-1)*SET(no).YSize;
        xofs = (r-1)*SET(no).XSize;
      else
        yofs = nan;
        xofs = nan;
      end
      cellx{z} = xofs;
      celly{z} = yofs;
    end
  otherwise
    cellx = num2cell(zeros(SET(no).ZSize,1));
    celly = num2cell(zeros(SET(no).ZSize,1));
end;


%-------------------------------------------------
function [spacedist,timedist] = calcmmodedists(no) %#ok<DEFNU>
%-------------------------------------------------
%Calculate distances in space and time between mmode lines
global SET NO

if nargin < 1
  no = NO;
end
spacedist = sqrt(...
  ((SET(no).Mmode.Lx*SET(no).Mmode.M1-SET(no).Mmode.Lx*SET(no).Mmode.M2)*SET(no).ResolutionY).^2+...
  ((SET(no).Mmode.Ly*SET(no).Mmode.M1-SET(no).Mmode.Ly*SET(no).Mmode.M2)*SET(no).ResolutionX).^2);
  
if ~isrectilinear(SET(no).TimeVector)
  timedist = abs(...
    SET(no).TimeVector(round(SET(no).Mmode.T1))-...
    SET(no).TimeVector(round(SET(no).Mmode.T2)));
else
  timedist = abs(...
    SET(no).Mmode.T1-...
    SET(no).Mmode.T2);
  timedist = (timedist*SET(no).TIncr);
end
timedist = 1000*timedist;

%---------------------------
function calcdatasetpreview %#ok<DEFNU>
%---------------------------
%Calculate thumbnails. They are stored in the variable
%DATA.DATASETPREVIEW. Size of the thumbnails is given by
%DATA.GUISettings.ThumbnailSize. It is stored as a RGB image.

global DATA SET 

%Reserve memory
DATA.DATASETPREVIEW = repmat(uint8(0),...
  [DATA.GUISettings.ThumbnailSize*length(SET) DATA.GUISettings.ThumbnailSize 3]);

for noloop=1:length(SET);
  
  %Remap
  if isempty(SET(noloop).Colormap)
    tempim = remapuint8(...
      SET(noloop).IM(:,:,round(SET(noloop).TSize/2),round(SET(noloop).ZSize/2)),...
      noloop,returnmapping(noloop,true));
  else
    tempim = remapuint8(...
      SET(noloop).IM(:,:,round(SET(noloop).TSize/2),round(SET(noloop).ZSize/2)),...
      noloop);
  end
  
  % zero padding, elegantly done, no? :) /JU
  sz=size(tempim);
  tempim=padarray(tempim,round((length(tempim)-sz(1:2))/2));

  
  tempim = imresize(tempim,DATA.GUISettings.ThumbnailSize*[1 1],'bilinear');

  %Downsample
%  xind = round(linspace(1,size(tempim,1),64));
%  yind = round(linspace(1,size(tempim,2),64));
%  tempim = tempim(xind,:);
%  tempim = tempim(:,yind);
  
  %Store, vertically
  DATA.DATASETPREVIEW((noloop-1)*DATA.GUISettings.ThumbnailSize+(1:DATA.GUISettings.ThumbnailSize),:,:) = tempim;
end;

%----------------------------------------------------------
function [outim,slicestoinclude] = calcmontageviewim(no,matrix,segmentedonly,cmap,c,b,im,tfs,oneextraslice,slices) %#ok<DEFNU>
%----------------------------------------------------------
%Calculate view image for montage view.
global DATA SET

%Update number of rows and columns
if nargin < 2
  [rows,cols] = calcrowscols(no); %.cols,.rows
  matrix = [rows cols];
end

if nargin < 3
  segmentedonly = false;
end

if nargin < 7 || isempty(im)
  im = SET(no).IM;
end
if nargin < 8 || isempty(tfs)
  tfs = 1:SET(no).TSize;
end
if nargin < 9 || isempty(oneextraslice)
  oneextraslice = true;
end

%Decide which slices to view
if nargin < 10 || isempty(slices)
if segmentedonly %and currentslice
   slicestoinclude = segment('getmontagesegmentedslices',no);
%   slicestoincludeendo = find(findfunctions('findslicewithendo',no,tfs))';
%   slicestoincludeepi = find(findfunctions('findslicewithepi',no,tfs))';
%   slicestoincludervendo = find(findfunctions('findslicewithrvendo',no,tfs))';
%   slicestoinclude = unique([slicestoincludeendo slicestoincludeepi slicestoincludervendo]);
%   if min(slicestoinclude) > 1 && oneextraslice
%     slicestoinclude = [min(slicestoinclude)-1 slicestoinclude];
%   end
%   if max(slicestoinclude) < SET(no).ZSize && oneextraslice
%     slicestoinclude = [slicestoinclude max(slicestoinclude)+1];
%   end
else
  slicestoinclude = 1:SET(no).ZSize;
end
else
  slicestoinclude = slices;
end
%Add papillary visualization
stateandicon=segment('iconson','hidelv');
if not(stateandicon{1}) %not(isempty(SET(no).PapillaryIM)) && isequal(get(DATA.Handles.hidepapicon,'state'),'off')
  im(SET(no).PapillaryIM)=DATA.GUISettings.PapilarColor;
end

%Create space
if nargin >= 6 && ~isempty(cmap) && ~isempty(c) && ~isempty(b) 
  outim = repmat(uint8(0),[SET(no).XSize*matrix(1) SET(no).YSize*matrix(2) SET(no).TSize 3]);
elseif isempty(SET(no).Colormap)
  outim = repmat(uint8(0),[SET(no).XSize*matrix(1) SET(no).YSize*matrix(2) SET(no).TSize]);
else
  outim = repmat(uint8(0),[SET(no).XSize*matrix(1) SET(no).YSize*matrix(2) SET(no).TSize 3]);
end

for tloop=1:SET(no).TSize
  for i=1:numel(slicestoinclude)
    zloop = slicestoinclude(i);
    col = 1+mod(i-1,matrix(2));
    r = ceil(i/matrix(2));
    if nargin >= 6 && ~isempty(cmap) && ~isempty(c) && ~isempty(b) 
      outim(...
        (1+(r-1)*SET(no).XSize):(r*SET(no).XSize),...
        (1+(col-1)*SET(no).YSize):(col*SET(no).YSize),tloop,:) = remapuint8( ...
        im(:,:,tloop,zloop),no,cmap,c,b);
    else
      outim(...
        (1+(r-1)*SET(no).XSize):(r*SET(no).XSize),...
        (1+(col-1)*SET(no).YSize):(col*SET(no).YSize),tloop,:) = remapuint8( ...
        im(:,:,tloop,zloop),no);
    end
  end;
end;

%--------------------------------------------------------------------------
function [linex,liney] = calcmontageviewline(no,matrix,inputx,inputy,tf,segmentedonly,tfs,oneextraslice) %#ok<DEFNU>
%--------------------------------------------------------------------------
%Calculate view lines for montage view.
global SET

if nargin < 6
  segmentedonly = false;
end
if nargin < 7
  tfs = 1:SET(no).TSize;
end
if nargin < 8
  oneextraslice = true;
end

if isempty(inputx)
  linex=nan;
  liney=nan;
  return;
end

%Decide which slices to view , modified by Klas to include RV endo
if segmentedonly
  slicestoincludeendo = find(findfunctions('findslicewithendo',no,tfs))';
  slicestoincludeepi = find(findfunctions('findslicewithepi',no,tfs))';
  slicestoincludervendo = find(findfunctions('findslicewithrvendo',no,tfs))';
  
  slicestoinclude = unique([slicestoincludeendo slicestoincludeepi slicestoincludervendo]);
  if min(slicestoinclude) > 1 && oneextraslice
    slicestoinclude = [min(slicestoinclude)-1 slicestoinclude];
  end
  if max(slicestoinclude) < SET(no).ZSize && oneextraslice
    slicestoinclude = [slicestoinclude max(slicestoinclude)+1];
  end
else
  slicestoinclude = 1:SET(no).ZSize;
end

sizeline = size(inputx);

%Create space
linex = nan((sizeline(1)+1)*length(slicestoinclude),1);
liney = linex;

for i=1:numel(slicestoinclude)
  zloop = slicestoinclude(i);
  c = 1+mod(i-1,matrix(2));
  r = ceil(i/matrix(2));
  s = (i-1)*(sizeline(1)+1)+1;
  e = (i-1)*(sizeline(1)+1)+sizeline(1);
  linex(s:e) = inputx(:,tf,zloop)+(1+(r-1)*SET(no).XSize);
  liney(s:e) = inputy(:,tf,zloop)+(1+(c-1)*SET(no).YSize);
end

%-------------------------------------------------------
function map = returnmapping(no,forcetruecolor)
%-------------------------------------------------------
%If SET(no).Colormap exists:
%  Returns truecolor colormap for this.
%If not:
%  Returns grayscale 1-D indexcolor map.
%Unless:
%  forcetruecolor is true, then 3-D truecolor grayscale is returned,
%  regardless of SET(no).Colormap
%  See REMAPUINT8 below.

persistent grayxi
global DATA SET

if isempty(grayxi)
  grayxi= (0:DATA.GUISettings.ColorMapSize-1)'/(DATA.GUISettings.ColorMapSize-1);
end

if (nargin==2) && forcetruecolor
  map = [grayxi grayxi grayxi];
else
  if isempty(SET(no).Colormap)
    map = grayxi;
  else
    map = SET(no).Colormap;
  end
end

%--------------------------------
function z = remapuint8(varargin)
%--------------------------------
%Same as remapuint8viewim, but squeezes result to make sure it can be
%viewed other than in the main window.
z = remapuint8viewim(varargin{:});
if length(size(z)) > 3
  z = squeeze(z);
end

%-----------------------------------------------
function z = remapuint8viewim(im,no,tempmap,c,b)
%-----------------------------------------------
% If SET(no).Colormap exists:
%   Remap from image data to true color using color lookup.
% If not:
%   Remap according to indexed grayscape
% Unless:
%   An external colormap has been supplied, then use that.
%   This is used to force truecolor grayscale, by passing 
%   true to RETURNMAPPING.
% c,b are contrast/brightness settings, used by contrast_Callback when
% doing realtime update during mouse motion.
%
%  See RETURNMAPPING above

global SET NO

if nargin == 1
  no = NO;
end

if nargin >= 3
  %tempmap is only used for realtime updating of contrast/brightness tool.
  cmap = tempmap;
else
  % This uses SET(no).Colormap if it exists, grayscale indexcolor if not.
  cmap = returnmapping(no);
end

cmap = uint8(255*cmap);
cmax = length(cmap);
cmin = 1;

if (nargin<5)
  c = SET(no).IntensityMapping.Contrast;
  b = SET(no).IntensityMapping.Brightness;
  if isempty(c)
    c = SET(SET(no).Parent).IntensityMapping.Contrast;
    b = SET(SET(no).Parent).IntensityMapping.Brightness;    
  end
end

if ~isempty(SET(no).IntensityMapping.Compression)
  d = SET(no).IntensityMapping.Compression;
else
  d=1;
end

if isempty(d) && ~isempty(SET(SET(no).Parent))
  if isempty(SET(SET(no).Parent).IntensityMapping.Compression)
    d = SET(SET(no).Parent).IntensityMapping.Compression;
  else
    d=1;
  end
end

if isempty(d)
    d=1;
end

if d<0.1
  d=0.1;
  SET(no).IntensityMapping.Compression=0.1;
end

if d>20
  d=0.1;
  SET(no).IntensityMapping.Compression=0.1;
end

sz = [size(im,1) size(im,2) size(im,3) size(im,4)];
if isa(im,'single')||isa(im,'double')
  %im = c*im(:)+(b-0.5);
%   try
if d~=1
  try
    im = c.*nthroot(im(:),d)+b-0.5;
  catch
    im = c.*im(:)+b-0.5;
  end
else
  im = c.*im(:)+b-0.5;
end%c*255.*(im(:)./255+b-0.5).^(1./d);
%   catch
%     %first lift im so that it is between 0 and 1
%     im=abs(min(im(:)))+im;
%     im = c.*nthroot(im(:),d)+b-0.5;%c*255.*(im(:)./255+b-0.5).^(1./d);
%   end
  z = max(min(round(cmax*im),cmax),cmin);
  
elseif isa(im,'int16'),
  if ~isfield(SET(no),'minValue')||isempty(SET(no).minValue)||...
      ~isfield(SET(no),'maxValue')||isempty(SET(no).maxValue)
    SET(no).minValue = single(min(SET(no).IM(:)));
    SET(no).maxValue = single(max(SET(no).IM(:)));
  end
   mi = SET(no).minValue;
   ma = SET(no).maxValue;
   
   im = single(im(:));
   %im=c*im+(b-0.5);%TODO: should scale
   
   z = floor((im-mi)/(ma-mi)*cmax);
   if d~=1
     z = uint8(c.*nthroot(z(:),d)+b-0.5)+1;%uint8(c*z+(b-0.5)*256)+1; %
   else
     z = uint8(c.*z(:)+b-0.5)+1;
   end;
   %z = max(min(round(c*z+(b-0.5)),cmax),cmin);   
else
  myfailed(dprintf('Segment does not (yet) support type:%s',class(im)));
  return;
end

if (size(cmap,2) == 3)
    % Extract r,g,b components
    z3 = repmat(uint8(0),[size(z) 3]);
    z3(:,1) = cmap(z,1);
    z3(:,2) = cmap(z,2);
    z3(:,3) = cmap(z,3);
    z = reshape(z3,[sz 3]);
else
    %Index color
    z3 = repmat(uint8(0),size(z));
    z3(:) = cmap(z(:),1);
    z = reshape(z3,sz);
end

%-------------------------------------------------------------------
function [xnew,ynew] = resampleclosedcurve(x,y,n,method,distributed)
%-------------------------------------------------------------------
%General helper fcn to resample a closed curve.
if nargin<4
	method='linear';%used to be linear*
  distributed = 'same';
end
if nargin < 5
  distributed = 'same';
end

switch distributed
  
  case 'same'
    xnew = interp1(x,linspace(1,size(x,1),n),method);
    ynew = interp1(y,linspace(1,size(y,1),n),method);
    
  case 'evenly'
    %Remve duplicate points  
    d = [0;cumsum(sqrt(diff(x).^2+diff(y).^2))];
    ind = [1;diff(d)]>1e-3;
    uniquex = x(ind);
    uniquey = y(ind);
    if sqrt((uniquex(1)-uniquex(end)).^2+(uniquey(1)-uniquey(end)).^2)>1e-3
      uniquex = uniquex(1:end-1);
      uniquey = uniquey(1:end-1);
    end
    %Find center,angle, sort on angle
    mx = mean(uniquex);
    my = mean(uniquey);
    rad = sqrt((uniquex-mx).^2+(uniquey-my).^2);
    ang = angle(complex(uniquex-mx,uniquey-my));
    %ensure unique values in angle to be used in interpolation
    [ang,anguniqueindex] = unique(ang(:));
    rad = rad(:);
    rad = rad(anguniqueindex);
    %resample and interpolate radius to get points evenly
    %distributed by angle
    aaa = [ang(:)+2*pi ; ang(:) ; ang(:)-2*pi];
    rrr = [rad(:) ; rad(:) ; rad(:)];
    resampleangle = linspace(pi/2,-3*pi/2,n+1)';
    newr = interp1(aaa,rrr,resampleangle(1:end-1),'linear');
    z = newr.*exp(1i*resampleangle(1:end-1));
    xnew = real(z) + mx;
    ynew = imag(z) + my;
end


%-------------------------------------------------------------
function [xnew,ynew] = resamplemodel(x,y,n,method,distributed)
%-------------------------------------------------------------
%Resample the a stack of contours (typically segmentation)
%to different number of points (n).

if nargin<4
	method='linear';%used to be linear*
  distributed = 'same';
end
if nargin < 5 
  distributed = 'same';
end

xnew = nan(n,size(x,2),size(x,3));
ynew = xnew;

for zloop=1:size(x,3)  
  for tloop=1:size(x,2)
    if not(isnan(x(1,tloop,zloop)))
      [tempx,tempy] = resampleclosedcurve(...
        x(:,tloop,zloop),...
        y(:,tloop,zloop),n,method,distributed);
      xnew(:,tloop,zloop) = tempx;
      ynew(:,tloop,zloop) = tempy;
    end
  end;
end;

 
%--------------------------------------------------------------------------
function [meanint,defectextent,varargout] = calcintensityanddefect(im,...
  tf,numpoints,numsectors,nprofiles,pos, ...
  no,endox,endoy,epix,epiy,sz,resolution,defect,numwidth,timeresolved) %#ok<DEFNU>
%--------------------------------------------------------------------------
%Calculates the mean intensity within sectors as a preperation to 
%generate a bullseye plot.
%
%-im                                     image
%-tf                                      timeframe
%-numpoint                        numpoints to evaluate in
%- numsectors:                   number of sectors
%- nprofiles:                         number of profiles
%- pos:                                  slices to use
%- endox,endoy,epix,epiy: myocardial borders
%sz:                                        size of image stack
%resolution:                         resolution of image stack
%defect:                                (EH: I do not know what it is, works if empty)
%numwidth:                        specify if it should calculate several sectors across 
%                                           wallthickness. If a scalar it
%                                           specifices the number of
%                                           layers in the myocardium. If a
%                                           two element vector assumes
%                                           endo and epi percentages cut.
%timeresolved:                   true if timeresolvde

%This function is in serious need of refactory as it is very complicated
%with multiple input arguments.  EH:

if nargin < 13
  myfailed('Too few input arguments.');
  return;
end
if nargin < 14
  defect = [];
end
if nargin < 15
  numwidth = 1;
end
if nargin < 16
  timeresolved = true;
end

varargout = {};
numslices = length(pos);
meanint = zeros(numsectors,numslices);

%Added EH: If numwidth is a vector then it is interpreted as 
%[perencentendo percent epi].
if length(numwidth)>1
  numlayers = length(numwidth)-1;
else
  numlayers = numwidth;
end;

%Check if we should return the mask and the sectors
if nargout>2
  storemask = true;
  maskcell = cell(numslices,numsectors,numlayers);
else
  storemask = false;
end;
if nargout>3
  storesector = true;
  sectorcell = cell(numslices,numsectors,numlayers,2);
else
  storesector = false;
end;

%Upsample model
if not(nprofiles == size(endox,1))
  [endox,endoy] = resamplemodel(endox,endoy,nprofiles);
  [epix,epiy] = resamplemodel(epix,epiy,nprofiles);
end

mask = zeros(sz(1),sz(2));
defectextent = zeros(numsectors,numslices);
for sloop = 1:numslices
  if not(isnan(epix(1,tf,pos(sloop))))

    %Find sectors
    [~,~,episectors] = findmeaninsectorslice('epi',numpoints,tf,pos(sloop),numsectors,no, ...
      endox(:,tf,pos(sloop)),endoy(:,tf,pos(sloop)),epix(:,tf,pos(sloop)),epiy(:,tf,pos(sloop)));
    if not(isnan(endox(1,tf,pos(sloop))))
    [~,~,endosectors] = findmeaninsectorslice('endo',numpoints,tf,pos(sloop),numsectors,no, ...
      endox(:,tf,pos(sloop)),endoy(:,tf,pos(sloop)),epix(:,tf,pos(sloop)),epiy(:,tf,pos(sloop)));
    else
      endosectors = episectors;
    end
    
    %calulate wallthickness for whole slice
    xdiffslice = abs(endox(:,tf,pos(sloop))-epix(:,tf,pos(sloop)));
    ydiffslice = abs(endoy(:,tf,pos(sloop))-epiy(:,tf,pos(sloop)));
    wallthicknessinslice = mean(sqrt((xdiffslice*resolution(1)).^2+(ydiffslice*resolution(2)).^2));

    N = nprofiles;
    for loop=1:numsectors
      if episectors(loop)<episectors(loop+1)
        tempepix = epix(episectors(loop):episectors(loop+1),tf,pos(sloop));
        tempepiy = epiy(episectors(loop):episectors(loop+1),tf,pos(sloop));
      elseif numsectors == 1
        tempepix = epix(:,tf,pos(sloop));
        tempepiy = epiy(:,tf,pos(sloop));
      elseif episectors(loop)==episectors(loop+1)
        if episectors(loop) == 1
          startsector = N;
          endsector = episectors(loop)+1;
        elseif episectors(loop) == N
          startsector = episectors(loop)-1;
          endsector = 1;
        else
          startsector = episectors(loop)-1;
          endsector = episectors(loop)+1;
        end
        tempepix = epix(episectors(loop),tf,pos(sloop))+ ...
          [0.5*(epix(startsector,tf,pos(sloop))-epix(episectors(loop),tf,pos(sloop))); ...
          0.5*(epix(endsector,tf,pos(sloop))-epix(episectors(loop),tf,pos(sloop)))];
        tempepiy = epiy(episectors(loop),tf,pos(sloop))+ ...
          [0.5*(epiy(startsector,tf,pos(sloop))-epiy(episectors(loop),tf,pos(sloop))); ...
          0.5*(epiy(endsector,tf,pos(sloop))-epiy(episectors(loop),tf,pos(sloop)))];
      else
        tempepix = [...
          epix(episectors(loop):N,tf,pos(sloop)) ; ...
          epix(1:episectors(loop+1),tf,pos(sloop))];
        tempepiy = [...
          epiy(episectors(loop):N,tf,pos(sloop)) ; ...
          epiy(1:episectors(loop+1),tf,pos(sloop))];
      end;
      
      if not(isnan(endox(1,tf,pos(sloop))))
        if endosectors(loop)<endosectors(loop+1)
          tempendox = endox(endosectors(loop):endosectors(loop+1),tf,pos(sloop));
          tempendoy = endoy(endosectors(loop):endosectors(loop+1),tf,pos(sloop));
        elseif numsectors == 1
          tempendox = endox(:,tf,pos(sloop));
          tempendoy = endoy(:,tf,pos(sloop));
        elseif endosectors(loop)==endosectors(loop+1)
          if endosectors(loop) == 1
            startsector = N;
            endsector = endosectors(loop)+1;
          elseif endosectors(loop) == N
            startsector = endosectors(loop)-1;
            endsector = 1;
          else
            startsector = endosectors(loop)-1;
            endsector = endosectors(loop)+1;
          end
          tempendox = endox(endosectors(loop),tf,pos(sloop))+ ...
            [0.5*(endox(startsector,tf,pos(sloop))-endox(endosectors(loop),tf,pos(sloop))); ...
            0.5*(endox(endsector,tf,pos(sloop))-endox(endosectors(loop),tf,pos(sloop)))];
          tempendoy = endoy(endosectors(loop),tf,pos(sloop))+ ...
            [0.5*(endoy(startsector,tf,pos(sloop))-endoy(endosectors(loop),tf,pos(sloop))); ...
            0.5*(endoy(endsector,tf,pos(sloop))-endoy(endosectors(loop),tf,pos(sloop)))];
        else
          tempendox = [...
            endox(endosectors(loop):N,tf,pos(sloop)) ; ...
            endox(1:endosectors(loop+1),tf,pos(sloop))];
          tempendoy = [...
            endoy(endosectors(loop):N,tf,pos(sloop)) ; ...
            endoy(1:endosectors(loop+1),tf,pos(sloop))];
        end;
      else
        tempendox = repmat(mean(epix(:,tf,pos(sloop))),[2 1]);
        tempendoy = repmat(mean(epiy(:,tf,pos(sloop))),[2 1]);
      end
      
      if (length(numwidth) == 1) && isequal(numwidth,1) %standard
        %Create mask image
        mask = roipoly(mask,...
          [tempendoy;flipud(tempepiy);tempendoy(1)],...
          [tempendox;flipud(tempepix);tempendox(1)]);
        if nansum(mask(:)) == 0
          mask(round(mean([tempendox;tempepix])),round(mean([tempendoy;tempepiy]))) = 1;
        end
        %the same number of points along the lines
        if length(tempepix) ~= length(tempendox)
          xi = linspace(0,1,length(tempendox));
          xx = linspace(0,1,length(tempepix));
          tempendox = interp1(xi,tempendox,xx)';
          tempendoy = interp1(xi,tempendoy,xx)';
        end
        %Store the mask
        if storemask
          maskcell{sloop,loop} = mask;
        end;
                
        %calculate the wallthickness for the current mask
        xdiff = abs(tempendox-tempepix); %min(abs(repmat(tempendox,[1 length(tempepix)])-repmat(tempepix',[length(tempendox) 1])),[],2);
        ydiff = abs(tempendoy-tempepiy); %min(abs(repmat(tempendoy,[1 length(tempepiy)])-repmat(tempepiy',[length(tempendoy) 1])),[],2);
        tempwallthickness = mean(sqrt((xdiff*resolution(1)).^2+(ydiff*resolution(2)).^2));
        if sum(sum(mask)) == 0 && tempwallthickness/wallthicknessinslice > 0.15
          %create a mask of size 1 pixel if the mask is empty and the wall
          %thickness is over 15% of mean wall thickness in this slice
          mask(round(mean([tempendoy;flipud(tempepiy)])),round(mean([tempendox;flipud(tempepix)]))) = 1;
        end
        %calculate mean intensity within the sector
        if nansum(mask(:)) == 0 || tempwallthickness/wallthicknessinslice < 0.15
          meanint(loop,sloop) = NaN;
        else
          if length(size(im)) == 3 && timeresolved == 0
            meanint(loop,sloop) = nansum(nansum(mask.*im(:,:,pos(sloop))))/(nansum(nansum(mask)));
          else
            meanint(loop,sloop) = nansum(nansum(mask.*im(:,:,tf,pos(sloop))))/(nansum(nansum(mask)));
          end
        end
        if length(size(defect)) == 3 && timeresolved == 0
          defectextent(loop,sloop) = 100*sum(sum(mask.*defect(:,:,pos(sloop))))/(sum(sum(mask))+eps);
        elseif length(size(defect)) == 4
          defectextent(loop,sloop) = 100*sum(sum(mask.*defect(:,:,tf,pos(sloop))))/(sum(sum(mask))+eps);
        else
          defectextent = [];
        end
        
      else %to divide the myocardium into several parts in the width
        
        %%% numwidth is either a scalar >1 or a two element vector
        
        if length(numwidth)==1
          %Is a scalar
          PercentFromEndo = (0:numwidth)/numwidth; 
        else          
          PercentFromEndo = numwidth;
        end;
        
        %the same number of points along the lines
        if length(tempepix) ~= length(tempendox)
          xi = linspace(0,1,length(tempendox));
          xx = linspace(0,1,length(tempepix));
          tempendox = interp1(xi,tempendox,xx)';
          tempendoy = interp1(xi,tempendoy,xx)';
        end
        for wloop = 1:(length(PercentFromEndo)-1) %was numwidth
          %Create mask image
          xin = tempepix+(tempendox-tempepix)*(1-PercentFromEndo(wloop));
          yin = tempepiy+(tempendoy-tempepiy)*(1-PercentFromEndo(wloop));
          xout = tempepix+(tempendox-tempepix)*(1-PercentFromEndo(wloop+1));
          yout = tempepiy+(tempendoy-tempepiy)*(1-PercentFromEndo(wloop+1));
          sectory = [yin;flipud(yout);yin(1)];
          sectorx = [xin;flipud(xout);xin(1)];
          mask = roipoly(mask,sectory,sectorx);
          if nansum(mask(:)) == 0
            sectory = round(mean([tempendox;tempepix]));
            sectorx = round(mean([tempendoy;tempepiy]));
            mask(sectory,sectorx) = 1;
          end
          %Store the mask
          if storemask
            maskcell{sloop,loop,wloop} = mask;
          end
          if storesector
            sectorcell{sloop,loop,wloop,1} = sectorx;
            sectorcell{sloop,loop,wloop,2} = sectory;
          end
          
          %calculate the wallthickness for the current mask
          xdiff = abs(tempendox-tempepix); %min(abs(repmat(tempendox,[1 length(tempepix)])-repmat(tempepix',[length(tempendox) 1])),[],2);
          ydiff = abs(tempendoy-tempepiy); %min(abs(repmat(tempendoy,[1 length(tempepiy)])-repmat(tempepiy',[length(tempendoy) 1])),[],2);
          tempwallthickness = mean(sqrt((xdiff*resolution(1)).^2+(ydiff*resolution(2)).^2));
          if sum(sum(mask)) == 0 && tempwallthickness/wallthicknessinslice > 0.15
            %create a mask of size 1 pixel if the mask is empty and the wall thickness is over 2 mm
            mask(round(mean([tempendoy;flipud(tempepiy)])),round(mean([tempendox;flipud(tempepix)]))) = 1;
          end
          %calculate mean intensity within the sector
          if nansum(mask(:)) == 0 || tempwallthickness/wallthicknessinslice < 0.15
            meanint(loop,sloop,wloop) = NaN;
          else
            if length(size(im)) == 3 && timeresolved == 0
              meanint(loop,sloop,wloop) = nansum(nansum(mask.*im(:,:,pos(sloop))))/(nansum(nansum(mask)));
            else
              meanint(loop,sloop,wloop) = nansum(nansum(mask.*im(:,:,tf,pos(sloop))))/(nansum(nansum(mask)));
            end
          end
          if length(size(defect)) == 3 && timeresolved == 0
            defectextent(loop,sloop,wloop) = 100*sum(sum(mask.*defect(:,:,pos(sloop))))/(sum(sum(mask))+eps);
          elseif length(size(defect)) == 4
            defectextent(loop,sloop,wloop) = 100*sum(sum(mask.*defect(:,:,tf,pos(sloop))))/(sum(sum(mask))+eps);
          else
            defectextent = [];
          end
        end
      end 
    end
  end
end

%If wanted return the extra output
if storemask
  varargout{1} = maskcell;
end;
if storesector
  varargout{2} = sectorcell;
end;

%----------------------------------------------------------------------
function [res,varargout] = findmeaninsector(type,inp,pos,numsectors,no,tf)
%----------------------------------------------------------------------
%Find indicies of points along the contour that corresponds to which
%sector. Then also calculated mean values of the input inp. 

global DATA SET NO

if nargin<4
  numsectors=6;
end;

if nargin<5
  no = NO;
end;

if nargin<6
  tf = 1:SET(no).TSize;
end

if strcmp(SET(no).ImageViewPlane,'Short-axis')
  distributed = 'evenly';
else
  distributed = 'same';  
end

%Parse input
numtf = length(tf);
numpoints = size(inp,1);
isendo = isequal(type,'endo');
nslices = length(pos);

%Upsample model
if not(numpoints==DATA.NumPoints)
  if ~isempty(SET(no).EndoX)
  [endox,endoy] = resamplemodel(SET(no).EndoX(:,tf,:),SET(no).EndoY(:,tf,:),numpoints,'linear',distributed);
  else
    endox = nan(numpoints,numtf,SET(no).ZSize);
    endoy = endox;
  end;
  
  if ~isempty(SET(no).EpiX)
    [epix,epiy] = resamplemodel(SET(no).EpiX(:,tf,:),SET(no).EpiY(:,tf,:),numpoints,'linear',distributed);
  else
    epix = nan(numpoints,numtf,SET(no).ZSize);
    epiy = epix;    
  end;
  
else
  if ~isempty(SET(no).EndoX)
    endox = SET(no).EndoX(:,tf,:);
    endoy = SET(no).EndoY(:,tf,:);
  else
    endox = nan(numpoints,numtf,SET(no).ZSize);
    endoy = endox;    
  end;
  if ~isempty(SET(no).EpiX)
    epix = SET(no).EpiX(:,tf,:);
    epiy = SET(no).EpiY(:,tf,:);
  else
    epix = nan(numpoints,numtf,SET(no).ZSize);
    epiy = epix;        
  end;
end;

meanx = nan(numtf,nslices);
meany = nan(numtf,nslices);
ang = nan([numpoints-1 numtf nslices]);

if isendo
  %Endocardium
  for sloop=1:nslices
    for tloop=1:numtf

      %Find center
      x = endox(1:(end-1),tloop,pos(sloop));
      y = endoy(1:(end-1),tloop,pos(sloop));
      if isnan(epix(1,tloop,pos(sloop)))||SET(no).EndoCenter
        meanx(tloop,pos(sloop)) = mean(x(:));
        meany(tloop,pos(sloop)) = mean(y(:));
      else
        meanx(tloop,pos(sloop)) = mean(epix(:,tloop,pos(sloop)));
        meany(tloop,pos(sloop)) = mean(epiy(:,tloop,pos(sloop)));
      end;

      %Find sectors
      ang(:,tloop,sloop) = angle(complex(...
        endoy(1:(end-1),tloop,pos(sloop))-meany(tloop,pos(sloop)),...
        endox(1:(end-1),tloop,pos(sloop))-meanx(tloop,pos(sloop))));
    end;
  end;
else
  %Epicardium
  for sloop=1:nslices
    for tloop=1:numtf
      %Find center
      if SET(no).EndoCenter
        x = endox(1:(end-1),tloop,pos(sloop));
        y = endoy(1:(end-1),tloop,pos(sloop));        
      else
        x = epix(1:(end-1),tloop,pos(sloop));
        y = epiy(1:(end-1),tloop,pos(sloop));
      end;
      meanx(tloop,pos(sloop)) = mean(x(:));
      meany(tloop,pos(sloop)) = mean(y(:));

      %Find sectors
      ang(:,tloop,sloop) = angle(complex(...
        epiy(1:(end-1),tloop,pos(sloop))-meany(tloop,pos(sloop)),...
        epix(1:(end-1),tloop,pos(sloop))-meanx(tloop,pos(sloop))));
    end;
  end;
end;

ang = ang*180/pi;
wantangles = linspace(-180,180,numsectors+1);
wantangles = wantangles+SET(no).SectorRotation;

ind = (wantangles>180);
wantangles(ind) = wantangles(ind)-360;

ind = (wantangles<-180);
wantangles(ind) = wantangles(ind)+360;

sectors = zeros(numsectors+1,nslices,numtf);
for tloop=1:numtf
  for sloop=1:nslices
    for loop=1:numsectors
      [~,inda] = min(abs(ang(:,tloop,sloop)-wantangles(loop)));
      sectors(loop,sloop,tloop) = inda;
    end;
  end;
end;
sectors(numsectors+1,:,:) = sectors(1,:,:);

%Calculate mean over each sector
res = zeros(numsectors,nslices,numtf);
maxres = res;
warnstat = warning;
warning off; %#ok<WNOFF>
for tloop=1:numtf
  for zloop=1:nslices
    for loop=1:numsectors
      if sectors(loop,zloop,tloop)<sectors(loop+1,zloop,tloop)
        idx = sectors(loop,zloop,tloop):sectors(loop+1,zloop,tloop);
      else
        idx = [...
          sectors(loop,zloop,tloop):(numpoints-1) ...
          1:sectors(loop+1,zloop,tloop)];
      end;
      res(loop,zloop,tloop) = mean(inp(idx,tloop,pos(zloop)));
      maxres(loop,zloop,tloop) = max(inp(idx,tloop,pos(zloop)));
    end;
  end;
end;
warning(warnstat);

varargout = cell(1,nargout-1);
if nargout>=2
  varargout{1} = meanx;
end;
if nargout>=3
  varargout{2} = meany;
end;
if nargout>=4
  varargout{3} = sectors;
end;
if nargout>=5
  varargout{4} = maxres;
end;

%--------------------------------------------------------------------------
function [varargout] = findmeaninsectorslice(type,numpoints,tf,slice,...
  numsectors,no,endox,endoy,epix,epiy,rotation)
%--------------------------------------------------------------------------
%Find indicies of points along the contour that corresponds to which
%sector. Used when analysing the myocardium of short axis slices. Operates
%on given slice.
global DATA SET NO

if nargin<5
  numsectors = 6;
end;

if nargin<6
  no = NO;
end;

if nargin < 7 || (nargin==11 && isempty(endox))
  if ~isempty(SET(no).EndoX)
    endox = SET(no).EndoX(:,tf,slice);
    endoy = SET(no).EndoY(:,tf,slice);
  else
    endox = nan(numpoints,length(tf),length(slice));
    endoy = endox;
  end
  if ~isempty(SET(no).EpiX)
    epix = SET(no).EpiX(:,tf,slice);
    epiy = SET(no).EpiY(:,tf,slice);
  else
    epix = nan(numpoints,length(tf),length(slice));
    epiy = epix;
  end
end

if nargin < 11
  rotation = SET(no).SectorRotation;
end

if isequal(type,'endo')
  isendo=1;
else
  isendo=0;
end;

%Upsample model
if not(numpoints==DATA.NumPoints)
  if ~isempty(endox)
    [endox,endoy] = resamplemodel(endox,endoy,numpoints);
  else
    endox = nan(numpoints,length(tf),length(slice));
    endoy = endox;
  end;
  if ~isempty(epix)
    [epix,epiy] = resamplemodel(epix,epiy,numpoints);
  else
    epix = nan(numpoints,length(tf),length(slice));
    epiy = epix;    
  end;
end

%Find center
if SET(no).EndoCenter && ~isempty(endox) && ~isnan(endox(1))
  x = mean(endox);
  y = mean(endoy);
else
  x = mean(epix);
  y = mean(epiy);
end;

meanx = mean(x(:));
meany = mean(y(:));

%Find sectors
if isendo
  ang = angle(complex(...
    endoy(1:(end-1))-meany,...
    endox(1:(end-1))-meanx));
else
  ang = angle(complex(...
    epiy(1:(end-1))-meany,...
    epix(1:(end-1))-meanx));
end;
ang = ang*180/pi;
wantangles = linspace(-180,180,numsectors+1);
wantangles = wantangles+rotation;

sectors = zeros(numsectors+1,1);
for loop=1:(numsectors+1)
  [trash1,inda] = min(abs(ang(:)-wantangles(loop)));
  [trash2,temp] = min(abs(ang(:)-wantangles(loop)-360));
  if trash2<trash1
    inda = temp;
    trash1 = trash2;
  end;
  [trash2,temp] = min(abs(ang(:)-wantangles(loop)+360));
  if trash2<trash1
    inda = temp;
  end;
  sectors(loop) = inda;
end;

varargout = cell(1,nargout-1);
if nargout>=1
  varargout{1} = meanx;
end;
if nargout>=2
  varargout{2} = meany;
end;
if nargout>=3
  varargout{3} = sectors;
end;

%-----------------------------------------------------------------------
function [pos,xdir,ydir,zdir] = getimageposandorientation(no,slice,type)
%-----------------------------------------------------------------------
%Return image orientation for current slice or second input argument.
%If two outputvariables, then zdir is also outputed.
global SET

if nargin < 3
  type = 'one';
end
if nargin<2
  slice = SET(no).CurrentSlice;
end;

if SET(no).Rotated
  %Rotated
  pos = SET(no).ImagePosition;
  rotdir = SET(no).ImageOrientation(4:6);
  r1dir = SET(no).ImageOrientation(1:3);
  r2dir = cross(rotdir,r1dir);

  %Calculate radius to center
  radius = SET(no).RotationCenter*SET(no).ResolutionY;
    
  %Calculate center of rotation
  center = pos+radius*r1dir;
    
  alpha = (pi/(SET(no).ZSize))*(slice-1);
  pos = center+(...
    radius*cos(alpha)*r1dir-...
    radius*sin(alpha)*r2dir);

  %Find orientation
  xdir = rotdir;
  ydir = -(pos-center); %new position - center of rotation
  ydir = ydir./sqrt(sum(ydir.^2));
  zdir = cross(ydir,xdir);
else
  %Normal cartesian coordinates
  xdir = SET(no).ImageOrientation(4:6);
  ydir = SET(no).ImageOrientation(1:3);  
  zdir = cross(ydir,xdir);

  
  if strcmp(type,'hla')
    temp = xdir;
    xdir = -zdir;
    zdir = -temp;
    pos = SET(no).ImagePosition-...
      zdir*SET(no).ResolutionX*(slice-1);
  elseif strcmp(type,'vla')
    temp = ydir;
    ydir = xdir;
    xdir = -zdir;
    zdir = temp;
    pos = SET(no).ImagePosition+...
      zdir*SET(no).ResolutionY*(slice-1);
  elseif strcmp(type,'gla')
    ang = SET(no).GLA.angle;
    temp = -xdir*cos(ang)+ydir*sin(ang);
    temp2 = ydir*cos(ang)+xdir*sin(ang);
    xdir = -zdir;
    ydir = temp2;%*sign(cos(ang));
    zdir = temp;
    pos = SET(no).ImagePosition-...
      zdir*SET(no).ResolutionX*(slice-1)*cos(ang)+ ...
      zdir*SET(no).ResolutionY*(slice-1)*sin(ang);
    pos = xyz2rlapfh(no,SET(no).GLA.x0,SET(no).GLA.y0,1);
  else
    %Take a slice depending offset.
    pos = SET(no).ImagePosition-...
      zdir*(SET(no).SliceThickness+SET(no).SliceGap)*(slice-1);
    
%     %for multiple slices (future improvement)
%     pos = repmat(SET(no).ImagePosition,length(slice),1)-...
%       (diag(zdir*(SET(no).SliceThickness+SET(no).SliceGap))*repmat(slice-1,3,1))';
    
  end
end;

%------------------------
function z = econvn(im,f)
%------------------------
%Function to filter image, uses fastest available convolver
switch ndims(im)
  case 4
    z = im;
    if size(im,3)~=1
      for loop=1:size(im,4);
        z(:,:,:,loop) = econv3(im(:,:,:,loop),f);
      end;
    else
      for loop=1:size(im,4);
        z(:,:,1,loop) = conv2(im(:,:,:,loop),f,'same');
      end;
    end;
  case 3
    z = econv3(im,f);
  case 2
    z = conv2(im,f,'same');
  otherwise
    error('Number of dimensions not supported.');
end;

%---------------------------------------
function [pos] = rlapfh2xyz(no,rl,ap,fh)
%---------------------------------------
%Convert from RL,AP,FH coordinates to Segment internal coordinate system.
global DATA SET

if SET(no).Rotated
  myfailed('Not yet implemented for rotated image stacks.',DATA.GUI.Segment);
else
  xdir = SET(no).ImageOrientation(4:6)';
  ydir = SET(no).ImageOrientation(1:3)';  
  zdir = cross(...
    SET(no).ImageOrientation(1:3),...
    SET(no).ImageOrientation(4:6))';    
  
  rl = rl(:)'-SET(no).ImagePosition(1); %Translate corner of box to 0,0,0
  ap = ap(:)'-SET(no).ImagePosition(2);
  fh = fh(:)'-SET(no).ImagePosition(3);
  
  pos = [xdir ydir -zdir]\[rl;ap;fh];
  pos(1,:) = pos(1,:)/SET(no).ResolutionX+1;
  pos(2,:) = pos(2,:)/SET(no).ResolutionY+1;
  pos(3,:) = pos(3,:)/(SET(no).SliceGap+SET(no).SliceThickness)+1;  
  
end;

%------------------------------------
function [pos] = xyz2rlapfh(no,x,y,z)
%------------------------------------
%Converts from segment coordinate system to RL,AP,FH coordinate system.
global DATA SET

if SET(no).Rotated
  pos = [];
  myfailed('Not yet implemented for rotated image stacks.',DATA.GUI.Segment);
else
  zdir = cross(...
    SET(no).ImageOrientation(1:3),...
    SET(no).ImageOrientation(4:6));  
  x = (x(:)-1)*SET(no).ResolutionX;
  y = (y(:)-1)*SET(no).ResolutionY;
  z = (z(:)-1)*(SET(no).SliceThickness+SET(no).SliceGap);
  pos = repmat(SET(no).ImagePosition,length(x),1)+...
    repmat(SET(no).ImageOrientation(4:6),length(x),1).*repmat(x,1,3)+...
    repmat(SET(no).ImageOrientation(1:3),length(x),1).*repmat(y,1,3)-... %was minus
    repmat(zdir,length(x),1).*repmat(z,1,3);    
end;

%-------------------------------------
function [impos,ormat] = calcormat(no) %#ok<DEFNU>
%-------------------------------------
%Calculate orientation matrix for stack number no
%Can be used for coordinate transformations:
%RLAPFH = impos + ormat * (XYZ - 1)
%XYZ = ormat \ (RLAPFH - impos) + 1
global SET

impos = SET(no).ImagePosition';
ormat = [SET(no).ResolutionX*SET(no).ImageOrientation(4:6)' ...
  SET(no).ResolutionY*SET(no).ImageOrientation(1:3)' ...
  (SET(no).SliceGap+SET(no).SliceThickness) * ...
  cross(SET(no).ImageOrientation(4:6),SET(no).ImageOrientation(1:3))'];

%------------------------------------------------------------------
function [x,y,z] = cyl2cart(xin,yin,slice,numslices,rotationcenter)
%------------------------------------------------------------------
%Convert from cylindrical to cartesian coordinates

z = xin;
rad = (yin-rotationcenter);
angle = (slice-1)*pi/numslices;
x = rad.*sin(angle)+rotationcenter;
y = rad.*cos(angle)+rotationcenter;

%--------------------------------------------------------
function [window,level] = con2win(contrast,brightness,no)
%--------------------------------------------------------
%Convert from contrast to window & level

global SET NO
if nargin < 3
  no = NO;
  if nargin < 2
    if nargin == 1
      no = contrast;
    end
    contrast = SET(no).IntensityMapping.Contrast;
    brightness = SET(no).IntensityMapping.Brightness;
  end
end

if isa(SET(no).IM,'single') || isa(SET(no).IM,'double')
  window = SET(no).IntensityScaling/contrast;
  level = SET(no).IntensityOffset + window*(1-brightness);
else
  window = (SET(no).maxValue-SET(no).minValue)/contrast;
  level = SET(no).minValue + window*(1-brightness);
end

%--------------------------------------------------------
function [contrast,brightness] = win2con(window,level,no)
%--------------------------------------------------------
%Inverse of con2win

global SET NO
if nargin < 3
  no = NO;
  if nargin == 1
    no = window;
  end
end

if isa(SET(no).IM,'single') || isa(SET(no).IM,'double')
  contrast = SET(no).IntensityScaling/window;
  brightness = 1-(level-SET(no).IntensityOffset)/window;
else
  contrast = (SET(no).maxValue-SET(no).minValue)/window;
  brightness = 1-(level-SET(no).minValue)/window;
end

%-----------------------------------
function [x,y] = sax2gla(xi,yi,zi,no)
%-----------------------------------
%Convert image coordinates from short-axis to GLA view
global SET
ang = SET(no).GLA.angle;
res = SET(no).ResolutionY*cos(ang)+SET(no).ResolutionX*abs(sin(ang));
x = zi;
y = 1/res*abs((yi-SET(no).GLA.y0)*cos(ang)*SET(no).ResolutionY + ...
  (xi-SET(no).GLA.x0)*sin(ang)*SET(no).ResolutionX);

%-----------------------------------
function [x,y,z] = gla2sax(xi,yi,no)
%-----------------------------------
%Convert image coordinates from GLA to short-axis view
global SET
ang = SET(no).GLA.angle;
res = SET(no).ResolutionY*cos(ang)+SET(no).ResolutionX*abs(sin(ang));
z = xi;
x = SET(no).GLA.x0 + sin(SET(no).GLA.angle)*yi*res/SET(no).ResolutionX;
y = SET(no).GLA.y0 + cos(SET(no).GLA.angle)*yi*res/SET(no).ResolutionY;


%-------------------------------------------
function [mantelarea, sliceMantelArea] = calcmantelarea(no)
%-------------------------------------------
%Calculate the mantel area in cm^2 of the LV endocardium in each timeframe. 
global SET

if sum(findfunctions('findslicewithendo',no))==0
  myfailed('No endocardial segmentation found.');
  mantelarea = 0;
  return;
end

mantelarea=zeros(SET(no).TSize,SET(no).ZSize);

for timeframe=1:SET(no).TSize %loop through all timeframes
  slices = find(findfunctions('findslicewithendo',no,timeframe)); %find slices with endo segmentation
  for slice=slices %loop through all slices with endo
    x=[SET(no).EndoX(:,timeframe, slice); SET(no).EndoX(1,timeframe, slice)];
    y=[SET(no).EndoY(:,timeframe, slice); SET(no).EndoY(1,timeframe, slice)];
    circumference=sum(sqrt((diff(x)*(SET(no).ResolutionX)).^2+(diff(y)*SET(no).ResolutionY).^2));
    mantelarea(timeframe,slice)=circumference*(SET(no).SliceThickness+SET(no).SliceGap)*0.01; %time 0.01 to transform from mm^2 cm^2
  end
end

sliceMantelArea=mantelarea';
mantelarea=sum(mantelarea,2)';


%--------------------
function calcflow(no) %#ok<DEFNU>
%--------------------
%calculates flow from ROIs

global DATA SET

nom = SET(no).Flow.MagnitudeNo;
nop = SET(nom).Flow.PhaseNo;

warnedempty = false;
outsize = [SET(nom).XSize SET(nom).YSize];
roinbr = SET(nom).RoiCurrent;

if length(SET(nom).Roi(roinbr).T)==1  %For non-time resolved Phase Contrast images
  SET(nom).TIncr=1;
end

for tloop = SET(nom).Roi(roinbr).T
  %Create mask
  mask = logical(segment('createmask',...
    outsize,...
    SET(nom).Roi(roinbr).Y(:,tloop),...
    SET(nom).Roi(roinbr).X(:,tloop)));
  
  %Extract phase image
  temp = SET(nop).IM(:,:,tloop,SET(nom).Roi(roinbr).Z);
  
  %If empty phasecorr, the do not add phase correction.
  if isempty(SET(nop).Flow.PhaseCorr)
    veldata = SET(nom).Roi(roinbr).Sign*...
      (temp-0.5)*2*SET(nop).VENC;
  else
    %Phase correction
    if SET(nop).Flow.PhaseCorrTimeResolved
      %Time resolved phase correction
      veldata = SET(nom).Roi(roinbr).Sign*...
        (temp-0.5-SET(nop).Flow.PhaseCorr(:,:,tloop,SET(nom).Roi(roinbr).Z))*2*SET(nop).VENC;
    else
      %Stationary phase correction
      veldata = SET(nom).Roi(roinbr).Sign*...
        (temp-0.5-SET(nop).Flow.PhaseCorr(:,:,1,SET(nom).Roi(roinbr).Z))*2*SET(nop).VENC;
    end;
  end;
  
  veldata = veldata(mask);
  if isempty(veldata)
    if not(warnedempty)
      mywarning('Empty ROI detected. Should not occur.',DATA.GUI.Segment);
    end;
    warnedempty = true;
  else
    posveldata = veldata(veldata>0);
    negveldata = veldata(veldata<0);
    SET(nom).Flow.Result(roinbr).velmean(tloop) = mean(veldata);
    SET(nom).Flow.Result(roinbr).velstd(tloop) = std(veldata);
    SET(nom).Flow.Result(roinbr).velmax(tloop) = max(veldata);
    SET(nom).Flow.Result(roinbr).velmin(tloop) = min(veldata);
    SET(nom).Flow.Result(roinbr).kenergy(tloop) = sum((veldata/100).^3/2*(SET(nom).ResolutionX*SET(nom).ResolutionY/1e6)*1060); %kg/m^3
    SET(nom).Flow.Result(roinbr).area(tloop) = SET(nom).ResolutionX*SET(nom).ResolutionY*sum(mask(:))/100; %cm^2
    SET(nom).Flow.Result(roinbr).netflow(tloop) = (10/1000)*SET(nom).ResolutionX*SET(nom).ResolutionY*sum(veldata); %cm^3
    SET(nom).Flow.Result(roinbr).posflow(tloop) = (10/1000)*SET(nom).ResolutionX*SET(nom).ResolutionY*sum(posveldata); %cm^3
    SET(nom).Flow.Result(roinbr).negflow(tloop) = (10/1000)*SET(nom).ResolutionX*SET(nom).ResolutionY*sum(negveldata); %cm^3
  end;
end;
SET(nom).Flow.Result(roinbr).diameter = NaN;
hr = (60/(SET(nom).TSize*SET(nom).TIncr));

netflow = SET(nom).Flow.Result(roinbr).netflow;
timeframes = SET(nom).Roi(roinbr).T;%1:SET(nom).TSize;
TIncr = SET(nom).TIncr;
%Sum
SET(nom).Flow.Result(roinbr).nettotvol = nansum(netflow(timeframes))*TIncr;
SET(nom).Flow.Result(roinbr).netforwardvol = nansum(netflow(timeframes).*(netflow(timeframes)>0))*TIncr;
SET(nom).Flow.Result(roinbr).netbackwardvol = nansum(netflow(timeframes).*(netflow(timeframes)<0))*TIncr;
SET(nom).Flow.Result(roinbr).regfrac = abs(100*SET(nom).Flow.Result(roinbr).netbackwardvol/SET(nom).Flow.Result(roinbr).netforwardvol);
SET(nom).Flow.Result(roinbr).sv = SET(nom).Flow.Result(roinbr).nettotvol/hr*60;


%---------------------------
function clearflow(no,roinbr) %#ok<DEFNU>
%---------------------------
%clear ROI flow from ROI with roinbr

global SET NO

if nargin == 0
  no = NO;
elseif nargin == 1
  roinbr = SET(no).RoiCurrent;
end
if ~isempty(SET(no).Flow)
  if ~isempty(SET(no).Flow.Result)
    if length(SET(no).Flow.Result) == 1
      clearallflow(no);
    else
      SET(no).Flow.Result(roinbr) = [];
    end
  end
end



%------------------------
function clearallflow(no) 
%------------------------
%clear all ROI flow from no
global SET

if nargin == 0
  no = NO;
end

if ~isempty(SET(no).Flow)
  SET(no).Flow.Result = [];
end
