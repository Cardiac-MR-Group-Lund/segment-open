function [varargout] = specttools(varargin)
%------------------------------------------
%tools for analysis of SPECT images
%
%written by Helen Soneson 2010-01-26

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%------------------------------------------------------
function [U,P] = globalsurfapprox1D_callback(Q,p,n) %#ok<DEFNU>
%------------------------------------------------------
%global surface approximation algorithm
%input
%Q:   point cloud
%p:   degree of the nonrational surface
%n:   (n+1) is the number of control points
%output
%U:   knots
%P:   final control points
%
%Global surface approximation with fixed numer of control points.
%The algorithm interpolates the 4 corner points of Q exactly and
%approximates the remaining Q. It uses repeated least squares curve fits by
%fitting the (r+1) rows of data, which resulting in the
%final control points for the surface.
%The implementation was done according to "The NURBS Book" by Piegl Tiller
%
%Written by Helen Soneson 2014-09-26

%(r+1)x(s+1) is the number of points in the point cloud Q
r = size(Q,1)-1;
uk = SurfMeshParams1D(r,Q);  %according to algorithm 9.3

%compute knots U
%define the internal knots according to eq 9.68 and 9.69
U = zeros(n+1+p+1,1);
U(end-p-1:end,1) = ones(p+2,1);
d = (r+1)/(n-p+1);
for j = 1:n-p
  i = floor(j*d);
  alpha = j*d-i;
  U(p+j) = (1-alpha)*uk(i)+alpha*uk(i+1);
end

%Compute Nu, the (r-1)x(n-1) matrix of scalars from the basis functions
Nu = zeros(n-1,r-1);
for i = 1:n-1
  for uloop = 1:r-1
    u = uk(uloop+1);
    tempNu = fastbasisfct(u,i,p,U);
    Nu(i,uloop) = tempNu;
  end
end
NuTNu = Nu*Nu';

%u-directional fits
P = zeros(n+1,3);  %control points in the u-direction
P(1,:) = Q(1,:);  %interpolates the boundary points
P(n+1,:) = Q(r+1,:);  %interpolates the boundary points
R = zeros(r-1,3);  %R-values according to eq 9.63
for k = 2:r
  R(k-1,:) = Q(k,:)-fastbasisfct(uk(k),0,p,U).*Q(1,:)-fastbasisfct(uk(k),n,p,U).*Q(r+1,:);
end
Ru(:,1) = (R(:,1)'*Nu')';  %R-vector according to 9.67
Ru(:,2) = (R(:,2)'*Nu')';
Ru(:,3) = (R(:,3)'*Nu')';
%solve the system with invers
P(2:n,1) = NuTNu\Ru(:,1);
P(2:n,2) = NuTNu\Ru(:,2);
P(2:n,3) = NuTNu\Ru(:,3);



%---------------------------------------
function uk = SurfMeshParams1D(r,Q)
%---------------------------------------
%calculate paramters for global surfacde approximation according to
%algorithm 9.3
%r:     (r+1) is the number of control points
%Q:     point cloud
%uk:    the calculated parameters

%uk
uk(1:r) = zeros(r,1);
uk(r+1) = 1;

total = 0;  %total chord length of row
cds = zeros(1,r);
for k = 2:r+1
  cds(k) = sqrt((Q(k,1)-Q(k-1,1))^2+(Q(k,2)-Q(k-1,2))^2+(Q(k,3)-Q(k-1,3))^2);
  total = total+cds(k);
end
if total == 0
  disp('zero')
else
  d = 0;
  for t = 2:r
    d = d+cds(t);
    uk(t) = uk(t)+d/total;
  end
end



%------------------------------------------------------
function [U,V,P] = globalsurfapprox_callback(Q,p,q,n,m) %#ok<DEFNU>
%------------------------------------------------------
%global surface approximation algorithm
%input
%Q:   point cloud
%p,q: degree of the nonrational surface
%n,m: (n+1)x(m+1) is the number of control points
%output
%U,V: knots
%P:   final control points
%
%Global surface approximation with fixed numer of control points.
%The algorithm interpolates the 4 corner points of Q exactly and
%approximates the remaining Q. It uses repeated least squares curve fits by
%starting to fit the (s+1) rows of data, which resulting in temporary
%control points. Then it fits across these control points to produce the
%final control points for the surface.
%The implementation was done according to "The NURBS Book" by Piegl Tiller
%
%Written by Helen Soneson 2009-12-01

%(r+1)x(s+1) is the number of points in the point cloud Q
r = size(Q,1)-1;
s = size(Q,2)-1;
[uk,vl] = SurfMeshParams(r,s,Q);  %according to algorithm 9.3

%compute knots U
%define the internal knots according to eq 9.68 and 9.69
U = zeros(n+1+p+1,1);
U(end-p-1:end,1) = ones(p+2,1);
d = (r+1)/(n-p+1);
for j = 1:n-p
  i = floor(j*d);
  alpha = j*d-i;
  U(p+j) = (1-alpha)*uk(i)+alpha*uk(i+1);
end

%compute knots V
%define the internal knots according to eq 9.68 and 9.69
V = zeros(m+1+q+1,1);
V(end-q-1:end,1) = ones(q+2,1);
d = (s+1)/(m-q+1);
for j = 1:m-q
  i = floor(j*d);
  alpha = j*d-i;
  V(q+j) = (1-alpha)*vl(i)+alpha*vl(i+1);
end

%Compute Nu, the (r-1)x(n-1) matrix of scalars from the basis functions
Nu = zeros(n-1,r-1);
for i = 1:n-1
  for uloop = 1:r-1
    u = uk(uloop+1);
    tempNu = basisfct(u,i,p,U);
    Nu(i,uloop) = tempNu;
  end
end
NuTNu = Nu*Nu';

%u-directional fits
tempP = zeros(n+1,s+1,3);  %temporal control points in the u-direction
for j = 1:s+1
  tempP(1,j,:) = Q(1,j,:);  %interpolates the boundary points
  tempP(n+1,j,:) = Q(r+1,j,:);  %interpolates the boundary points
  R = zeros(r-1,3);  %R-values according to eq 9.63
  for k = 2:r
    R(k-1,:) = Q(k,j,:)-basisfct(uk(k),0,p,U).*Q(1,j,:)-basisfct(uk(k),n,p,U).*Q(r+1,j,:);
  end
  Ru(:,1) = (R(:,1)'*Nu')';  %R-vector according to 9.67
  Ru(:,2) = (R(:,2)'*Nu')';
  Ru(:,3) = (R(:,3)'*Nu')';
  %solve the system with invers
  tempP(2:n,j,1) = NuTNu\Ru(:,1);
  tempP(2:n,j,2) = NuTNu\Ru(:,2);
  tempP(2:n,j,3) = NuTNu\Ru(:,3);
end

%Compute Nv, the (s-1)x(m-1) matrix of scalars from the basis functions
Nv = zeros(m-1,s-1);
for i = 1:m-1
  for vloop = 1:s-1
    v = vl(vloop+1);
    tempNv = basisfct(v,i,q,V);
    Nv(i,vloop) = tempNv;
  end
end
NvTNv = Nv*Nv';

%v-directional fits
P = zeros(n+1,m+1,3);  %the final control points
for i = 1:n+1
  P(i,1,:) = tempP(i,1,:);
  P(i,m+1,:) = tempP(i,s+1,:);
  R = zeros(s-1,3);
  for l = 2:s
    R(l-1,:) = tempP(i,l,:)-basisfct(vl(l),0,q,V).*tempP(i,1,:)-basisfct(vl(l),m,q,V).*tempP(i,s+1,:);
  end
  Rv(:,1) = (R(:,1)'*Nv')';
  Rv(:,2) = (R(:,2)'*Nv')';
  Rv(:,3) = (R(:,3)'*Nv')';
  P(i,2:m,1) = NvTNv\Rv(:,1);
  P(i,2:m,2) = NvTNv\Rv(:,2);
  P(i,2:m,3) = NvTNv\Rv(:,3);
end

%---------------------------------------
function [uk,vl] = SurfMeshParams(r,s,Q)
%---------------------------------------
%calculate paramters for global surfacde approximation according to
%algorithm 9.3
%r,s:   (r+1)x(s+1) is the number of control points
%Q:     point cloud
%uk,vl: the calculated parameters

%uk
num = s+1;
uk(1:r) = zeros(r,1);
uk(r+1) = 1;
for l = 1:s+1
  total = 0;  %total chord length of row
  cds = zeros(1,r);
  for k = 2:r+1
    cds(k) = sqrt((Q(k,l,1)-Q(k-1,l,1))^2+(Q(k,l,2)-Q(k-1,l,2))^2+(Q(k,l,3)-Q(k-1,l,3))^2);
    total = total+cds(k);
  end
  if total == 0
    num = num-1;
  else
    d = 0;
    for t = 2:r
      d = d+cds(t);
      uk(t) = uk(t)+d/total;
    end
  end
end
if num == 0
  disp('error');
  return;
end
for q = 2:r
  uk(q) = uk(q)/num;
end
    
%vl
%equal distribuated along the long-axis
vl = 0:1/s:1;


%--------------------------------------------------------------
function [U,V,O,P] = globalsurfapprox3D_callback(Q,p,q,g,n,m,h) %#ok<DEFNU>
%--------------------------------------------------------------
%global surface approximation algorithm
%input
%Q:     point cloud
%p,q,g: degree of the nonrational surface
%n,m,h: (n+1)x(m+1)x(h+1) is the number of control points
%output
%U,V,O: knots
%P:     final control points
%
%Global surface approximation with fixed numer of control points.
%The algorithm interpolates the 4 corner points of Q exactly and
%approximates the remaining Q. It uses repeated least squares curve fits by
%starting to fit the (s+1) rows of data, which resulting in temporary
%control points. Then it fits across these control points to produce the
%final control points for the surface.
%The implementation was done according to "The NURBS Book" by Piegl Tiller
%
%Written by Helen Soneson 2009-12-01

%(r+1)x(s+1)x(t+1) is the number of points in the point cloud Q
r = size(Q,1)-1;
s = size(Q,2)-1;
t = size(Q,3)-1;
[uk,vl,xl] = SurfMeshParams3D(r,s,t,Q);  %uk according to chord length, vl and xl equally distribuated

%compute knots U
%define the internal knots according to eq 9.68 and 9.69
U = zeros(n+1+p+1,1);
U(end-p-1:end,1) = ones(p+2,1);
d = (r+1)/(n-p+1);
for j = 1:n-p
  i = floor(j*d);
  alpha = j*d-i;
  U(p+j) = (1-alpha)*uk(i)+alpha*uk(i+1);
end

%compute knots V
V = zeros(m+1+q+1,1);
V(end-q-1:end,1) = ones(q+2,1);
d = (s+1)/(m-q+1);
for j = 1:m-q
  i = floor(j*d);
  alpha = j*d-i;
  V(q+j) = (1-alpha)*vl(i)+alpha*vl(i+1);
end

%compute knots O
O = zeros(h+1+g+1,1);
O(end-g-1:end,1) = ones(g+2,1);
d = (t+1)/(h-g+1);
for j = 1:h-g
  i = floor(j*d);
  alpha = j*d-i;
  O(g+j) = (1-alpha)*xl(i)+alpha*xl(i+1);
end

%Compute Nx, the (t-1)x(h-1) matrix of scalars from the basis functions
Nx = zeros(h-1,t-1);
for i = 1:h-1
  for xloop = 1:t-1
    x = xl(xloop+1);
    tempNx = basisfct(x,i,g,O);
    Nx(i,xloop) = tempNx;
  end
end
NxTNx = Nx*Nx';

%Compute Nu, the (r-1)x(n-1) matrix of scalars from the basis functions
Nu = zeros(n-1,r-1);
for i = 1:n-1
  for uloop = 1:r-1
    u = uk(uloop+1);
    tempNu = basisfct(u,i,p,U);
    Nu(i,uloop) = tempNu;
  end
end
NuTNu = Nu*Nu';

%Compute Nv, the (s-1)x(m-1) matrix of scalars from the basis functions
Nv = zeros(m-1,s-1);
for i = 1:m-1
  for vloop = 1:s-1
    v = vl(vloop+1);
    tempNv = basisfct(v,i,q,V);
    Nv(i,vloop) = tempNv;
  end
end
NvTNv = Nv*Nv';

%x-directional fits
tempP = zeros(r+1,s+1,h+1,3);  %temporal control points in the u-direction
for j = 1:s+1
  for i = 1:r+1
    tempP(i,j,1,:) = Q(i,j,1,:);
    tempP(i,j,h+1,:) = Q(i,j,t+1,:);
    R = zeros(t-1,3);
    for f = 2:t
      R(f-1,:) = Q(i,j,f,:)-basisfct(xl(f),0,g,O).*Q(i,j,1,:)-basisfct(xl(f),h,g,O).*Q(i,j,t+1,:);
    end
    Rx(:,1) = (R(:,1)'*Nx')';
    Rx(:,2) = (R(:,2)'*Nx')';
    Rx(:,3) = (R(:,3)'*Nx')';
    tempP(i,j,2:h,1) = NxTNx\Rx(:,1);
    tempP(i,j,2:h,2) = NxTNx\Rx(:,2);
    tempP(i,j,2:h,3) = NxTNx\Rx(:,3);
  end
end

%u-directional fits
temp2P = zeros(n+1,s+1,h+1,3);  %temporal control points in the u-direction
for time = 1:h+1
  for j = 1:s+1
    temp2P(1,j,time,:) = tempP(1,j,time,:);  %interpolates the boundary points
    temp2P(n+1,j,time,:) = tempP(r+1,j,time,:);  %interpolates the boundary points
    R = zeros(r-1,3);  %R-values according to eq 9.63
    for k = 2:r
      R(k-1,:) = tempP(k,j,time,:)-basisfct(uk(k),0,p,U).*tempP(1,j,time,:)-basisfct(uk(k),n,p,U).*tempP(r+1,j,time,:);
    end
    Ru(:,1) = (R(:,1)'*Nu')';  %R-vector according to 9.67
    Ru(:,2) = (R(:,2)'*Nu')';
    Ru(:,3) = (R(:,3)'*Nu')';
    %solve the system with invers
    temp2P(2:n,j,time,1) = NuTNu\Ru(:,1);
    temp2P(2:n,j,time,2) = NuTNu\Ru(:,2);
    temp2P(2:n,j,time,3) = NuTNu\Ru(:,3);
  end
end

%v-directional fits
P = zeros(n+1,m+1,h+1,3);  %the final control points
for time = 1:h+1
  for i = 1:n+1
    P(i,1,time,:) = temp2P(i,1,time,:);
    P(i,m+1,time,:) = temp2P(i,s+1,time,:);
    R = zeros(s-1,3);
    for l = 2:s
      R(l-1,:) = temp2P(i,l,time,:)-basisfct(vl(l),0,q,V).*temp2P(i,1,time,:)-basisfct(vl(l),m,q,V).*temp2P(i,s+1,time,:);
    end
    Rv(:,1) = (R(:,1)'*Nv')';
    Rv(:,2) = (R(:,2)'*Nv')';
    Rv(:,3) = (R(:,3)'*Nv')';
    P(i,2:m,time,1) = NvTNv\Rv(:,1);
    P(i,2:m,time,2) = NvTNv\Rv(:,2);
    P(i,2:m,time,3) = NvTNv\Rv(:,3);
  end
end



%----------------------------------------------
function [uk,vl,xl] = SurfMeshParams3D(r,s,t,Q)
%----------------------------------------------
%calculate paramters for global surfacde approximation according to
%algorithm 9.3
%r,s,t:     (r+1)x(s+1)x(t+1) is the number of control points
%Q:         point cloud
%uk,vl,xl:  the calculated parameters

%uk
num = s+1;
uk(1:r) = zeros(r,1);
uk(r+1) = 1;
for l = 1:s+1
  total = 0;  %total chord length of row
  cds = zeros(1,r);
  for k = 2:r+1
    cds(k) = sqrt((Q(k,l,1)-Q(k-1,l,1))^2+(Q(k,l,2)-Q(k-1,l,2))^2+(Q(k,l,3)-Q(k-1,l,3))^2);
    total = total+cds(k);
  end
  if total == 0
    num = num-1;
  else
    d = 0;
    for i = 2:r
      d = d+cds(i);
      uk(i) = uk(i)+d/total;
    end
  end
end
if num == 0
  disp('error');
  return;
end
for q = 2:r
  uk(q) = uk(q)/num;
end
    
%vl
%equal distribuated along the long-axis
vl = 0:1/s:1;

%xl
%equal distribuated over time
xl = 0:1/t:1;


%-------------------------------------
function N = basisfct(u,i,degree,grid)
%-------------------------------------
%calculate the B-spline basis function
%u:       current knot
%i:       index
%degree:  degree of the B-spline 
%grid:    the grid

% Add two additional knots at the ends of the knot sequence
extgrid = [grid(1); grid; grid(size(grid,2))];

% Degree = 0
if (degree == 0)
   if u >= extgrid(i+1) && u < extgrid(i+2)
       N = 1;
   else
       N = 0;
   end
% Degree != 0
else
   cterm1 = extgrid(i+1+degree)-extgrid(i+1);
   cterm2 = u-extgrid(i+1);
   if (cterm1 == 0) || (cterm2 == 0)
       term1 = 0;
   else
       term1 = cterm2/cterm1*basisfct(u,i,degree-1,grid);
   end
   cterm3 = extgrid(i+2+degree)-extgrid(i+2);
   cterm4 = extgrid(i+2+degree)-u;
   if (cterm3 == 0) || (cterm4 == 0)
       term2 = 0;
   else
       term2 = cterm4/cterm3*basisfct(u,i+1,degree-1,grid);
   end
   N = term1 + term2;  
end

%---------------------------------------------------------------------
function S = surfacepoints3D_Callback(p,q,g,P,U,V,O,nbrpU,nbrpV,nbrpO) %#ok<DEFNU>
%---------------------------------------------------------------------
%construct a surface from control points
%input
%p,q:         degree of the nonrational surface
%P:           control points
%U,V:         knot vectors
%nbrpU,nbrpV: number of points in the constructed surface
%output
%S:           points on the resulting surface

stepsizeU = (U(end)-U(1))/(nbrpU-1);
stepsizeV = (V(end)-V(1))/(nbrpV-1);
stepsizeO = (O(end)-O(1))/(nbrpO-1);
tempcp = zeros(size(P,3),3);
newcp = zeros(size(P,1),3);
S = zeros(nbrpU,nbrpV,nbrpO,3);
tempP = zeros(size(P,2),size(P,4));
uloop = 1;
for u = U(1):stepsizeU:(U(end))
  vloop = 1;
  for v =  V(1):stepsizeV:(V(end))
    xloop = 1;
    for x = O(1):stepsizeO:(O(end))
      for i = 0:size(P,1)-1
        for j = 0:size(P,3)-1
          tempP(:,:) = P(i+1,:,j+1,:);
          tempcp(j+1,:) = deBoor(tempP,V,v,q);
        end
        newcp(i+1,:) = deBoor(tempcp,O,x,g);
      end
      S(uloop,vloop,xloop,:) = deBoor(newcp,U,u,p);
      xloop = xloop+1;
    end    
    vloop = vloop+1;
  end
  uloop = uloop+1;
end


%---------------------------------------------------------
function S = surfacepoints_Callback(p,q,P,U,V,nbrpU,nbrpV) %#ok<DEFNU>
%---------------------------------------------------------
% construct a surface from control points
% input
% p,q:         degree of the nonrational surface
% P:           control points
% U,V:         knot vectors
% nbrpU,nbrpV: number of points in the constructed surface
% output
% S:           points on the resulting surface

stepsizeU = (U(end)-U(1))/(nbrpU-1);
stepsizeV = (V(end)-V(1))/(nbrpV-1);
loop = 1;
b = zeros(nbrpU*nbrpV,3);
newcp = zeros(size(P,1),3);
uloop = 1;
S = zeros(nbrpU,nbrpV,3);
tempP = zeros(size(P,2),size(P,3));
for u = U(1):stepsizeU:(U(end))
  vloop = 1;
  for v =  V(1):stepsizeV:(V(end))
    for i = 0:size(P,1)-1
      tempP(:,:) = P(i+1,:,:);
      newcp(i+1,:) = deBoor(tempP,V,v,q);
    end
    b(loop,:) = deBoor(newcp,U,u,p);
    S(uloop,vloop,:) = b(loop,:);
    loop = loop+1;
    vloop = vloop+1;
  end
  uloop = uloop+1;
end

%--------------------------------------
function bsplinecurve = deBoor(P,U,u,n)
%--------------------------------------
%deBoor algorithm for calculating a point on the surface
%P:   control points
%U:   knot vector
%u:   current knot
%n:   degree of the nonrational surface

%determine the position of u in the knot sequence
if u == U(end)
  posu = find(u == U,1,'first')-1;
else
  posu = find(u < U,1,'first')-1;
end
%find the corresponding control points
d = P(posu-n+1:posu+1,:);

for r = 1:n
  for i = 1:n-r+1
    te = U(posu+(i-1)+1);
    ts = U(posu+r+(i-1)-n);
    alpha = (u-ts)/(te-ts);
    d(i,:) = (1-alpha).*d(i,:)+(alpha).*d(i+1,:);
  end
end
bsplinecurve = d(1,:);



%------------------------------------------------------------------------
function compIM = compensatecounts(no,mask,t,startoutflow,endoutflow,compIM, ...
  startslice,endslice,percentil95maxcount,percentil100maxcount,cropim,lumenmask) %#ok<DEFNU>
%------------------------------------------------------------------------
%compensate the counts in the image due to thin basal and apical walls
%used in the segmentation of MaR and ischemia (rest-stress change) in SPECT images
%INPUT
%no:            thumbnail number
%mask:          myocardial mask
%startoutflow:  start pixel for the outflow tract in the basal slices
%endoutflow:    end pixel for the outflow tract in the basal slices
%compIM:        original intensity in the image stack
%startslice:    first slice with lv segmentation
%endslice:      last slice with lv segmentation
%percentil95maxcount:   95 percentile of maximal count in the image
%percentil100maxcount:  100 percentile of maximal count in the image
%cropim:        
%OUTPUT
%compIM:        the compensated intensity in the image stack
%
%Written by Helen Soneson 2011-08-24

global SET

%%%%% calculate the compensation values %%%%%
%compensation values for basal slices and around the outflow tract
percentile95 = zeros(1,length(startoutflow)+2);
nbrofbasal = min(endslice-startslice-2,length(startoutflow));
for basalloop = startslice:startslice+nbrofbasal+1
  if ndims(compIM) == 3
    immask = mask(:,:,basalloop).*compIM(:,:,basalloop);
  else
    immask = mask(:,:,t,basalloop).*compIM(:,:,t,basalloop);
  end
  sortmyocounts = sort(immask(:));
  sortmyocounts = sortmyocounts(sortmyocounts>0);
  percentile95(basalloop-startslice+1) = sortmyocounts(round(0.95*end));
end
if length(percentile95) == 1
  percentile95(2) = mean([percentile95 percentil95maxcount]);
end
percentile95percent = min(1,percentile95/percentil95maxcount);
fittedline = polyfit(1:length(percentile95percent),percentile95percent,1);
compensate1 = fittedline(1)+fittedline(2);
compensate2 = min([1 fittedline(1)*length(percentile95percent)+fittedline(2)]);
basalcompensate = [compensate1 compensate2];
if (nbrofbasal+2) ~= length(basalcompensate)
  %resample the basalcompensate
  x = linspace(0,1,length(basalcompensate));
  xi = linspace(0,1,nbrofbasal+2);
  basalcompensate = interp1(x,basalcompensate,xi);
end
%compensation values for the inferior part
halfnbrofslices = round((endslice-startslice+1)*0.75);
if strfind(lower(SET(no).PatientInfo.Sex),'f')
  inferiorcompensate = 1./interp1(linspace(0,1,3),[0.91 0.91 1],linspace(0,1,halfnbrofslices));
else
  inferiorcompensate = 1./interp1(linspace(0,1,3),[0.86 0.86 1],linspace(0,1,halfnbrofslices));
end
inferiorcompangle = [4.5*pi/4 7.5*pi/4];

%compensate whole basal slices, the region surrounding the outflow tract
%and in the inferior part
nbrofcompslices = ceil(15/(SET(no).SliceThickness+SET(no).SliceGap));
compensationvaluesbasal = 1./percentile95percent;
startangle = zeros(1,startslice+nbrofbasal+nbrofcompslices-1);
endangle = startangle;
for sliceloop = startslice:startslice+nbrofbasal-1  %the basal slices
  [startangle(sliceloop) endangle(sliceloop)] = calcangle(startoutflow(sliceloop-startslice+1), ...
    endoutflow(sliceloop-startslice+1),sliceloop+cropim(end,1)-1,no,t);
  compmaskinf = calccompmask(startangle(sliceloop)-pi/4,startangle(sliceloop), ...
    [1 mean([compensationvaluesbasal(sliceloop-startslice+1) 1]) compensationvaluesbasal(sliceloop-startslice+1)],sliceloop+cropim(3,1)-1,no,t);
  compmaskant = calccompmask(endangle(sliceloop),endangle(sliceloop)+pi/4, ...
    [compensationvaluesbasal(sliceloop-startslice+1) mean([compensationvaluesbasal(sliceloop-startslice+1) 1]) 1],sliceloop+cropim(3,1)-1,no,t);
  compmaskinferior = calccompmask(inferiorcompangle(1),inferiorcompangle(2), ...
    [1 inferiorcompensate(sliceloop-startslice+1) inferiorcompensate(sliceloop-startslice+1) 1],sliceloop+cropim(3,1)-1,no,t);
  maxcompvalue = max(compensationvaluesbasal(sliceloop-startslice+1));
  compmask = compmaskinf+compmaskant+compmaskinferior-2;
  compmask(compmask > maxcompvalue) = maxcompvalue;
  compmask = compmask(cropim(1,1):cropim(1,2),cropim(2,1):cropim(2,2));
  if ndims(compIM) == 3
    temporalIM = compIM(:,:,sliceloop).*compmask/basalcompensate(sliceloop-startslice+1);
    temporalIM(temporalIM > percentil100maxcount) = percentil100maxcount;
    compIM(:,:,sliceloop) = temporalIM;
  else
    temporalIM = compIM(:,:,t,sliceloop).*compmask/basalcompensate(sliceloop-startslice+1);
    temporalIM(temporalIM > percentil100maxcount) = percentil100maxcount;
    compIM(:,:,t,sliceloop) = temporalIM;
  end
end

%compensate in the region surrounding the outflow tract in the slices following the basal slices
%calculate the compensation value based on the lumen counts
if nargin >= 12
  nbrofslices = endslice-startslice+1;
  lumenslices = round(startslice+nbrofslices*0.25):round(startslice+nbrofslices*0.75);
  lumenmasknan = lumenmask;
  lumenmasknan(lumenmasknan==0) = NaN;
  lumencounts = lumenmasknan(:,:,lumenslices).*compIM(:,:,lumenslices);
  lowestcompvalue = mynanmean(lumencounts(:));
  tempcompensate = interp1([0 1],[lowestcompvalue 1],linspace(0,1,nbrofcompslices));
else
  %OLD compensation value
  tempcompensate = interp1(linspace(0,1,length(percentile95percent)),percentile95percent,linspace(0,1,nbrofcompslices+1));
  tempcompensate = tempcompensate(1:end-1);
end
compensationvalues = 1./tempcompensate;
for sliceloop = startslice+nbrofbasal:min([endslice,startslice+nbrofbasal+nbrofcompslices-1,startslice+length(inferiorcompensate)-1])
  slice = min([length(startoutflow) sliceloop-startslice]);
  [startangle(sliceloop) endangle(sliceloop)] = calcangle(startoutflow(slice),endoutflow(slice),sliceloop+cropim(3,1)-1,no,t);
  %the slices following the basal slices
  compmaskoutflow = calccompmask(startangle(sliceloop)-pi/3,endangle(sliceloop)+pi/3, ...
    [1 compensationvalues(sliceloop-startslice-nbrofbasal+1) 1],sliceloop+cropim(end,1)-1,no,t); 
  compmaskinferior = calccompmask(inferiorcompangle(1),inferiorcompangle(2), ...
    [1 inferiorcompensate(sliceloop-startslice+1) inferiorcompensate(sliceloop-startslice+1) 1],sliceloop+cropim(3,1)-1,no,t);
  compmask = compmaskoutflow+compmaskinferior-1;
  compmask = compmask(cropim(1,1):cropim(1,2),cropim(2,1):cropim(2,2));
  if ndims(compIM) == 3
    if sliceloop <= startslice+nbrofbasal+1
      temporalIM = compIM(:,:,sliceloop).*compmask/basalcompensate(sliceloop-startslice+1);
      temporalIM(temporalIM > percentil100maxcount) = percentil100maxcount;
      compIM(:,:,sliceloop) = temporalIM;
    else
      temporalIM = compIM(:,:,sliceloop).*compmask;
      temporalIM(temporalIM > percentil100maxcount) = percentil100maxcount;
      compIM(:,:,sliceloop) = temporalIM;
    end
  else
    if sliceloop <= startslice+nbrofbasal+1
      temporalIM = compIM(:,:,t,sliceloop).*compmask/basalcompensate(sliceloop-startslice+1);
      temporalIM(temporalIM > percentil100maxcount) = percentil100maxcount;
      compIM(:,:,t,sliceloop) = temporalIM;
    else
      temporalIM = compIM(:,:,t,sliceloop).*compmask;
      temporalIM(temporalIM > percentil100maxcount) = percentil100maxcount;
      compIM(:,:,t,sliceloop) = temporalIM;
    end
  end
end

%compensate inferior
for sliceloop = min(endslice,startslice+nbrofbasal+nbrofcompslices-1)+1:startslice+halfnbrofslices-1
  compmask = calccompmask(inferiorcompangle(1),inferiorcompangle(2), ...
    [1 inferiorcompensate(sliceloop-startslice+1) inferiorcompensate(sliceloop-startslice+1) 1],sliceloop+cropim(3,1)-1,no,t);
%   maxcompvalue = max(compensationvaluesbasal(sliceloop-startslice+1));
%   compmask(compmask > maxcompvalue) = maxcompvalue;
  compmask = compmask(cropim(1,1):cropim(1,2),cropim(2,1):cropim(2,2));
  if ndims(compIM) == 3
    temporalIM = compIM(:,:,sliceloop).*compmask;
    temporalIM(temporalIM > percentil100maxcount) = percentil100maxcount;
    compIM(:,:,sliceloop) = temporalIM;
  else
    temporalIM = compIM(:,:,t,sliceloop).*compmask;
    temporalIM(temporalIM > percentil100maxcount) = percentil100maxcount;
    compIM(:,:,t,sliceloop) = temporalIM;
  end
end

%compensation for low intensity in apex
%second next most apical slice
if ndims(compIM) == 3
  apicalcounts2 = compIM(:,:,endslice-2).*mask(:,:,endslice-2);
  apicalcounts1 = compIM(:,:,endslice-1).*mask(:,:,endslice-1);
  apicalcounts0 = compIM(:,:,endslice).*mask(:,:,endslice);
else
  apicalcounts2 = compIM(:,:,t,endslice-2).*mask(:,:,t,endslice-2);
  apicalcounts1 = compIM(:,:,t,endslice-1).*mask(:,:,t,endslice-1);
  apicalcounts0 = compIM(:,:,t,endslice).*mask(:,:,t,endslice);
end
apicalcounts2 = apicalcounts2(apicalcounts2 > 0);
meanapicalcounts2 = mean(apicalcounts2(:));
%next most apical slice
apicalcounts1 = apicalcounts1(apicalcounts1 > 0);
meanapicalcounts1 = mean(apicalcounts1(:));
%most apical slice
apicalcounts0 = apicalcounts0(apicalcounts0 > 0);
meanapicalcounts0 = mean(apicalcounts0(:));
%compensation
apicalcompensate = [mean([percentile95percent(1) percentile95percent(2)]) ...
  mean([percentile95percent(2) percentile95percent(min(end,3))])];
apicalcompensate1 = max([apicalcompensate(2) meanapicalcounts1/meanapicalcounts2]);
apicalcompensate0 = max([apicalcompensate(1) meanapicalcounts0/meanapicalcounts1]);
if ndims(compIM) == 3
  compIM(:,:,endslice-1) = compIM(:,:,endslice-1)/apicalcompensate1;
  compIM(:,:,endslice) = compIM(:,:,endslice)/apicalcompensate0;
else
  compIM(:,:,t,endslice-1) = compIM(:,:,t,endslice-1)/apicalcompensate1;
  compIM(:,:,t,endslice) = compIM(:,:,t,endslice)/apicalcompensate0;
end

%-------------------------------------------------------------------------
function compmask = calccompmask(startangle,endangle,compvalue,slice,no,t)
%-------------------------------------------------------------------------
%calculate the compensation mask for the current slice
%startangle:  start angle for the compensation region
%endangle:    end angle for the compensation region
%compvalue:   compensation values from startangle to endangle
%slice:       current slice
%no:          current image stack number
%t:           current time frame
%
%written by Helen Soneson 2009-03-31

global SET

%find the center point of LV in the current slice
centerx = round(mean(SET(no).EndoX(:,t,slice)))-1;
centery = round(mean(SET(no).EndoY(:,t,slice)))-1;
nx = floor(SET(no).XSize/2);  %number of points in the radial direction
ny = floor(SET(no).YSize/2);  %number of points in the radial direction
movex = nx-centerx;
movey = ny-centery;
compmask = zeros(SET(no).XSize,SET(no).YSize);

%calculate the angular matrix
[x,y] = ndgrid(linspace(-nx,nx,2*nx+1),linspace(-ny,ny,2*ny+1));
x = x+movex;
y = y+movey;
ang = mod(angle(complex(y,x))-pi-startangle,2*pi);  %angle
%correct the size
ang = ang(1:SET(no).XSize,1:SET(no).YSize);  
%define the compensation region
ang(ang > mod(endangle-startangle,2*pi)) = 0;
%define the compensation values
ang = round(ang*100);
maxangle = max(ang(:));
compensationvalues = interp1(linspace(1,maxangle,length(compvalue)),compvalue,1:maxangle);
for comploop = 1:length(compensationvalues)
  compindex = find(ang == comploop);
  compmask(compindex) = compensationvalues(comploop); %#ok<FNDSB>
end
compmask(compmask == 0) = 1;

%--------------------------------------------------------------------------
function [startangle endangle] = calcangle(startoutflow,endoutflow,slice,no,t)
%--------------------------------------------------------------------------
%calculate the angle for the compensation region based on the outflow tract
%startoutflow:  start point of the outflow tract
%endoutflow:    end point of the outflow tract
%slice:         number of the current basal slice
%
%written by Helen Soneson 2009-03-31

global SET

centerx = mean(SET(no).EndoX(:,:,slice));  %the center point of the LV
centery = mean(SET(no).EndoY(:,:,slice));  %the center point of the LV
%the start and end pixels for the outflow tract
startepix = SET(no).EpiX(startoutflow,t,slice);  
startepiy = SET(no).EpiY(startoutflow,t,slice);
endepix = SET(no).EpiX(endoutflow,t,slice);
endepiy = SET(no).EpiY(endoutflow,t,slice);

%find the start and end angle of the outflow tract
% startangle = atan((startepiy-centery)/(startepix-centerx+eps));
startangle = atan((startepix-centerx)/(startepiy-centery+eps));
if (startepiy-centery) > 0 && (startepix-centerx+eps) > 0
  startangle = pi+startangle;
elseif (startepiy-centery) < 0 && (startepix-centerx+eps) > 0
  startangle = 2*pi+startangle;
% elseif (startepiy-centery) < 0 && (startepix-centerx+eps) < 0
%   startangle = startangle;
elseif (startepiy-centery) > 0 && (startepix-centerx+eps) < 0
  startangle = pi+startangle;
end

% endangle = atan((endepiy-centery)/(endepix-centerx+eps));
endangle = atan((endepix-centerx)/(endepiy-centery+eps));
if (endepiy-centery) > 0 && (endepix-centerx+eps) > 0
  endangle = pi+endangle;
elseif (endepiy-centery) < 0 && (endepix-centerx+eps) > 0
  endangle = 2*pi+endangle;
% elseif (endepiy-centery) < 0 && (endepix-centerx+eps) < 0
%   endangle = endangle;
elseif (endepiy-centery) > 0 && (endepix-centerx+eps) < 0
  endangle = pi+endangle;
end



%------------------------------------
function newim = upsampleslices(f,im) %#ok<DEFNU>
%------------------------------------
%Upsamples along last dimension.
%Limitation input must be single or double

global DATA

%Find new size
switch ndims(im)
  case 4
    temp = squeeze(im(1,1,1,:));
  case 3
    temp = squeeze(im(1,1,:));    
  case 2
    temp = squeeze(im(1,:));    
  case 1
    temp = im;    
end

xi = linspace(1,length(temp),round(length(temp)*f));
newsize = length(interp1(temp,xi,'nearest'));

if newsize<1
  myfailed('Too few resulting slices. Aborting.',DATA.GUI.Segment);
  return;
end

%Reserve memory
zeroel = 0;
switch class(im)
  case 'double'
    zeroel = 0;
  case 'single'
    zeroel = single(0);
  case 'int16'
    zeroel = int16(0);
end;

switch ndims(im)
  case 4
    %x*y*t*z
      newim = repmat(zeroel,[...
        size(im,1) ...
        size(im,2) ...
        newsize ...
        size(im,3)]);
      %Loop over image volume
      for xloop=1:size(im,1)
        for yloop=1:size(im,2)
          newim(xloop,yloop,:,:) = interp1(...
            squeeze(im(xloop,yloop,:,:))',...
            xi,'linear');
        end
      end
      newim = permute(newim,[1 2 4 3]);
  case 3
    %a*b*z
      newim = repmat(zeroel,[...
        size(im,1) ...
        newsize ...
        size(im,2)]);
      %Loop over image volume     
      for xloop=1:size(im,1)
        newim(xloop,:,:) = interp1(...
          squeeze(im(xloop,:,:))',...
          xi,'linear');
      end    
      newim = permute(newim,[1 3 2]);
  case 2
    %a*z
    newim = interp1(im',xi,'linear')';
  case 1
    %z
    newim = interp1(im(:),xi,'linear');    
end


%----------------------------------------
function newim = upsamplevolume(f1,f2,im) %#ok<DEFNU>
%----------------------------------------
%Helper function to upsample an image spatially
%input
%f1:    upsample factor in x-direction
%f2:    upsample factor in y-direction
%im:    image to upsample
%output
%newim: upsampled image

%Find new size
newsize = size(imresize(im(:,:,1,1),[round(size(im,1)*f1) round(size(im,2)*f2)],'nearest'));

%Reserve memory
if isa(im,'single')
  newim = repmat(single(0),[newsize(1) newsize(2) size(im,3) size(im,4)]);
elseif isa(im,'int16')
  newim = repmat(int16(0),[newsize(1) newsize(2) size(im,3) size(im,4)]);
else
  newim = zeros(newsize(1),newsize(2),size(im,3),size(im,4));
end

%Loop over image volume
for tloop = 1:size(im,3)
  for zloop = 1:size(im,4)
    if isa(im,'single')
      newim(:,:,tloop,zloop) = single(imresize(im(:,:,tloop,zloop),[round(size(im,1)*f1) round(size(im,2)*f2)],'bicubic'));
    elseif isa(im,'int16')
      newim(:,:,tloop,zloop) = int16(imresize(im(:,:,tloop,zloop),[round(size(im,1)*f1) round(size(im,2)*f2)],'bicubic'));
    else
      newim(:,:,tloop,zloop) = imresize(im(:,:,tloop,zloop),[round(size(im,1)*f1) round(size(im,2)*f2)],'bicubic');
    end
  end
end


%--------------------------------------
function newline = upsampleline(f,line) %#ok<DEFNU>
%--------------------------------------
%Helper function to resample image stacks.
%factor 2 => x' = 2*x-0.5
%factor 3 => x' = 3*x-1
%factor 4 => x' = 4*x-1.5
%factor 5 => x' = 5*x-2
%factor 2.5 => x' = 2.5*x-1
%factor 3.5 => x' = 3.5*x-1.5
%factor 0.5 => x' = x'*0.5 (a odd)
%factor 0.5 => x' = x'*0.5+0.5 (a even)

if f == 1
    newline = line;
else
  if f>0
    d = (ceil(f)-1)/2;
  else
    d = 0;
  end

  if isa(line,'double')
    newline = line*f-d;
  else
    newline = zeros(size(line));
    for xloop = 1:size(line,1)
      for yloop = 1:size(line,2)
        newline{xloop,yloop} = line{xloop,yloop}*f-d;
      end
    end
  end
end


%--------------------------------------------------------------------------
function [endox,endoy,epix,epiy] = centerlinemethod(slice,no,t,nbrofpoints,nbrofminpoints, ...
  endoX,endoY,epiX,epiY) %#ok<DEFNU>
%--------------------------------------------------------------------------
%calculate the transmural corresponding endo, and epicardial points
%slice:           the current slice number
%no:              the current image stack number
%t:               the current time frame
%nbrofpoints:     the number of points along the mid-mural centerline
%nbrofminpoints:  the number of points to search
%
%written by Helen Soneson 2009-04-02

global SET

if nargin < 6
  %check so endo- and epicardium exist
  if isempty(SET(no).EndoX)
    myfailed('No endocardium exist');
    return;
  end
  if isempty(SET(no).EpiX)
    myfailed('No LV epicardium available.');
  end
  endoX = SET(no).EndoX;
  endoY = SET(no).EndoY;
  epiX = SET(no).EpiX;
  epiY = SET(no).EpiY;
end

%close the curves
tempendox = endoX(:,t,slice);
tempendoy = endoY(:,t,slice);
tempepix = epiX(:,t,slice);
tempepiy = epiY(:,t,slice);

if isnan(tempendox)
  tempendox = repmat(mean(tempepix),[size(tempendox,1) 1]);
  tempendoy = repmat(mean(tempepiy),[size(tempendoy,1) 1]);
end

centerlinex = mynanmean([tempendox';tempepix']);
centerliney = mynanmean([tempendoy';tempepiy']);
%resample to 100 points along the line
centerlinex = interp1(linspace(0,1,length(centerlinex)),centerlinex,linspace(0,1,nbrofpoints));
centerliney = interp1(linspace(0,1,length(centerliney)),centerliney,linspace(0,1,nbrofpoints));
%find the chords perpendicular to the centerline
tempk = (circshift(centerliney,[0 1])-circshift(centerliney,[0 -1]))./ ...
  (circshift(centerlinex,[0 1])-circshift(centerlinex,[0 -1]));
k = -1./tempk;
m = centerliney-k.*centerlinex;
%resample to 200 points along the endo- and epicardium
originalendox = interp1(linspace(0,1,size(tempendox,1)),tempendox,linspace(0,1,2*nbrofpoints));
originalendoy = interp1(linspace(0,1,size(tempendoy,1)),tempendoy,linspace(0,1,2*nbrofpoints));
originalepix = interp1(linspace(0,1,size(tempepix,1)),tempepix,linspace(0,1,2*nbrofpoints));
originalepiy = interp1(linspace(0,1,size(tempepiy,1)),tempepiy,linspace(0,1,2*nbrofpoints));
%find the endo- and epicardium point on the line perpendicular to the
%centerline
%endo
calcyendo = k'*originalendox+repmat(m,2*nbrofpoints,1)';
diffendo = abs(calcyendo-repmat(originalendoy,nbrofpoints,1));  %error from the perpendicular line
closeendo = sqrt((repmat(centerlinex,2*nbrofpoints,1)'-repmat(originalendox,nbrofpoints,1)).^2+ ...
  (repmat(centerliney,2*nbrofpoints,1)'-repmat(originalendoy,nbrofpoints,1)).^2);  %distance from the centerline
[~,minindexendo] = sort(closeendo,2);
endoindex = minindexendo(:,1:nbrofminpoints);  %find the 10 points closest to the centerline
%epi
calcyepi = k'*originalepix+repmat(m,2*nbrofpoints,1)';
diffepi = abs(calcyepi-repmat(originalepiy,nbrofpoints,1));  %error from the perpendicular line
closeepi = sqrt((repmat(centerlinex,2*nbrofpoints,1)'-repmat(originalepix,nbrofpoints,1)).^2+ ...
  (repmat(centerliney,2*nbrofpoints,1)'-repmat(originalepiy,nbrofpoints,1)).^2);  %distance from the centerline
[~,minindexepi] = sort(closeepi,2);
epiindex = minindexepi(:,1:nbrofminpoints);  %find the 10 points closest to the centerline
%find the endo- and epicardium point close to the each centerlinepoint
endoindexvalues = zeros(size(epiindex,1),nbrofminpoints);
epiindexvalues = zeros(size(epiindex,1),nbrofminpoints);
for minindexloop = 1:size(endoindex,1)
  endoindexvalues(minindexloop,:) = diffendo(minindexloop,endoindex(minindexloop,:));
  epiindexvalues(minindexloop,:) = diffepi(minindexloop,epiindex(minindexloop,:));
end
[~,tempindexendo] = min(endoindexvalues,[],2);
[~,tempindexepi] = min(epiindexvalues,[],2);
indexendo = zeros(1,size(endoindex,1));
indexepi = zeros(1,size(endoindex,1));
for indexloop = 1:size(endoindex,1)
  indexendo(indexloop) = endoindex(indexloop,tempindexendo(indexloop));
  indexepi(indexloop) = epiindex(indexloop,tempindexepi(indexloop));
end
%caculate the transmural corresponding endo- and epicardium
endox = originalendox(indexendo);
endoy = originalendoy(indexendo);
epix = originalepix(indexepi);
epiy = originalepiy(indexepi);


%--------------------------------------------------------------------------
function manuallyconfirm(no,originalno,resampledIM,startslice,endslice,midslice, ...
  LVcenterx,LVcentery, ...
  manualslices,maxcount,minnbrofslices,maxROIdiameter,samplefactor,segmentationof,silent) %#ok<DEFNU>
%--------------------------------------------------------------------------
%manual confirmation of the LV center point and the angle to RV center

global DATA SET

% %resample the midslice to circumferential coordinates
% slice = midslice-startslice+1;
% sz = size(resampledIM);
% nbrofradius = 80;
% maxradius = floor(min([sz(1)-LVcenterx(slice) sz(2)-LVcentery(slice) LVcenterx(slice) LVcentery(slice)]));
% rad = repmat(linspace(0,maxradius,maxradius)',1,nbrofradius);
% alpha = repmat(linspace(pi/2,5*pi/2,nbrofradius),maxradius,1);
% xi = LVcenterx(slice)+rad.*cos(alpha);  %interpolation x-points for the resampling process
% yi = LVcentery(slice)+rad.*sin(alpha);  %interpolation y-points for the resampling process
% inim = sum(resampledIM(:,:,1,max(startslice,midslice-4):midslice),4);
% inim = interp2(inim,yi,xi);  %the resampled image
% inim = double(inim);
% %use dijkstra method for finding the LV midmural line
% rigidity = 0.02;
% elasticity = 0.02;
% outim = dijkstra(inim, 1, rigidity, elasticity, 2);  %calculate a line through maximum intensity in inim
% switch segmentationof
%   case 'LV'
%     minradius = 20;
%     addradius = 5;
%   case {'MaR','severity','plotbullseye','Ischemia'}
%     minradius = round(20/mean([SET(no).ResolutionX SET(no).ResolutionY]));
%     addradius = round(5/mean([SET(no).ResolutionX SET(no).ResolutionY]));
%   otherwise
%     myfailed('Undefined segmentation object');
%     return;
% end
% meanradius = round(max(minradius,mean(outim(30:50))));
% %find the RV center outside of LV midmural line, defined by a region of low
% %counts to the left of the LV
% RVsearch = inim(meanradius+addradius:min(size(inim,1)-addradius,meanradius+addradius+minradius),30:50);
% f = fspecial('gaussian',5,0.75);
% RVsearchcount = conv2(RVsearch,f,'valid');
% if isempty(RVsearchcount)
%   RVsearchcount = RVsearch;
% end
% [~,tempindex] = min(RVsearchcount(:));
% [RVcenterindexx,RVcenterindexy] = ind2sub(size(RVsearchcount),tempindex);
% RVcenterindexx = RVcenterindexx+meanradius+addradius-1;
% RVcenterindexy = RVcenterindexy+29;
% RVcenterx = xi(RVcenterindexx,RVcenterindexy);
% RVcentery = yi(RVcenterindexx,RVcenterindexy);

slice = midslice-startslice+1;
nbrofradius = 80;
RVcenterx = LVcenterx(slice);
RVcentery = max(2,LVcentery(slice)-60/SET(no).ResolutionY);

%automatically confirmation
SET(no).Spect.RVcenter = [RVcenterx RVcentery];
return;


%open window for confirmation of LV and RV center points
DATA.GUI.Spectconfirmcenter = mygui(['+spect' filesep 'spectconfirmcenter.fig']);
gui = DATA.GUI.Spectconfirmcenter;
%save data for later use
gui.handles.no = no;
gui.handles.originalno = originalno;
gui.handles.im = resampledIM;
gui.handles.startslice = startslice;
gui.handles.endslice = endslice;
gui.handles.midslice = midslice;
gui.handles.LVcenterx = LVcenterx;
gui.handles.LVcentery = LVcentery;
gui.handles.RVcenterx = RVcenterx;
gui.handles.RVcentery = RVcentery;
gui.handles.nbrofradius = nbrofradius;
gui.handles.manualslices = manualslices;
gui.handles.maxcount = maxcount;
gui.handles.minnbrofslices = minnbrofslices;
gui.handles.maxROIdiameter = maxROIdiameter;
gui.handles.samplefactor = samplefactor;
gui.handles.segmentationof = segmentationof;
gui.handles.silent = silent;

%plot the image
gui.handles.image = imagesc(resampledIM(:,:,1,midslice),'parent',gui.handles.axesim);
hold(gui.handles.axesim,'on');
gui.handles.LVcenterpoint = plot(gui.handles.axesim,LVcentery(slice),LVcenterx(slice),'yo');
gui.handles.RVcenterpoint = plot(gui.handles.axesim,RVcentery,RVcenterx,'ro');
hold(gui.handles.axesim,'off');
colormap('spect');
axis(gui.handles.axesim,'off','equal')

%added to be able to run maketest, confirms centerpoint directly
if DATA.Testing %isequal(popfrombuffer('KeyStroke'),'confirm');
  confirmcenter;
else
  %Blocks execution until figure is closed.
  uiwait(gui.fig);
end


%---------------------
function confirmcenter
%---------------------
%callback from spectconfirmcenter.fig

global DATA SET

gui = DATA.GUI.Spectconfirmcenter;

no = gui.handles.no;
originalno = gui.handles.originalno;
im = gui.handles.im;
startslice = gui.handles.startslice;
endslice = gui.handles.endslice;
midslice = gui.handles.midslice;
LVcenterx = gui.handles.LVcenterx;
LVcentery = gui.handles.LVcentery;
RVcenterx = gui.handles.RVcenterx;
RVcentery = gui.handles.RVcentery;
nbrofradius = gui.handles.nbrofradius;
manualslices = gui.handles.manualslices;
maxcount = gui.handles.maxcount;
minnbrofslices = gui.handles.minnbrofslices;
maxROIdiameter = gui.handles.maxROIdiameter;
samplefactor = gui.handles.samplefactor;
segmentationof = gui.handles.segmentationof;
silent = gui.handles.silent;
SET(no).Spect.RVcenter = [RVcenterx RVcentery];

try
  DATA.GUI.Spectconfirmcenter = close(DATA.GUI.Spectconfirmcenter);
catch %#ok<CTCH>
  close(gcf)
end

if isequal(segmentationof,'LV')
  %define endo and epi with spline fitting
  spect.spectlvsegmentation('myocardialborders',no,originalno,im,startslice,endslice,midslice, ...
    LVcenterx,LVcentery,RVcenterx,RVcentery,nbrofradius,manualslices, ...
    maxcount,minnbrofslices,maxROIdiameter,samplefactor,silent);
elseif isequal(segmentationof,'MaR')
  spect.spectmarsegmentation('segmentmar',no,LVcenterx,LVcentery,RVcenterx,RVcentery);
elseif isequal(segmentationof,'severity')
  spect.spectmarsegmentation('findseverity',no,LVcenterx,LVcentery,RVcenterx,RVcentery);
elseif isequal(segmentationof,'plotbullseye')
  spect.spectmarsegmentation('plotbullseyeTPD',no,LVcenterx,LVcentery,RVcenterx,RVcentery);
elseif isequal(segmentationof,'Ischemia')
  %perfusion analysis
else
  myfailed('Error in RV center definition');
  return;
end

%-------------------------
function setmanuallycenter %#ok<DEFNU>
%-------------------------
%callback from spectconfirmcenter.fig

global DATA

gui = DATA.GUI.Spectconfirmcenter;
set(gui.handles.textdescription,'String',...
  'Use the right mouse click to set the RV center point');
%set the callback for the mouseclick in the image
set(get(gui.handles.axesim,'children'),'ButtonDownFcn','spect.specttools(''updatecenter'')');

%--------------------
function updatecenter %#ok<DEFNU>
%--------------------
%update center points

global DATA

gui = DATA.GUI.Spectconfirmcenter;

%catch mouse click
type = get(gui.fig,'SelectionType');
switch type
  case 'normal'  %left button on the mouse = LV center point
%     point = round(get(gui.handles.axesim,'CurrentPoint'));
%     gui.handles.LVcenterx = point(1,2);
%     gui.handles.LVcentery = point(1,1);
%     set(gui.handles.LVcenterpoint,'xdata',gui.handles.LVcenterx);
%     set(gui.handles.LVcenterpoint,'ydata',gui.handles.LVcentery);
  case 'alt'  %right button on the mouse = RV center point
    point = round(get(gui.handles.axesim,'CurrentPoint'));
    gui.handles.RVcenterx = point(1,2);
    gui.handles.RVcentery = point(1,1);
    set(gui.handles.RVcenterpoint,'ydata',gui.handles.RVcenterx);
    set(gui.handles.RVcenterpoint,'xdata',gui.handles.RVcentery);
end


%-------------------------------------------
function no = uniqueimagestack(no,imagetype)
%-------------------------------------------
global DATA
if length(no) > 1
  textno = '';
  for j = 1:length(no)
    textno = sprintf('%s%d ',textno,no(j));
  end
  text = (['Found more than one Perfusion ',imagetype,' image stack. Define which image stack out of ',num2str(textno),'to use.']);
  %Find correct rest image stack
  notemp = inputdlg({text},[imagetype,' NO'],1,{sprintf('%d',no(1))});
  if isempty(notemp)
    myfailed('Invalid image stack.',DATA.GUI.Segment);
    uniqueimagestack(no,imagetype);
  else
    [no,ok] = str2num(notemp{1}); %#ok<ST2NM>
    if not(ok)
      myfailed('Invalid image stack.',DATA.GUI.Segment);
      uniqueimagestack(no,imagetype);
    end
  end
end





% %--------------------------------------------------------------------------
% function meanint = calcmeanint(im,numsectors,nprofiles,pos,endox,endoy,epix,epiy,sz,resolution) %#ok<DEFNU>
% %--------------------------------------------------------------------------
% %Calculates the mean intensity within sectors as a preperation to generate
% %a bullseye
% %Based on the function viability(calctransmuralityarea)
% 
% %numsectors:            number of sectors
% %nprofiles:             number of profiles
% %pos:                   slices to use
% %endox,endoy,epix,epiy: myocardial borders
% %sz:                    size of image stack
% %res:                   resolution of image stack
% 
% if nargin < 6
%   myfailed('Too few input arguments');
%   return;
% end
% 
% numslices = length(pos);
% meanint = zeros(numsectors,numslices);
% 
% %Upsample model
% if not(nprofiles == size(endox,1))
%   [endox,endoy] = segment('resamplemodel',endox,endoy,nprofiles);
%   [epix,epiy] = segment('resamplemodel',epix,epiy,nprofiles);
% end
% 
% tf = 1;
% 
% mask = zeros(sz(1),sz(2));
% for sloop = 1:numslices
%   if not(isnan(endox(1,tf,pos(sloop))))&&...
%       not(isnan(epix(1,tf,pos(sloop))))
% 
%     %Find sectors
%     [meanx,meany,episectors] = findmeaninsectorslice('epi',numsectors, ...
%       endox(:,tf,pos(sloop)),endoy(:,tf,pos(sloop)),epix(:,tf,pos(sloop)),epiy(:,tf,pos(sloop)));
%     [meanx,meany,endosectors] = findmeaninsectorslice('endo',numsectors, ...
%       endox(:,tf,pos(sloop)),endoy(:,tf,pos(sloop)),epix(:,tf,pos(sloop)),epiy(:,tf,pos(sloop)));
% 
%     N = nprofiles;
%     for loop=1:numsectors
%       if episectors(loop)<episectors(loop+1)
%         tempepix = epix(episectors(loop):episectors(loop+1),tf,pos(sloop));
%         tempepiy = epiy(episectors(loop):episectors(loop+1),tf,pos(sloop));
%       else
%         tempepix = [...
%           epix(episectors(loop):N,tf,pos(sloop)) ; ...
%           epix(1:episectors(loop+1),tf,pos(sloop))];
%         tempepiy = [...
%           epiy(episectors(loop):N,tf,pos(sloop)) ; ...
%           epiy(1:episectors(loop+1),tf,pos(sloop))];
%       end;
% 
%       if endosectors(loop)<endosectors(loop+1)
%         tempendox = endox(endosectors(loop):endosectors(loop+1),tf,pos(sloop));
%         tempendoy = endoy(endosectors(loop):endosectors(loop+1),tf,pos(sloop));
%       else
%         tempendox = [...
%           endox(endosectors(loop):N,tf,pos(sloop)) ; ...
%           endox(1:endosectors(loop+1),tf,pos(sloop))];
%         tempendoy = [...
%           endoy(endosectors(loop):N,tf,pos(sloop)) ; ...
%           endoy(1:endosectors(loop+1),tf,pos(sloop))];
%       end;
% 
%       %Create mask image
%       mask = roipoly(mask,...
%         [tempendoy;flipud(tempepiy)],...
%         [tempendox;flipud(tempepix)]);
%       
%       %the same number of points along the lines
%       if length(tempepix) ~= length(tempendox)
%         xi = linspace(0,1,length(tempepix));
%         xx = linspace(0,1,length(tempendox));
%         tempepix = interp1(xi,tempepix,xx)';
%         tempepiy = interp1(xi,tempepiy,xx)';
%       end
%       %calculate the wallthickness fo the current mask
%       tempwallthickness = mean((tempendox-tempepix).^2./resolution(1)+(tempendoy-tempepiy).^2./resolution(2));
%       if sum(sum(mask)) == 0 && tempwallthickness > 2
%         %create a mask of size 1 pixel if the mask is empty and the wall thickness is over 2 mm
%         mask(round(mean([tempendoy;flipud(tempepiy)])),round(mean([tempendox;flipud(tempepix)]))) = 1;
%       end
%       %calculate mean intensity within the sector
%       if sum(mask(:)) == 0
%         meanint(loop,sloop) = NaN;
%       else
%         meanint(loop,sloop) = sum(sum(mask.*im(:,:,pos(sloop))))/(sum(sum(mask)));
%       end
%     end
%   end
% end
% 
% %-------------------------------------------------------------------------------
% function [varargout] = findmeaninsectorslice(type,numsectors,endox,endoy,epix,epiy)
% %-------------------------------------------------------------------------------
% %Find indicies of points along the contour that corresponds to which
% %sector. Used when analysing the myocardium of short axis slices. Operates
% %on given slice.
% 
% if isequal(type,'endo')
%   isendo=1;
% else
%   isendo=0;
% end;
% 
% %Find center
% x = mean(epix);
% y = mean(epiy);
% meanx = mean(x(:));
% meany = mean(y(:));
% 
% %Find sectors
% if isendo
%   ang = angle(complex(...
%     endoy(1:(end-1))-meany,...
%     endox(1:(end-1))-meanx));
% else
%   ang = angle(complex(...
%     epiy(1:(end-1))-meany,...
%     epix(1:(end-1))-meanx));
% end;
% ang = ang*180/pi;
% wantangles = linspace(-180,180,numsectors+1);
% 
% sectors = zeros(numsectors+1,1);
% for loop=1:(numsectors+1)
%   [trash1,inda] = min(abs(ang(:)-wantangles(loop)));
%   [trash2,temp] = min(abs(ang(:)-wantangles(loop)-360));
%   if trash2<trash1
%     inda = temp;
%     trash1 = trash2;
%   end;
%   [trash2,temp] = min(abs(ang(:)-wantangles(loop)+360));
%   if trash2<trash1
%     inda = temp;
%   end;
%   sectors(loop) = inda;
% end;
% 
% varargout = cell(1,nargout-1);
% if nargout>=1
%   varargout{1} = meanx;
% end;
% if nargout>=2
%   varargout{2} = meany;
% end;
% if nargout>=3
%   varargout{3} = sectors;
% end;
