function [setnew,r,myoT] = alignSlides(setold,t0, t1, slices, done)
% Registers images from timeframe t0 to t1 in slices s, 
%inputs
%int t0 initial timeframe.
%int t1 final timeframe
%index array slices, slices to register, set slices=[] to automaticly detect s where epi and endo is present
%
% done logical(number of timeframes x number of slices) to true to not
% register an image, default is false

%Written by Daniel
%
%Modified (slightly) for Segment coding standard and debugging by Einar Heiberg.
  
global SETNEW OPT_PAR
%parameters

% faster and less robust parameter example
OPT_PAR.nControl=10; %desired number of control points, range 2-80 increase it to improve results at the cost of running time 
OPT_PAR.weigthsRadius=10;% was 5; % raidius around epi and endo to register, in pixels. 
OPT_PAR.scale=1; %scales for gaussian low pass filter
OPT_PAR.laplace=10; % smoothing parameter 

%Simulated anealing parameters
OPT_PAR.META.boundTPS=0.1; %parameter bounds
OPT_PAR.META.boundAffine=0.1;
OPT_PAR.META.tMax=1; %maximum time per image in minutes, should not be reached
OPT_PAR.META.TEMP_START=10; %initial temperaure, higher gives more random steps, zero for downhill simplex
OPT_PAR.META.TEMP_END=0;
OPT_PAR.META.MAX_ITER_TOTAL=10; %decrease for faster runtime at the cost of worse results (and increase to improve results)
OPT_PAR.META.MAX_ITER_FIRST=20;
OPT_PAR.META.COOL_RATE=0.95;
OPT_PAR.META.RATIO=.99;
calcResults=false; % get diagnostic results (dice) if true

% % slower and more robust parameter example
% %     
% %parameters
% OPT_PAR.nControl=20; %desired number of control points, range 2-80 increase it to improve results at the cost of running time 
% OPT_PAR.weigthsRadius =5; % raidius around epi and endo to register, in pixels. 
% OPT_PAR.scale=[1]; %scales for gaussian low pass filter
% OPT_PAR.laplace=10; % smoothing parameter
% 
% %Simulated anealing parameters
% OPT_PAR.META.boundTPS=0.1; %parameter bounds
% OPT_PAR.META.boundAffine=0.1;
% OPT_PAR.META.tMax=1;    %maximum time per image in minutes, should not be reached
% OPT_PAR.META.TEMP_START=100; %initial temperaure, higher gives more random steps, zero for downhill simplex
% OPT_PAR.META.TEMP_END=0;
% OPT_PAR.META.MAX_ITER_TOTAL=30; %decrease for faster runtime at the cost of worse results (and increase to improve results)
% OPT_PAR.META.MAX_ITER_FIRST=20;
% OPT_PAR.META.COOL_RATE=0.95;
% OPT_PAR.META.RATIO=.99;
% calcResults=false; % get diagnostic results (dice) if true

SETNEW = setold;
if nargin<5
    done=logical(false(SETNEW.TSize,SETNEW.ZSize)); %default, register all images in range
end

if nargin<4
    slices=[];
end

OPT_PAR.sz = size(SETNEW.IM);

%Make sure it is a 4 element vector
if length(OPT_PAR.sz)<4
  OPT_PAR.sz = [OPT_PAR.sz 1 1];
  OPT_PAR.sz = OPT_PAR.sz(1:4);
end

v = ~isnan(SETNEW.EpiX); %true in images where epi is present 
if isempty(slices)
      slices=find(sum(sum(v,2),1)); %if slice not specified set where endo/epi is availible
end
    try 
    referenceFrame=SETNEW.referenceFrame; %reference frame already specified, ie, we are running again to refine results or add more timeframes
    %[OPT_PAR.k,OPT_PAR.z,referenceFrame,slices]=init(referenceFrame,s); %% initilizes weigths and control points
    catch
   
        for i=1:size(slices)
            referenceFrame(slices(i))=find(sum(v(:,:,slices(i)),1),1,'first');
        end
    %[OPT_PAR.k,OPT_PAR.z,referenceFrame,slices]=init([],s); %% initilizes weigths and control points, sets reference frame if its not specified 
    end

r=zeros(OPT_PAR.sz(3),OPT_PAR.sz(4),4);
SETNEW.referenceFrame=referenceFrame;
wbstep = 1/(length(slices)*length(t0:t1));
wbval = 0;
h = mywaitbarstart(length(slices)*length(t0:t1),dprintf('Generating image stacks'));
    for i=1:length(slices)
     s=slices(i);
     [OPT_PAR.k,OPT_PAR.z]=init(referenceFrame(s),s); 
     for t=t0:t1
        a=[0 0;eye(2)];
        c=zeros(size(OPT_PAR.z,1),2);
        param=[a(:); c(:)]; %initial parameters
        
        
            if ~(done(t,s) || (t==referenceFrame(s)))
                tic;
                [SETNEW.IM(:,:,t,s),param,r(t,s,1)]=register(SETNEW.IM(:,:,referenceFrame(s),s),SETNEW.IM(:,:,t,s), param,referenceFrame(s),s);
                
                if calcResults && ~isnan(SETNEW.EpiY(1,t,s)) 
                    
                    [r(t,s,2:3),myoT(t-t0+1,:,:)]=dice(param,OPT_PAR.z,referenceFrame(s),t,s); 
                    
                end
               
                r(t,s,4)=toc; % time per image
                SETNEW.EpiY(:,t,s)=SETNEW.EpiY(:,referenceFrame(s),s); 
                SETNEW.EpiX(:,t,s)=SETNEW.EpiX(:,referenceFrame(s),s);
                SETNEW.EndoY(:,t,s)=SETNEW.EndoY(:,referenceFrame(s),s);
                SETNEW.EndoX(:,t,s)=SETNEW.EndoX(:,referenceFrame(s),s);
                wbval = wbval + wbstep;
%                 waitbar(wbval,h);           
                h = mywaitbarupdate(h);

            end
      end 

   end
          
        

    setnew=SETNEW; %return
    clear SETNEW
    mywaitbarclose(h);
%     close(h);
    clear OPT_PAR
     

    function [dice,myoT,myoR]=dice(param,z,baseframe,t,s,sr)
        % dice similarity measure, myocardium overlap coefficient for
        % registred and unregistred images, requires endo/epi in t as well as baseframe. 
        % probably only useful to developer
        % not used in registration
        global SETNEW OPT_PAR
        sz=OPT_PAR.sz;
        A=reshape(param(1:6),3,2);
        A=[[1 0 0]'  A];
        sr=2;%sampling rate 
        c=reshape(param(7:end),(length(param)-6)/2,2);
        EpiR=[SETNEW.EpiY(:,baseframe,s) SETNEW.EpiX(:,baseframe,s)];
        EndoR=[SETNEW.EndoY(:,baseframe,s) SETNEW.EndoX(:,baseframe,s)];
        Epi=[SETNEW.EpiY(:,t,s) SETNEW.EpiX(:,t,s)];
        Endo=[SETNEW.EndoY(:,t,s) SETNEW.EndoX(:,t,s)];
        [X,Y]=meshgrid(1:1:sz(2)*sr,1:1:sz(1)*sr);
        myoR=logical(false(sz(1)*sr,sz(2)*sr));
        [IN] = inpolygon(X,Y,EpiR(:,1)*sr,EpiR(:,2)*sr);
        myoR(IN)=true;
        [IN] = inpolygon(X,Y,EndoR(:,1)*sr,EndoR(:,2)*sr);
        myoR(IN)=false;
        
        myo=logical(false(sz(1)*sr,sz(2)*sr));
        [IN] = inpolygon(X,Y,Epi(:,1)*sr,Epi(:,2)*sr);
        myo(IN)=true;
        
        [IN] = inpolygon(X,Y,Endo(:,1)*sr,Endo(:,2)*sr);
        myo(IN)=false;
    
        myoT=transformIm(myo,A,c,[],z, 'final')>0.5;
        dice(1)=2*(sum(myoT(:) .* myoR(:)))/(sum(myoT(:))+sum(myoR(:)));

        dice(2)=2*(sum(myo(:) .* myoR(:)))/(sum(myo(:))+sum(myoR(:)));
   
 %-----------------------------------------------
function perf= perfusion(t,s)
%-----------------------------------------------
%not used for registration
global SETNEW OPT_PAR
imt=SETNEW.IM(:,:,t,s); 

endopct = (100-20)/100;
epipct = 80/100;

if endopct < 0 || endopct > 1 || epipct < 0 || epipct > 1 || (1-endopct) > epipct
  myfailed('Invalid input of percentage')
end



  
  %Check if valid segmentation exists
  if not(isnan(SETNEW.EndoX(1,t,s)))&&not(isnan(SETNEW.EpiX(1,t,s)))
    
    %Make sure segmentation is correct
    segment('checkconsistency',t,s);
    
    endox = SETNEW.EndoX(:,t,s);
    endoy = SETNEW.EndoY(:,t,s);
    epix = SETNEW.EpiX(:,t,s);
    epiy = SETNEW.EpiY(:,t,s);
    xin = epix+(endox-epix)*(1-endopct);
    yin = epiy+(endoy-epiy)*(1-endopct);
    xout = epix+(endox-epix)*epipct;
    yout = epiy+(endoy-epiy)*epipct;
    
  end %valid segmentation exist
  [X,Y]=meshgrid(1:1:OPT_PAR.sz(1)*1,1:1:OPT_PAR.sz(2)*1);
  myoR=logical(false(OPT_PAR.sz(1),OPT_PAR.sz(2)));
  [IN] = inpolygon(X,Y,yin,xin);
  myoR(IN)=true;
  [IN] = inpolygon(X,Y,yout,xout);
  myoR(IN)=false;

  
  perf(1)=mean(mean(imt(myoR)));
  mx=mean(endox);
  my=mean(endoy);
  mrad=mean(sqrt((endox-mx).^2+(endoy-my).^2));
  theta=1:pi/40:2*pi;
  %halvLV=mrad/2*[sin(theta) cos(theta)];
  LV=logical(false(size(myoR)));
  LV(inpolygon(X,Y,mrad/2*sin(theta)+my,mrad/2*cos(theta)+mx))=true;
  perf(2)=mean(mean(imt(LV)));
  
  
  

  

        

        


function curve=transformC(X,A,c,z) 
    %transforms a curve, 
    % not used for registration
    global OPT_PAR;
    
    maX=max(max(OPT_PAR.sz(1:2))); 
    % scale coordinates
    X=X./maX-0.5;    
    K=ctps_gen_K(X,z);
    X=[ones(size(X,1),1) X];
    U=X*A+K*[zeros(1,size(c,1))' c];
    U=[U(:,2)./U(:,1) U(:,3)./U(:,1)];
    %descale
    curve=(U+0.5).*maX;

function [k,z,ref,s]=init(ref,s)
    %initializes constant parameters before optimization
    
    global SETNEW OPT_PAR
    OPT_PAR.sz=size(SETNEW.IM);
    Epi=[SETNEW.EpiY(:,ref,s) SETNEW.EpiX(:,ref,s)];
    Endo=[SETNEW.EndoY(:,ref,s) SETNEW.EndoX(:,ref,s)];
    w = zeros(OPT_PAR.sz(1:2));
    if isempty(SETNEW.EpiX)
        %this should never happen, if it did anyway someone probably edited the
        %main alignslides function
        disp('Warning:No Epi/Endo in slice')
        kant=round(0.2*sz(1:2));
        w(kant(1):end-kant(1),kant(2):end-kant(2))=1;
    else 
        %set weigths based on endo/epi
        ind=[ [round(SETNEW.EpiX(:,ref,s)) ; round(SETNEW.EndoX(:,ref,s))] [round(SETNEW.EpiY(:,ref,s)); round(SETNEW.EndoY(:,ref,s))]];
    
    
        OPT_PAR.weightsRadius=5; 
    for i =1:length(ind)
        x=ind(i,1);
        y=ind(i,2);
        w(min(OPT_PAR.sz(1),max(1,x-OPT_PAR.weightsRadius)):min(OPT_PAR.sz(1),max(1,x+OPT_PAR.weightsRadius)),min(OPT_PAR.sz(2),max(1,y-OPT_PAR.weightsRadius)):min(OPT_PAR.sz(2),max(1,y+OPT_PAR.weightsRadius)))=1; % check if in bounds and set weigths
    end
   rectanguarGrid=false;
   end
   
   if rectanguarGrid 
       %rectangular control point grid, not currently used 
        [col,row] = find(w);
        R2=10;
        bounds=[min(row)-R2 max(row)+R2;min(col)-R2 max(col)+R2];
        nControl=round(sqrt(OPT_PAR.nControl));     %total number of control points is nControl^2
        gridsize=(bounds(:,2)-bounds(:,1))/nControl;      
        [X,Y] = meshgrid(gridsize(1)/2:gridsize(1):gridsize(1)*(nControl-0.5),gridsize(1)/2:gridsize(1):gridsize(1)*(nControl-0.5));
        
        z=[X(:)+bounds(1,1) Y(:)+bounds(2,1)];  %control points
    
    else % 
       % circular control point grid, replace with Helen's code
      %%% Use code from Helen (centerline method). 
       z=Endo/2+Epi/2;
       z=z(1:round(length(z)/OPT_PAR.nControl):length(z),:);   
       %z=[z; mean(z)]; 
    end
w=w(:);
k=w>0;
[x,y]=meshgrid(1:OPT_PAR.sz(2),1:OPT_PAR.sz(1));
X=[x(:) y(:)];
maX=max(max(OPT_PAR.sz(1:2)));
OPT_PAR.X=X(k,:);
OPT_PAR.x=x;
OPT_PAR.y=y;
OPT_PAR.X=(OPT_PAR.X)./maX-0.5; %scale image points
z=z./maX-0.5;%scale control points



function  [Ir,param,r] = register(It, Ib, param,t0,s)
    
    global OPT_PAR;
    
    META=OPT_PAR.META;
    a=[0 0;eye(2)];
    c=zeros(size(OPT_PAR.z,1),2);
    param0=[a(:); c(:)]; %initial parameters
    
    A=eye(3,3);
    [~,OPT_PAR.K]=transformIm(It,A,c,[],OPT_PAR.z, 'set_K'); %this will set K
    for i =1:length(OPT_PAR.scale) %only one scale currently used
        h = fspecial('gaussian',[21 21], OPT_PAR.scale(i));
        IT = imfilter(It , h,'replicate'); 
        IB = imfilter(Ib, h,'replicate');
        OPT_PAR.IT=reshape(modNGF(IT),OPT_PAR.sz(1:2));
        OPT_PAR.IT=single((OPT_PAR.IT(OPT_PAR.k)));
        OPT_PAR.IB=single(reshape(modNGF(IB),OPT_PAR.sz(1:2)));
        alfa=META.boundTPS*ones(size(param0));
        alfa(1:6)=META.boundAffine;
        ub=param0+alfa; 
        lb=param0-alfa;
        f0=cost_SIMPSA(param0);
        options = SIMPSASET('MAX_TIME',1,'TEMP_START',META.TEMP_START,'TEMP_END',0,'MAX_ITER_TOTAL', META.MAX_ITER_TOTAL,'MAX_ITER_TEMP_FIRST', META.MAX_ITER_FIRST,'COOL_RATE',META.COOL_RATE,'TOLFUN',0.001,'DISPLAY','none','TOLX', 10^-3,'INITIAL_ACCEPTANCE_RATIO',META.RATIO);
        [param_final,FVAL,EXITFLAG]= SIMPSA('cost_SIMPSA',param0,lb,ub,options);%update parameters for each scale space/pyramid
        f1=cost_SIMPSA(param_final);
        r=f1/f0; %return for diagnostics
    end
    param=param_final; %parameters to return
    A=reshape(param(1:6),3,2);
    A=[[1 0 0]'  A];
    c=reshape(param(7:end),(length(param)-6)/2,2);
    [Ir]=transformIm(Ib,A,c,[],OPT_PAR.z, 'final');

    
    
function [K] =  ctps_gen_K(x,z)

% calulates TPS kernel matrix
[n, M] = size (x); 
[m, N] = size (z);
dim    = M;

% calc. the K matrix.
% 2D: K = r^2 * log r
% 3D: K = -r
K= zeros (n,m);

for it_dim=1:dim
  tmp = x(:,it_dim) * ones(1,m) - ones(n,1) * z(:,it_dim)';
  tmp = tmp .* tmp;
  K = K + tmp;
  
end
if dim == 2
      mask = K < 1e-10; % to avoid singularity.
      K = 0.5*K .* log(K + mask) .* (K>1e-10);
else
  K = - sqrt(K);
end


function [I1,K,reg,motion_field]=transformIm(I0,A,c,K,z, method)
    %image transformation I1=T(I0) 
    
    %inputs
    %I0= image to be transformed 
    %A = affine transformation matrix 3 x 3
    %c = TPS weigths  2 x (number of control points)
    %K =TPS kernel (number of control points) x (number of pixels)
    %K =[] sets K, do not re-set K if unless pixels or control points
    %change
    %z = TPS control points
    % method = option string
        %set_K do not transform, sets K
        
        % "opt" faster using bilinear interpolation for optimization
        %     only transforms points w>0 K should already be set
        
        % "final" whole image transformation 
   %outputs
   % I1 =transformed image
   % K =TPS kernel should only be set with "set_K" or "final"
   % reg = laplacian regularization only set with "opt"
   % motion_field only set by "final"
    
    
    global OPT_PAR
    
    I1=zeros(size(I0));

    maX=max(max(OPT_PAR.sz(1:2)));
    % scale coordinates
    
    if isempty(K)
        
        [K] = ctps_gen_K(OPT_PAR.X,z);
        
    end
    if strcmp('set_K', method) %do nothing, we only want to set K
%       


    elseif strcmp('opt', method)   %faster used for optimization
        
        U= inverseFunc(OPT_PAR.X,A,c,K);
       
        reg=sum(sum((OPT_PAR.X-U).^2)); 
        %descale
        U=(U+0.5).*maX;
        
        I1 = interp2(OPT_PAR.x,OPT_PAR.y,I0,U(:,1), U(:,2));

    elseif strcmp('final', method) % slower, used only for final transformation
        [x,y]=meshgrid(1:size(I0,2),1:size(I0,1));
        maX=max(size(I0));
        X=[x(:) y(:)];
        X=X./maX-0.5;
        [K] = ctps_gen_K(X,z);
        U= inverseFunc(X,A,c,K);
        %descale
        U=(U+0.5).*maX;
        motion_field=(U-X);
        I1 = interp2(x,y,I0,U(:,1), U(:,2));
        I1(isnan(I1))=0;
        I1=reshape(I1,size(I0));
    end
    


function U= inverseFunc(X,A,c,K)
    %inverse TPS transformation
    
    X=[ones(size(X,1),1) X];
    U=(X-K*[zeros(1,size(c,1))' c])/(A);
    %U=U(:,1:2);
    U=[U(:,2)./U(:,1) U(:,3)./U(:,1)];
    %diff=U-U1;


function E= cost_SIMPSA(param)
    %cost function
global  OPT_PAR
    J=OPT_PAR.laplace;
    A=reshape(param(1:6),3,2);
    A=[[1 0 0]'  A];
    c=reshape(param(7:end),(length(param)-6)/2,2);
    [Itrans,~,reg]= transformIm(OPT_PAR.IB,A,c,OPT_PAR.K,OPT_PAR.z, 'opt');
    E=sum(sum((OPT_PAR.IT-Itrans).^2))+reg*J;
        





    
function G=modNGF(I,noise)
    % calculates amplitude of modified normalized gradient field
    [FX,FY] = gradient((I));
    
    n=numel(I);
    gradnorm2=((FX.^2+FY.^2));
    noise_diam=15;
    if nargin<2
        noise=noise_est((gradnorm2),ones(noise_diam)); %local noise estimation
        %noise=std((gradnorm2(:))); %global noise estimation only for test
    end
    eps= noise+sum(sum(gradnorm2))/(noise_diam^2);
    adjnorm=(gradnorm2+eps.^2).^(1/2);
    x=FX./adjnorm;
    y=FY./adjnorm;
    G=x.^2+y.^2;
    
function J = noise_est(varargin)
% modified version of matlabs local noise estimation, classifies edges as
% outliers
%STDFILT Local standard deviation of image.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision.2 $  $Date: 2006/06/15 20:09:36 $

[I, h] = ParseInputs(varargin{:});

if (~isa(I,'double'))
    I = double(I);
end


n = sum(h(:));

% If n = 1 then return default J (all zeros) to avoid the divideByZero warning.
% Otherwise, calculate standard deviation. The formula for standard deviation
% can be rewritten in terms of the theoretical definition of
% convolution. However, in practise, use correlation in IMFILTER to avoid a
% flipped answer when NHOOD is asymmetric.
% conv1 = imfilter(I.^2,h,'symmetric') / (n-1); 
% conv2 = imfilter(I,h,'symmetric').^2 / (n*(n-1));
% std = sqrt(conv1-conv2).  
% These equations can be further optimized for speed.
BW = edge(I,'canny');
%BW=zeros(size(BW));
n=imfilter(double(~BW), h , 'symmetric');
n1 = n - 1;
if n ~= 1
  conv1 = imfilter(I.^2, h , 'symmetric')./(n.*n1);
  conv2 = imfilter(I, h, 'symmetric').^2 ./ (n.*n1);
  J = sqrt(max((conv1 - conv2),0));
else
  J = zeros(size(I));
end

%%%%%%%%%%%%%%%ParseInputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,H] = ParseInputs(varargin)

narginchk(1,2);

validateattributes(varargin{1},{'numeric','logical'},{'real','nonsparse'}, ...
              mfilename, 'I',1);
I = varargin{1};

if nargin == 2
  validateattributes(varargin{2},{'logical','numeric'},{'nonsparse'}, ...
                mfilename,'NHOOD',2);
  H = varargin{2};
  
  eid = sprintf('Images:%s:invalidNeighborhood',mfilename);
  
  % H must contain zeros and/or ones.
  bad_elements = (H ~= 0) & (H ~= 1);
  if any(bad_elements(:))
    msg = 'NHOOD must be a matrix that contains zeros and/or ones.';
    error(eid,'%s',msg);
  end
  
  % H's size must be a factor of 2n-1 (odd).
  sizeH = size(H);
  if any(floor(sizeH/2) == (sizeH/2) )
    msg = 'NHOOD must have a size that is odd in each dimension.';
    error(eid,'%s',msg);
  end

  if ~isa(H,'double')
    H = double(H);
  end

else
  H = ones(3);
end


    
function [X,FVAL,EXITFLAG,OUTPUT] = SIMPSA(FUN,X0,LB,UB,OPTIONS,varargin)
%SIMPSA finds a minimum of a function of several variables using an algorithm 
% that is based on the combination of the non-linear smplex and the simulated 
% annealing algorithm (the SIMPSA algorithm, Cardoso et al., 1996). 
% In this paper, the algorithm is shown to be adequate for the global optimi-
% zation of an example SETNEW of unconstrained and constrained NLP functions.
%
%   SIMPSA attempts to solve problems of the form:
%       min F(X) subject to: LB <= X <= UB
%        X
%
% Algorithm partly based on section 10.4 and 10.9 in "Numerical Recipes in C",
% ISBN 0-521-43108-5, and the paper of Cardoso et al, 1996.
%                                                                             
%   X=SIMPSA(FUN,X0) start at X0 and finds a minimum X to the function FUN. 
%   FUN accepts input X and returns a scalar function value F evaluated at X.
%   X0 may be a scalar, vector, or matrix.
%   
%   X=SIMPSA(FUN,X0,LB,UB) defines a SETNEW of lower and upper bounds on the 
%   design variables, X, so that a solution is found in the range 
%   LB <= X <= UB. Use empty matrices for LB and UB if no bounds exist. 
%   SETNEW LB(i) = -Inf if X(i) is unbounded below; SETNEW UB(i) = Inf if X(i) is 
%   unbounded above.
%   
%   X=SIMPSA(FUN,X0,LB,UB,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument 
%   created with the SIMPSASETNEW function. See SIMPSASETNEW for details. 
%   Used options are TEMP_START, TEMP_END, COOL_RATE, INITIAL_ACCEPTANCE_RATIO,
%   MIN_COOLING_FACTOR, MAX_ITER_TEMP_FIRST, MAX_ITER_TEMP_LAST, MAX_ITER_TEMP,
%   MAX_ITER_TOTAL, MAX_TIME, MAX_FUN_EVALS, TOLX, TOLFUN, DISPLAY and OUTPUT_FCN.
%   Use OPTIONS = [] as a place holder if no options are SETNEW.
%   
%   X=SIMPSA(FUN,X0,LB,UB,OPTIONS,varargin) is used to supply a variable 
%   number of input arguments to the objective function FUN.
%   
%   [X,FVAL]=SIMPSA(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%   
%   [X,FVAL,EXITFLAG]=SIMPSA(FUN,X0,...) returns an EXITFLAG that describes the 
%   exit condition of SIMPSA. Possible values of EXITFLAG and the corresponding 
%   exit conditions are:
%   
%     1  Change in the objective function value less than the specified tolerance.
%     2  Change in X less than the specified tolerance.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Maximum time exceeded.
%   
%   [X,FVAL,EXITFLAG,OUTPUT]=SIMPSA(FUN,X0,...) returns a structure OUTPUT with 
%   the number of iterations taken in OUTPUT.nITERATIONS, the number of function
%   evaluations in OUTPUT.nFUN_EVALS, the temperature profile in OUTPUT.TEMPERATURE,
%   the simplexes that were evaluated in OUTPUT.SIMPLEX and the best one in 
%   OUTPUT.SIMPLEX_BEST, the costs associated with each simplex in OUTPUT.COSTS and 
%   from the best simplex at that iteration in OUTPUT.COST_BEST, the amount of time 
%   needed in OUTPUT.TIME and the options used in OUTPUT.OPTIONS.
% 
%   See also SIMPSASETNEW, SIMPSAGET



% Copyright (C) 2006 Brecht Donckels, BIOMATH, brecht.donckels@ugent.be
% 
% inspired by:
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005 Henning Schmidt, FCC, henning@fcc.chalmers.se
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.



% handle variable input arguments

if nargin < 5
    OPTIONS = [];
    if nargin < 4
        UB = 1e5;
        if nargin < 3
            LB = -1e5;
        end
    end
end

% check input arguments

if ~ischar(FUN)
    error('''FUN'' incorrectly specified in ''SIMPSA''');
end
if ~isfloat(X0)
    error('''X0'' incorrectly specified in ''SIMPSA''');
end
if ~isfloat(LB)
    error('''LB'' incorrectly specified in ''SIMPSA''');
end
if ~isfloat(UB)
    error('''UB'' incorrectly specified in ''SIMPSA''');
end
if length(X0) ~= length(LB)
    error('''LB'' and ''X0'' have incompatible dimensions in ''SIMPSA''');
end
if length(X0) ~= length(UB)
    error('''UB'' and ''X0'' have incompatible dimensions in ''SIMPSA''');
end

% declaration of global variables

global NDIM nFUN_EVALS TEMP YBEST PBEST

% SETNEW EXITFLAG to default value

EXITFLAG = -2;

% determine number of variables to be optimized

NDIM = length(X0);

% seed the random number generator

rand('state',sum(100*clock));

% SETNEW default options

DEFAULT_OPTIONS = SIMPSASET('TEMP_START',[],...  % starting temperature (if none provided, an optimal one will be estimated)
             'TEMP_END',1,...                    % end temperature
             'COOL_RATE',10,...                  % small values (<1) means slow convergence,large values (>1) means fast convergence
             'INITIAL_ACCEPTANCE_RATIO',0.95,... % when initial temperature is estimated, this will be the initial acceptance ratio in the first round
             'MIN_COOLING_FACTOR',0.9,...        % minimum cooling factor (<1)
             'MAX_ITER_TEMP_FIRST',50,...        % number of iterations in the preliminary temperature loop
             'MAX_ITER_TEMP_LAST',50,...         % number of iterations in the last temperature loop (pure simplex)
             'MAX_ITER_TEMP',10,...              % number of iterations in the remaining temperature loops
             'MAX_ITER_TOTAL',2500,...           % maximum number of iterations tout court
             'MAX_TIME',2500,...                 % maximum duration of optimization
             'MAX_FUN_EVALS',2500,...            % maximum number of function evaluations
             'TOLX',1e-6,...                     % maximum difference between best and worst function evaluation in simplex
             'TOLFUN',1e-3,...                   % maximum difference between the coordinates of the vertices
             'DISPLAY','none',...                % 'iter' or 'none' indicating whether user wants feedback
             'OUTPUT_FCN',[]);                   % string with output function name

% update default options with supplied options

OPTIONS = SIMPSASET(DEFAULT_OPTIONS,OPTIONS);

% store options in OUTPUT

OUTPUT.OPTIONS = OPTIONS;

% initialize simplex
% ------------------

% create empty simplex matrix p (location of vertex i in row i)
P = zeros(NDIM+1,NDIM);
% create empty cost vector (cost of vertex i in row i)
Y = zeros(NDIM+1,1);
% SETNEW best vertex of initial simplex equal to initial parameter guess
PBEST = X0(:)';
% calculate cost with best vertex of initial simplex
YBEST = CALCULATE_COST(FUN,PBEST,LB,UB,varargin{:});

% initialize temperature loop
% ---------------------------

% SETNEW temperature loop number to one
TEMP_LOOP_NUMBER = 1;

% if no TEMP_START is supplied, the initial temperature is estimated in the first
% loop as described by Cardoso et al., 1996 (recommended)

% therefore, the temperature is SETNEW to YBEST*1e5 in the first loop
if isempty(OPTIONS.TEMP_START)
    TEMP = abs(YBEST)*1e5;
else
    TEMP = OPTIONS.TEMP_START;
end

% initialize OUTPUT structure
% ---------------------------

OUTPUT.TEMPERATURE = zeros(OPTIONS.MAX_ITER_TOTAL,1);
OUTPUT.SIMPLEX = zeros(NDIM+1,NDIM,OPTIONS.MAX_ITER_TOTAL);
OUTPUT.SIMPLEX_BEST = zeros(OPTIONS.MAX_ITER_TOTAL,NDIM);
OUTPUT.COSTS = zeros(OPTIONS.MAX_ITER_TOTAL,NDIM+1);
OUTPUT.COST_BEST = zeros(OPTIONS.MAX_ITER_TOTAL,1);

% initialize iteration data
% -------------------------

% start timer
tic
% SETNEW number of function evaluations to one
nFUN_EVALS = 1;
% SETNEW number of iterations to zero
nITERATIONS = 0;

% temperature loop: run SIMPSA till stopping criterion is met
% -----------------------------------------------------------

while 1
    
    % detect if termination criterium was met
    % ---------------------------------------
    
    % if a termination criterium was met, the value of EXITFLAG should have changed
    % from its default value of -2 to -1, 0, 1 or 2
    
    if EXITFLAG ~= -2
        break
    end
    
    % SETNEW MAXITERTEMP: maximum number of iterations at current temperature
    % --------------------------------------------------------------------
    
    if TEMP_LOOP_NUMBER == 1
        MAXITERTEMP = OPTIONS.MAX_ITER_TEMP_FIRST*NDIM;
        % The initial temperature is estimated (is requested) as described in 
        % Cardoso et al. (1996). Therefore, we need to store the number of 
        % successful and unsuccessful moves, as well as the increase in cost 
        % for the unsuccessful moves.
        if isempty(OPTIONS.TEMP_START)
            [SUCCESSFUL_MOVES,UNSUCCESSFUL_MOVES,UNSUCCESSFUL_COSTS] = deal(0);
        end
    elseif TEMP < OPTIONS.TEMP_END
        TEMP = 0;
        MAXITERTEMP = OPTIONS.MAX_ITER_TEMP_LAST*NDIM;
    else
        MAXITERTEMP = OPTIONS.MAX_ITER_TEMP*NDIM;
    end
    
    % construct initial simplex
    % -------------------------
    
    % 1st vertex of initial simplex
    P(1,:) = PBEST;
    Y(1) = CALCULATE_COST(FUN,P(1,:),LB,UB,varargin{:});
    
    % if output function given then run output function to plot intermediate result
    if ~isempty(OPTIONS.OUTPUT_FCN)
        feval(OPTIONS.OUTPUT_FCN,P(1,:),Y(1));
    end
    
    % remaining vertices of simplex
    for k = 1:NDIM
        % copy first vertex in new vertex
        P(k+1,:) = P(1,:);
        % alter new vertex
        P(k+1,k) = LB(k)+rand*(UB(k)-LB(k));
        % calculate value of objective function at new vertex
        Y(k+1) = CALCULATE_COST(FUN,P(k+1,:),LB,UB,varargin{:});
    end
    
    % store information on what step the algorithm just did
    ALGOSTEP = 'initial simplex';
    
    % add NDIM+1 to number of function evaluations
    nFUN_EVALS = nFUN_EVALS + NDIM;
    
    % note:
    %  dimensions of matrix P: (NDIM+1) x NDIM
    %  dimensions of vector Y: (NDIM+1) x 1
    
    % give user feedback if requested
    if strcmp(OPTIONS.DISPLAY,'iter')
        if nITERATIONS == 0
            disp(' Nr Iter  Nr Fun Eval    Min function       Best function        TEMP           Algorithm Step');
        else
            disp(sprintf('%5.0f      %5.0f       %12.6g     %15.6g      %12.6g       %s',nITERATIONS,nFUN_EVALS,Y(1),YBEST,TEMP,'best point'));
        end
    end

    % run full metropolis cycle at current temperature
    % ------------------------------------------------
    
    % initialize vector COSTS, needed to calculate new temperature using cooling
    % schedule as described by Cardoso et al. (1996)
    COSTS = zeros((NDIM+1)*MAXITERTEMP,1);
    
    % initialize ITERTEMP to zero
    
    ITERTEMP = 0;
    
    % start

    for ITERTEMP = 1:MAXITERTEMP
        
        % add one to number of iterations
        nITERATIONS = nITERATIONS + 1;
        
        % Press and Teukolsky (1991) add a positive logarithmic distributed variable,
        % proportional to the control temperature T to the function value associated with 
        % every vertex of the simplex. Likewise,they subtract a similar random variable 
        % from the function value at every new replacement point.
        % Thus, if the replacement point corresponds to a lower cost, this method always
        % accepts a true down hill step. If, on the other hand, the replacement point 
        % corresponds to a higher cost, an uphill move may be accepted, depending on the
        % relative COSTS of the perturbed values.
        % (taken from Cardoso et al.,1996)
        
        % add random fluctuations to function values of current vertices
        YFLUCT = Y+TEMP*abs(log(rand(NDIM+1,1)));
        
        % reorder YFLUCT, Y and P so that the first row corresponds to the lowest YFLUCT value
        help = sortrows([YFLUCT,Y,P],1);
        YFLUCT = help(:,1);
        Y = help(:,2);
        P = help(:,3:end);
        
        % store temperature at current iteration
        OUTPUT.TEMPERATURE(nITERATIONS) = TEMP;

        % store information about simplex at the current iteration
        OUTPUT.SIMPLEX(:,:,nITERATIONS) = P;
        OUTPUT.SIMPLEX_BEST(nITERATIONS,:) = PBEST;
        
        % store cost function value of best vertex in current iteration
        OUTPUT.COSTS(nITERATIONS,:) = Y;
        OUTPUT.COST_BEST(nITERATIONS) = YBEST;
        
        if strcmp(OPTIONS.DISPLAY,'iter')
            disp(sprintf('%5.0f      %5.0f       %12.6g     %15.6g      %12.6g       %s',nITERATIONS,nFUN_EVALS,Y(1),YBEST,TEMP,ALGOSTEP));
        end
        
        % if output function given then run output function to plot intermediate result
        if ~isempty(OPTIONS.OUTPUT_FCN)
            feval(OPTIONS.OUTPUT_FCN,P(1,:),Y(1));
        end
        
        % end the optimization if one of the stopping criteria is met
        %% 1. difference between best and worst function evaluation in simplex is smaller than TOLFUN 
        %% 2. maximum difference between the coordinates of the vertices in simplex is less than TOLX
        %% 3. no convergence,but maximum number of iterations has been reached
        %% 4. no convergence,but maximum time has been reached
            
        if (abs(max(Y)-min(Y)) < OPTIONS.TOLFUN) && (TEMP_LOOP_NUMBER ~= 1)
            if strcmp(OPTIONS.DISPLAY,'iter')
                disp('Change in the objective function value less than the specified tolerance (TOLFUN).')
            end
            EXITFLAG = 1;
            break;
        end
        
        if (max(max(abs(P(2:NDIM+1,:)-P(1:NDIM,:)))) < OPTIONS.TOLX) && (TEMP_LOOP_NUMBER ~= 1)
            if strcmp(OPTIONS.DISPLAY,'iter')
                disp('Change in X less than the specified tolerance (TOLX).')
            end
            EXITFLAG = 2;
            break;
        end
        
        if (nITERATIONS >= OPTIONS.MAX_ITER_TOTAL*NDIM) || (nFUN_EVALS >= OPTIONS.MAX_FUN_EVALS*NDIM*(NDIM+1))
            if strcmp(OPTIONS.DISPLAY,'iter')
                disp('Maximum number of function evaluations or iterations reached.');
            end
            EXITFLAG = 0;
            break;
        end
        
        if toc/60 > OPTIONS.MAX_TIME
            if strcmp(OPTIONS.DISPLAY,'iter')
                disp('Exceeded maximum time.');
            end
            EXITFLAG = -1;
            break;
        end
        
        % begin a new iteration
        
        %% first extrapolate by a factor -1 through the face of the simplex
        %% across from the high point,i.e.,reflect the simplex from the high point
        [YFTRY,YTRY,PTRY] = AMOTRY(FUN,P,-1,LB,UB,varargin{:});
        
        %% check the result
        if YFTRY <= YFLUCT(1)
            %% gives a result better than the best point,so try an additional
            %% extrapolation by a factor 2
            [YFTRYEXP,YTRYEXP,PTRYEXP] = AMOTRY(FUN,P,-2,LB,UB,varargin{:});
            if YFTRYEXP < YFTRY
                P(end,:) = PTRYEXP;
                Y(end) = YTRYEXP;
                ALGOSTEP = 'reflection and expansion';
            else
                P(end,:) = PTRY;
                Y(end) = YTRY;
                ALGOSTEP = 'reflection';
            end
        elseif YFTRY >= YFLUCT(NDIM)
            %% the reflected point is worse than the second-highest, so look
            %% for an intermediate lower point, i.e., do a one-dimensional
            %% contraction
            [YFTRYCONTR,YTRYCONTR,PTRYCONTR] = AMOTRY(FUN,P,-0.5,LB,UB,varargin{:});
            if YFTRYCONTR < YFLUCT(end)
                P(end,:) = PTRYCONTR;
                Y(end) = YTRYCONTR;
                ALGOSTEP = 'one dimensional contraction';
            else
                %% can't seem to get rid of that high point, so better contract
                %% around the lowest (best) point
                X = ones(NDIM,NDIM)*diag(P(1,:));
                P(2:end,:) = 0.5*(P(2:end,:)+X);
                for k=2:NDIM
                    Y(k) = CALCULATE_COST(FUN,P(k,:),LB,UB,varargin{:});
                end
                ALGOSTEP = 'multiple contraction';
            end
        else
            %% if YTRY better than second-highest point, use this point
            P(end,:) = PTRY;
            Y(end) = YTRY;
            ALGOSTEP = 'reflection';
        end
        
        % the initial temperature is estimated in the first loop from 
        % the number of successfull and unsuccesfull moves, and the average 
        % increase in cost associated with the unsuccessful moves
        
        if TEMP_LOOP_NUMBER == 1 && isempty(OPTIONS.TEMP_START)
            if Y(1) > Y(end)
                SUCCESSFUL_MOVES = SUCCESSFUL_MOVES+1;
            elseif Y(1) <= Y(end)
                UNSUCCESSFUL_MOVES = UNSUCCESSFUL_MOVES+1;
                UNSUCCESSFUL_COSTS = UNSUCCESSFUL_COSTS+(Y(end)-Y(1));
            end
        end

    end

    % stop if previous for loop was broken due to some stop criterion
    if ITERTEMP < MAXITERTEMP
        break;
    end
    
    % store cost function values in COSTS vector
    COSTS((ITERTEMP-1)*NDIM+1:ITERTEMP*NDIM+1) = Y;
    
    % calculated initial temperature or recalculate temperature 
    % using cooling schedule as proposed by Cardoso et al. (1996)
    % -----------------------------------------------------------
    
    if TEMP_LOOP_NUMBER == 1 && isempty(OPTIONS.TEMP_START)
        TEMP = -(UNSUCCESSFUL_COSTS/(SUCCESSFUL_MOVES+UNSUCCESSFUL_MOVES))/log(((SUCCESSFUL_MOVES+UNSUCCESSFUL_MOVES)*OPTIONS.INITIAL_ACCEPTANCE_RATIO-SUCCESSFUL_MOVES)/UNSUCCESSFUL_MOVES);
    elseif TEMP_LOOP_NUMBER ~= 0
        STDEV_Y = std(COSTS);
        COOLING_FACTOR = 1/(1+TEMP*log(1+OPTIONS.COOL_RATE)/(3*STDEV_Y));
        TEMP = TEMP*min(OPTIONS.MIN_COOLING_FACTOR,COOLING_FACTOR);
    end
    
    % add one to temperature loop number
    TEMP_LOOP_NUMBER = TEMP_LOOP_NUMBER+1;
    
end

% return solution
X = PBEST;
FVAL = YBEST;

% store number of function evaluations
OUTPUT.nFUN_EVALS = nFUN_EVALS;

% store number of iterations
OUTPUT.nITERATIONS = nITERATIONS;

% trim OUTPUT data structure
OUTPUT.TEMPERATURE(nITERATIONS+1:end) = [];
OUTPUT.SIMPLEX(:,:,nITERATIONS+1:end) = [];
OUTPUT.SIMPLEX_BEST(nITERATIONS+1:end,:) = [];
OUTPUT.COSTS(nITERATIONS+1:end,:) = [];
OUTPUT.COST_BEST(nITERATIONS+1:end) = [];

% store the amount of time needed in OUTPUT data structure
OUTPUT.TIME = toc;

return

% ==============================================================================

% AMOTRY FUNCTION
% ---------------

function [YFTRY,YTRY,PTRY] = AMOTRY(FUN,P,fac,LB,UB,varargin)
% Extrapolates by a factor fac through the face of the simplex across from 
% the high point, tries it, and replaces the high point if the new point is 
% better.

global NDIM TEMP

% calculate coordinates of new vertex
psum = sum(P(1:NDIM,:))/NDIM;
PTRY = psum*(1-fac)+P(end,:)*fac;

% evaluate the function at the trial point.
YTRY = CALCULATE_COST(FUN,PTRY,LB,UB,varargin{:});
% substract random fluctuations to function values of current vertices
YFTRY = YTRY-TEMP*abs(log(rand(1)));

return

% ==============================================================================

% COST FUNCTION EVALUATION
% ------------------------

function [YTRY] = CALCULATE_COST(FUN,PTRY,LB,UB,varargin)

global YBEST PBEST NDIM nFUN_EVALS

for i = 1:NDIM
    % check lower bounds
    if PTRY(i) < LB(i)
        YTRY = 1e12+(LB(i)-PTRY(i))*1e6;
        return
    end
    % check upper bounds
    if PTRY(i) > UB(i)
        YTRY = 1e12+(PTRY(i)-UB(i))*1e6;
        return
    end
end

% calculate cost associated with PTRY
YTRY = feval(FUN,PTRY,varargin{:});

% add one to number of function evaluations
nFUN_EVALS = nFUN_EVALS + 1;

% save the best point ever
if YTRY < YBEST
    YBEST = YTRY;
    PBEST = PTRY;
end

return

function options = SIMPSASET(varargin)
%SIMPSASETNEW Create/alter simpsa optimization OPTIONS structure.
%   OPTIONS = SIMPSASETNEW('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates an
%   optimization options structure OPTIONS in which the named parameters have
%   the specified values.  Any unspecified parameters are SETNEW to [] (parameters
%   with value [] indicate to use the default value for that parameter when
%   OPTIONS is passed to the optimization function). It is sufficient to type
%   only the leading characters that uniquely identify the parameter.  Case is
%   ignored for parameter names.
%   NOTE: For values that are strings, the complete string is required.
%
%   OPTIONS = SIMPSASETNEW(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of OLDOPTS
%   with the named parameters altered with the specified values.
%
%   OPTIONS = SIMPSASETNEW(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS.  Any parameters in NEWOPTS
%   with non-empty values overwrite the corresponding old parameters in
%   OLDOPTS.
%
%   SIMPSASETNEW with no input arguments and no output arguments displays all
%   parameter names and their possible values, with defaults shown in {}
%   when the default is the same for all functions that use that option -- use
%   SIMPSASETNEW(OPTIMFUNCTION) to see options for a specific function.
%
%   OPTIONS = SIMPSASETNEW (with no input arguments) creates an options structure
%   OPTIONS where all the fields are SETNEW to [].

%   See also SIMPSAGET, SIMPSA
% 
% Copyright (C) 2006 Brecht Donckels, BIOMATH, brecht.donckels@ugent.be
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.



% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf('                   TEMP_START: [ positive scalar ]\n');
    fprintf('                     TEMP_END: [ positive scalar ]\n');
    fprintf('                    COOL_RATE: [ positive scalar ]\n');
    fprintf('     INITIAL_ACCEPTANCE_RATIO: [ positive scalar < 1 {0.95} ]\n');
    fprintf('           MIN_COOLING_FACTOR: [ positive scalar < 1 {0.9}]\n');
    fprintf('          MAX_ITER_TEMP_FIRST: [ positive scalar {100} ]\n');
    fprintf('           MAX_ITER_TEMP_LAST: [ positive scalar {20} ]\n');
    fprintf('               MAX_ITER_TOTAL: [ positive scalar {2500} ]\n');
    fprintf('                     MAX_TIME: [ positive scalar {2500} ]\n');
    fprintf('                MAX_FUN_EVALS: [ positive scalar {2500} ]\n');
    fprintf('                         TOLX: [ positive scalar {1e-6} ]\n');
    fprintf('                       TOLFUN: [ positive scalar {1e-6} ]\n');
    fprintf('                      DISPLAY: [ ''iter'' or ''none'' {''iter''} ]\n');
    fprintf('                   OUTPUT_FCN: [ function_handle ]\n');
    fprintf('\n');
return;
end

Names = [
    'TEMP_START               '
    'TEMP_END                 '
    'COOL_RATE                '
    'INITIAL_ACCEPTANCE_RATIO '
    'MIN_COOLING_FACTOR       '
    'MAX_ITER_TEMP_FIRST      '
    'MAX_ITER_TEMP_LAST       '
    'MAX_ITER_TEMP            '
    'MAX_ITER_TOTAL           '
    'MAX_TIME                 '
    'MAX_FUN_EVALS            '
    'TOLX                     '
    'TOLFUN                   '
    'DISPLAY                  '
    'OUTPUT_FCN               '
    ];

m = size(Names,1);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in odeSETNEW(o1,o2,...).
options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end
i = 1;
while i <= nargin
  arg = varargin{i};
  if ischar(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error('MATLAB:odeSETNEW:NoPropNameOrStruct',...
            ['Expected argument %d to be a string property name ' ...
                     'or an options structure\ncreated with SIMANSETNEW.'], i);
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        options.(deblank(Names(j,:))) = val;
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error('MATLAB:odeSETNEW:ArgNameValueMismatch',...
        'Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~ischar(arg)
      error('MATLAB:odeSETNEW:NoPropName',...
            'Expected argument %d to be a string property name.', i);
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error('MATLAB:odeSETNEW:InvalidPropName',...
            'Unrecognized property name ''%s''.', arg);
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subSETNEWs of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
          msg = [msg ', ' deblank(Names(k,:))];
        end
        msg = sprintf('%s).', msg);
        error('MATLAB:odeSETNEW:AmbiguousPropName', msg);
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error('MATLAB:odeSETNEW:NoValueForProp',...
        'Expected value for property ''%s''.', arg);
end



    