function [varargout] = autogeneratesax(varargin)
%-----------------------------------------------
%%% Automatically generate short-axis view from transversal 3D CT stack.

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%---------------------------------
function autogenerate_init(set_no)
%---------------------------------
% Check if images are dimensions 512 x 512

global SET NO net_p1 net_p2 net_p0

if nargin < 1
  set_no = NO;
end

dims = [size(SET(set_no).IM,1),size(SET(set_no).IM,2)];

if ~isequal(dims,[512,512])
    logdisp('Wrong image dimensions! Must be 512 x 512')
    return
end

logdisp('Loading Networks');
normImgStack = uint8(NormalizeImgStack(SET(set_no).IM)*255);
if isequal(net_p0,[])
    net_p0 = load(['+ct' filesep 'ResNet50_p0.mat']); net_p0 = net_p0.net;
end
if isequal(net_p1,[])
    net_p1 = load(['+ct' filesep 'ResNet50_p1.mat']); net_p1 = net_p1.net;
end
if isequal(net_p2,[])
    net_p2 = load(['+ct' filesep 'ResNet50_p2.mat']); net_p2 = net_p2.net2;
end
logdisp('Predicting Points - Phase 1');

%%% Begin Automatic Short-Axis view generation and LV Segmentation.
%%% Predicts the most suitable slice, and use that for annotation.

% Find suitable slice
delz = predictslice(set_no,normImgStack);
% delz =  round(size(SET(set_no).IM,4)/2);
normImg = normImgStack(:,:,delz);

pred_p1 = double(net_p1.predict(normImg));

% Store predictions in setstruct
% if pred_p1(3) > pred_p1(1) 
SET(set_no).CT.avx = [pred_p1(3),pred_p1(1)];
SET(set_no).CT.avy = [pred_p1(4),pred_p1(2)];
% else
%     SET(set_no).CT.avx = [pred_p1(1),pred_p1(3)];
%     SET(set_no).CT.avy = [pred_p1(2),pred_p1(4)];
% end
SET(set_no).CT.apexx = pred_p1(5);
SET(set_no).CT.apexy = pred_p1(6);
SET(set_no).CT.delz = delz;
SET(set_no).CT.RVpulmix = nan;
SET(set_no).CT.RVpulmiy = nan;
SET(set_no).CT.RVpulmiz = nan;
logdisp('Generating Orthogonal Image');
[orth_img,orthimageorientation,orthimageposition] = GenerateOrthogonalImageWithNewPoints(set_no);
orth_img(isnan(orth_img)) = 0;

% Store scaling factor
x_scale = size(orth_img,1)/640;
y_scale = size(orth_img,2)/426;

% Resize image to fit network 2
orth_img = imresize(orth_img,[640,426]);

% Normalize 
norm_orth_img = uint8(NormalizeImgStack(orth_img)*255);
logdisp('Predicting points - Phase 2');
pred_p2 = net_p2.predict(norm_orth_img);

SET(set_no).CT.orthavx = double([pred_p2(3),pred_p2(1)]*x_scale);
SET(set_no).CT.orthavy = double([pred_p2(4),pred_p2(2)]*y_scale);
SET(set_no).CT.orthapexx = double(pred_p2(5)*x_scale);
SET(set_no).CT.orthapexy = double(pred_p2(6)*y_scale);

SET(set_no).CT.orthimageorientation = orthimageorientation;
SET(set_no).CT.orthimageposition = orthimageposition;

newno = numel(SET)+1;
SET(set_no).CT.generatedsaxno = newno;
logdisp('Generating SAX stack');
ct.imagestacksct('autogeneratesax_Callback',true,set_no)
logdisp('SAX generated')

%----------------------------------------------------
function slice_no = predictslice(set_no,normImgStack)
%----------------------------------------------------
%%% Predicts a suitable slice from a normalized image stack using CNN.

global SET net_p0
% Number of slices
N = size(normImgStack,3);
% Resize and reshape to fit network
images = imresize(normImgStack,[128,128]);
images = reshape(images,[128 128 1 N]);
% Only predict on the middle 50 % of slices
lower_idx = round(N*0.25);
upper_idx = round(N*0.75);
prediction = net_p0.classify(images(:,:,1,lower_idx:upper_idx));
slices = find(prediction == "Good");
if slices
  slice_no = round(mean(slices))+ lower_idx;
else
  disp("Couldn't find suitable slice using network. Will take middle slice!")
  slice_no = round(N/2);
end

%_-------------------------------------------------
function normImgStack = NormalizeImgStack(imgstack)
%--------------------------------------------------
%%% Normalize image stack (or single image) by subtracting mean and
%%% dividing with standard deviation.
%%% Input stack must have size (d1,d2,1,n) where "d1", "d2" are image
%%% dimensions, "n" is number of 2D images.

N = size(imgstack,4);
dim1 = size(imgstack,1);
dim2 = size(imgstack,2);
normImgStack = zeros(dim1,dim2,N);

for i = 1:N  
  tmp = single(imgstack(:,:,i));
  tmp(isnan(tmp)) = 0;
  mu      = nanmean(tmp(:));
  normImg = (tmp-mu)/std(tmp(:));
  normImgStack(:,:,i) = normImg;
end

%--------------------------------------------------------------------------
function [orthimg,orthimageorientation,orthimageposition] = GenerateOrthogonalImageWithNewPoints(set_count)
%--------------------------------------------------------------------------
%%%     Creates orthogonal image. All code is from Segment

global SET NO

tmp     = SET(set_count);
delz    = tmp.CT.delz;
avx     = tmp.CT.avx; avy = tmp.CT.avy;
apexx   = tmp.CT.apexx; apexy = tmp.CT.apexy;
XSize   = tmp.XSize; YSize = tmp.YSize; ZSize = tmp.ZSize;
rawim   = tmp.IM; rawim = reshape(rawim,size(rawim,1),size(rawim,2),size(rawim,4));

avvec = [avx(2) avy(2)]-[avx(1) avy(1)];
apexvec = [apexx apexy]-[avx(1) avy(1)];
apexproj = [avx(1) avy(1)] + dot(avvec,apexvec)/norm(avvec)^2 * avvec;
apexlinex = apexx + (apexproj(1)-apexx)*[-1000 0 1 1000];
apexliney = apexy + (apexproj(2)-apexy)*[-1000 0 1 1000];

apexpoint = [apexlinex(2) apexliney(2)];
apexvec = [apexlinex(3) apexliney(3)]-apexpoint;

%find limits (image borders)
t0 = (0.5-apexpoint)./apexvec;
tmax = ([XSize YSize]+0.5-apexpoint)./apexvec;
tall = sort([t0 tmax]);
tlims = tall(2:3);
xbasevec = [apexvec/norm(apexvec) 0];
xmax = (tlims(2)-tlims(1))*norm(apexvec);
newxsz = xmax; %norm(xbasevec.*scalevec);

ybasevec = [0 0 1];
ymax = ZSize-1;
newysz = ymax; %norm(ybasevec.*scalevec);
zbasevec = -cross(xbasevec,ybasevec);

imagepos = [apexpoint+tlims(1)*apexvec 1];%upper corner beyond apex
%Store orthogonal image position and orientation in coordinates from raw
%image
orthimageposition = imagepos;
orthimageorientation = [xbasevec; ybasevec; zbasevec];

[LINSPX,LINSPY] = ndgrid(linspace(0,xmax,newxsz),linspace(0,ymax,newysz));
X = imagepos(1) + LINSPX * xbasevec(1) + LINSPY * ybasevec(1);
Y = imagepos(2) + LINSPX * xbasevec(2) + LINSPY * ybasevec(2);
Z = imagepos(3) + LINSPX * xbasevec(3) + LINSPY * ybasevec(3);
orthimg = interp3(rawim(:,:,:), Y,X,Z,'nearest');

annopts = [avx apexx; avy apexy; delz*[1 1 1]];
vecs = orthimageorientation*(annopts-repmat(imagepos',1,3));

