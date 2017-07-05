function smoothedImage=smooth3D(image,smoothradius, ResolutionX, ResolutionY, ResolutionZ)

XSize=size(image,1);
YSize=size(image,2);
TSize=size(image,3);
ZSize=size(image,4);
timeframes=1:TSize;

%smoothing filter
n = 9;

%X direction
x = linspace(-n*ResolutionX,n*ResolutionX,2*n+1);
f = exp(-(x.^2)/(smoothradius.^2));
f = f./sum(f(:));
f = single(f);
fx = f(:);

%Y direction
x = linspace(-n*ResolutionY,n*ResolutionY,2*n+1);
f = exp(-(x.^2)/(smoothradius.^2));
f = f./sum(f(:));
f = single(f);
fy = f(:)';

%Z direction
x = linspace(-n*ResolutionZ,n*ResolutionZ,2*n+1);
f = exp(-(x.^2)/(smoothradius.^2));
f = f./sum(f(:));
f = single(f);
fz = reshape(f,[1 1 length(f)]);

%smooth
temp=repmat(single(0),[XSize YSize TSize ZSize]);
if ZSize>1
  %Smooth in x-dir
  for tloop=timeframes
    temp(:,:,tloop,:) = econv3(...
      single(squeeze(image(:,:,tloop,:))),fx);
  end;

  %Smooth in y-dir
  for tloop=timeframes
    temp(:,:,tloop,:) = econv3(squeeze(temp(:,:,tloop,:)),fy);
  end;

  %Smooth in z-dir
  for tloop=timeframes
    temp(:,:,tloop,:) = econv3(squeeze(temp(:,:,tloop,:)),fz);
  end;
else
  %Smooth in x-dir
  temp = econv3(...
    single(squeeze(image)),fx);

  %Smooth in y-dir
  temp = econv3(squeeze(temp),fy);
end
smoothedImage=temp;
