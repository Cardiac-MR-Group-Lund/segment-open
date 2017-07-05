function imissingle=classcheckim(nos)
%This function check if image is single and if not asks if to convert.

%Jane Sjögren, Erik Södervall.

global NO DATA

imissingle=true;

if nargin==0
  nos=NO;
end

int16nos=[];
for j=nos
  if checkifint16(j)
    int16nos=[int16nos j]; %#ok<AGROW>
  end
end

if not(isempty(int16nos))
  s=dprintf('This function is not supported for images loaded as integers. To proceed you must convert image number ');
  for j=int16nos
    s=sprintf('%s %d, ',s,j);
  end
  s=s(1:end-1);
  dprintf('%s into float. Do you want to do that now?',s);
  convert=yesno(s,'',DATA.GUI.Segment);
  if convert
    for j=int16nos
      convertimtosingle(j);
    end
  else
    imissingle=false;
  end
end

%------------------------------
function check=checkifint16(no)
%------------------------------
global SET NO

if nargin==0
	no=NO;
end
check=isa(SET(no).IM,'int16');

%-----------------------------
function convertimtosingle(no)
%-----------------------------
global SET NO

if nargin==0
  no=NO;
end
im=single(SET(no).IM);
maxValue=max(im(:));
minValue=min(im(:));
im=(im-minValue)/(maxValue-minValue);
SET(no).IntensityScaling=maxValue-minValue;
SET(no).IntensityOffset=minValue;
SET(no).IM=im;