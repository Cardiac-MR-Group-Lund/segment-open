function cxt1map2segment()
global SET NO DATA
%OBS: the CX T1 map variable in the selected .mat file needs to be 
%the same size (and from the same subject) as theselected image stack.  

%Instructions for use: 
% Added functionality for practical usage: 
% 
% Importing custom T1 maps for a subject: 
% 
% instructions: 
% 1. Load in some MOLLI images (either just the map, one TI, all TI images or all TI images + T1 map)
% 
% 2. Draw your ROIs
% 
% 3. Push 'ctrl-alt-space' and select the .mat file containing your custom T1 map. 
% 
% 4. Done. The custom T1 map is displayed as a new image stack.  

%Load CX T1 map. 
[filename,pathname] = myuigetfile('*.mat','Please select one of christos T1-maps');
if (filename(1)~=0)
load([pathname, filename]);
variableinfo = whos('-file',[pathname, filename]);
%Search for mathcing image in file: It will take the first one it finds.
nbrofvariables = length(variableinfo);

if nbrofvariables>0
sizematches = ones(nbrofvariables,1);

%Mark current timeframe: 
ind = ((1:SET(NO).TSize)==SET(NO).CurrentTimeFrame);

%Duplicate imagestack: 
tools('duplicateimagestack_Callback')

%Remove all timeframes except one. 
DATA.Silent = true;
tools('removetimeframes',ind);
DATA.Silent = false;
%Extract image size: 
setimsize = size(SET(NO).IM);
for n = 1:nbrofvariables
	%Value should be zero if the images matches: 
	sizematches(n) = nnz(bsxfun(@minus,variableinfo(n).size,setimsize));
end
% T1map = T1_map_mean_NORM;
k = find(sizematches==0,1,'first');
if ~isempty(k)
	%assign T1map: 
	eval(['T1map = ', variableinfo(k).name]);
	T1map = single(T1map);
	
	%This assumes that the T1 map has the same variable name as the
	%filename: 
	SET(NO).IntensityOffset = min(T1map(:));
	T1maprescaled = T1map-SET(NO).IntensityOffset;
	SET(NO).IntensityScaling = max(T1maprescaled(:));
	
	%Scale T1 map: 
	SET(NO).IM = T1maprescaled/SET(NO).IntensityScaling;
	SET(NO).VENC = [];
else
	myfailed('Error: T1map and segment image dimensions do not agree. Nothing was changed. ');
end

calcpreview = true;
segment('makeviewim',DATA.CurrentPanel,NO);
segment('updatemodeldisplay');
drawfunctions('drawimageno');
drawfunctions('drawthumbnails',calcpreview);
end
end