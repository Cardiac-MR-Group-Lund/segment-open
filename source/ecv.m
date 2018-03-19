function [varargout] = ecv(varargin)
%-----------------------------------
% Calculate ECV from MOLLI T1 stacks pre and post Gd

if nargin == 0
  varargin{1} = 'init_Callback';
end

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard


%---------------------
function init_Callback %#ok<DEFNU>
%---------------------
%initialization of the ECV analysis

global SET NO

%find the thumbnail number of preT1map and postT1map
prenos = strcmp(strtrim({SET.ImageType}),'T1 map Pre');
postnos = strcmp(strtrim({SET.ImageType}),'T1 map Post');
prenos = find(prenos);
postnos = find(postnos);
  
%find image stacks with ROI blood
preno = [];
for no = prenos
  islvendo = existfunctions('existendoinslices',no);
  islvepi = existfunctions('existepiinslices',no);
  if SET(no).RoiN >= 1 && any(strcmp({SET(no).Roi.Name},'Blood')) && sum(islvendo+islvepi);
    preno = [preno no];
  end
end
postno = [];
for no = postnos
  islvendo = existfunctions('existendoinslices',no);
  islvepi = existfunctions('existepiinslices',no);
  if SET(no).RoiN >= 1 && any(strcmp({SET(no).Roi.Name},'Blood')) && sum(islvendo+islvepi);
    postno = [postno no];
  end
end

multislice = false;
singslicenos = [SET.ZSize] == 1;
if (numel(preno) > 1 && any(not(singslicenos(preno)))) || ...
    (numel(postno) > 1 && any(not(singslicenos(postno)))) 
  if ismember(NO,[preno postno])
    if ismember(NO,preno)
      preno = NO;
      if numel(postno) > 1
        for k =1:length(postno),postnostri{k}=postno(k);end
        postnonbr = mymenu('Select post T1 map NO',postnostri);
        if isempty(postnonbr)
          mywarning('Not correctly defined image stack number');
          return;
        else
          postno = postno(postnonbr);
        end
      end
    else
      postno = NO;
      if numel(preno) > 1
        for k =1:length(preno),prenostri{k}=preno(k);end
        prenonbr = mymenu('Select post T1 map NO',prenostri);
        if isempty(prenonbr)
          mywarning('Not correctly defined image stack number');
          return;
        else
          preno = preno(prenonbr);
        end
      end
    end
  else
    myfailed('Multi-slice data. Select the image stack to analyze.');
    return;
  end
end
%check if multi-slice data, if then only use one image stack pair
if any(not(singslicenos(preno))) || any(not(singslicenos(postno)))
  multislice = true;
  if length(preno) > 1 || length(postno) > 1
    %if multi-slice data, use only NO if that is applicable
    if ismember(NO,[preno postno]);
      if ismember(NO,preno)
        prenonbr = find(ismember(NO,preno));
        suggestedpostno = postno(min(prenonbr,length(postno)));        
        if length(postno) > 1
          for k =1:length(postno),postnostri{k}=postno(k);end
          postnonbr = mymenu('Select post T1 map NO',postnostri);
          if isempty(postnonbr)
            mywarning('Not correctly defined image stack number');
            return;
          else
            postno = postno(postnonbr);
          end
        else
          postno = suggestedpostno;
        end  
        preno = NO;
      elseif ismember(NO,postno)
        postnonbr = find(ismember(NO,postno));
        suggestedpreno = preno(min(postnonbr,length(preno)));
        if length(preno) > 1
          for k =1:length(preno),prenostri{k}=preno(k);end
          prenonbr = mymenu('Select pre T1 map NO',prenostri);
          if isempty(prenonbr)
            mywarning('Not correctly defined image stack number');
            return;
          else
            preno = preno(prenonbr);
          end
        else
          preno = suggestedpreno;
        end
        postno = NO;
      end
    else
      myfailed('Multi-slice data. Select the image stack to analyze.');
      return;
    end
  end
end

%error checks
if isempty(preno) || isempty(postno)
  myfailed(['Could not find T1 map pre or T1 map post image stacks with' ...
    ' ROI''s labelled ''Blood''.']);
  return
end
if length(preno) ~= length(postno) %SET(preno).ZSize ~= SET(postno).ZSize
  mywarning('Number of image stacks in pre- and post-T1 map must be equal');
  return;
end
% if length(preno) > 3
%   mywarning('Only supports up to 3 T1 map image stack pairs.');
%   return;
% end

%ask for the hematocrit value
while true
  hem = inputdlg('Hematocrit [0 1]','Input value');
  if isempty(hem)
    return
  else
    hem = str2double(hem{1});
    if hem < 0 && hem > 1
      mywarning('The hematocrit value must be within [0 1]');
    end
    if hem >= 0 && hem <= 1
      break
    end
  end
end

%find the pre T1 map image, ask if there are more than one image (tf>1)
askno = [];
loop = 1;
for sliceloop = 1:length(preno)
  if SET(preno(sliceloop)).TSize > 1
    askno = [askno preno(sliceloop)];
    if loop==1
      askstri{loop} = dprintf('Time frame for\n\npre T1 map NO %s',num2str(preno(sliceloop)));
    else
      askstri{loop} = dprintf('pre T1 map NO %s',num2str(preno(sliceloop)));
    end
    loop = loop+1;
  end
end
if ~isempty(askno)
  %initial guess on NO based on existens of LV segmentation
  isLV=[];
  for j=1:length(preno)
    if ~isempty(SET(preno(j)).EndoX)
      isLV(:,j)=sum(~isnan(SET(preno(j)).EndoX(1,:,:)),3)';
    end
  end
  for j=1:length(preno)
    if ~isempty(isLV) && sum(isLV(:,j))>1
      [~,isLVind] = max(isLV(:,j));
      guess{j} = num2str(isLVind);
    else
      guess{j}='';
    end
  end
  %find time frame
  while true
    premaptftemp = inputdlg(askstri,'Pre T1 map',1,guess);
    if isempty(premaptftemp)
      return;
    end
    ok = [];
    premaptf = zeros(1,length(preno));
    for k = 1:length(premaptftemp)
      premaptf(k) = str2double(premaptftemp{k});
      if isnan(premaptf(k)) || isempty(premaptf(k)) || premaptf(k)<1 || premaptf(k)>SET(askno(k)).TSize
        ok = [ok 0];
      else
        ok = [ok 1];
      end
    end
    if sum(ok) == length(premaptf)
      break
    else
      mywarning('Not correctly defined time frame number');
    end
  end
else
  premaptf = ones(1,length(preno));
end

%find the post T1 map image, ask if there are more than one image (tf>1)
askno = [];
askstri = {};
loop = 1;
for sliceloop = 1:length(postno)
  if SET(postno(sliceloop)).TSize > 1
    askno = [askno postno(sliceloop)];
    if sliceloop == 1
      askstri{loop} = dprintf('Time frame for \n\npost T1 map NO %s',num2str(postno(sliceloop)));
    else
      askstri{loop} = dprintf('post T1 map NO %s',num2str(postno(sliceloop)));
    end
    loop = loop+1;
  end
end
if ~isempty(askno)
  %initial guess on NO based on existens of LV segmentation
  isLV=[];
  for j=1:length(postno)
    if ~isempty(SET(postno(j)).EndoX)
      isLV(:,j)=sum(~isnan(SET(postno(j)).EndoX(1,:,:)),3)';
    end
  end
  for j=1:length(postno)
    if ~isempty(isLV)
      [~,isLVind] = max(isLV(:,j));
      guess{j} = num2str(isLVind);
    else
      guess{j}='';
    end
  end
  while true
    postmaptftemp = inputdlg(askstri,'Post T1 maps',1,guess);
    if isempty(postmaptftemp)
      return;
    end
    ok = [];
    postmaptf = zeros(1,length(postno));
    for k = 1:length(postmaptftemp)
      postmaptf(k) = str2double(postmaptftemp{k});
      if isnan(postmaptf(k)) || isempty(postmaptf(k)) || postmaptf(k)<1 || postmaptf(k)>SET(askno(k)).TSize
        ok = [ok 0];
      else
        ok = [ok 1];
      end
    end
    if sum(ok) == length(postmaptf)
      break
    else
      mywarning('Not correctly defined time frame number');
    end
  end
else
  postmaptf = ones(1,length(postno));
end

%If more than one time frame, ask so the pre- and post images are correctly
%coupled
if length(preno) > 1
  for k=1:length(postno)
    guess{k} = num2str(postno(k));
    if k==1
      askstri{1} = dprintf('Coupled post T1 map for\n\npre T1 map NO %s',num2str(preno(k)));
    else
      askstri{k} = dprintf('pre T1 map NO %s',num2str(preno(k)));
    end
  end
  while true
    coupledpostnotemp = inputdlg(askstri,'Pre T1 maps',1,guess);
    if isempty(coupledpostnotemp)
      return;
    end
    ok = [];
    coupledpostno = zeros(1,length(preno));
    for k = 1:length(coupledpostnotemp)
      coupledpostno(k) = str2double(coupledpostnotemp{k});
      if isnan(coupledpostno(k)) || isempty(coupledpostno(k)) || ~ismember(coupledpostno(k),postno)
        ok = [ok 0];
      else
        ok = [ok 1];
      end
    end
    if sum(ok) == length(coupledpostno)
      postno = coupledpostno; %sort the postno according to input
      break
    else
      mywarning('Not correctly defined post T1 map NO');
    end
  end
end

%check so the images have the same resolution
notok = 0;
for k = 1:length(preno)
  if abs(SET(preno(k)).ResolutionX - SET(postno(k)).ResolutionX) > 1e-2 || ...
      abs(SET(preno(k)).ResolutionY - SET(postno(k)).ResolutionY) > 1e-2
    notok = 1;
  end
end
if notok
  mywarning('The pre T1 map and post T1 map must have the same image resolution');
  return;
end

if multislice
  %find slices with LV segmentation
  endosegslicepre = find(existfunctions('existendoinslices',preno,premaptf));
  endosegslicepost = find(existfunctions('existendoinslices',postno,postmaptf));
  episegslicepre = find(existfunctions('existepiinslices',preno,premaptf));
  episegslicepost = find(existfunctions('existepiinslices',postno,postmaptf));
  lvsegslicepre = intersect(endosegslicepre,episegslicepre);
  lvsegslicepost = intersect(endosegslicepost,episegslicepost);
  %need to be the LV segmentation in the same number of slices in pre and
  %post
  if length(lvsegslicepre) ~= length(lvsegslicepost)
    myfailed('Need to have LV segmentation in the same number of slices in the pre and post T1 maps');
    return;
  end
else
  %check so there is LV segmentation in both images
  for k = 1:length(preno)
    if isempty(SET(preno(k)).EndoX) || isempty(SET(preno(k)).EpiX) || ...
        all(isnan(SET(preno(k)).EndoX(1,premaptf(k),:))) || all(isnan(SET(preno(k)).EpiX(1,premaptf(k),:)))
      notok = 1;
    end
  end
  if notok
    mywarning('No LV segmentation in pre T1 map found');
    return
  end
  lvsegslicepre = 1;
  for k = 1:length(preno)
    if isempty(SET(postno(k)).EndoX) || isempty(SET(postno(k)).EpiX) || ...
        all(isnan(SET(postno(k)).EndoX(1,postmaptf(k),:))) || all(isnan(SET(postno(k)).EpiX(1,postmaptf(k),:)))
      notok = 1;
    end
  end
  if notok
    mywarning('No LV segmentation in post T1 map found');
    return
  end
  lvsegslicepost = 1;
end

%open the GUI for visualization
opengui(preno,postno,hem,premaptf,postmaptf,multislice,lvsegslicepre,lvsegslicepost);

%------------------------
function roiinit_Callback %#ok<DEFNU>
%------------------------
% Calculate ECV from MOLLI T1 stacks pre and post Gd within predefined ROIs

global SET 

rois = {};

prenos = strcmp(strtrim({SET.ImageType}),'T1 map Pre');
postnos = strcmp(strtrim({SET.ImageType}),'T1 map Post');
% prenos = strcmp(strtrim({SET.SeriesDescription}),'MOLLI pre-Gd');
% postnos = strcmp(strtrim({SET.SeriesDescription}),'MOLLI post-Gd');
singframenos = [SET.TSize] == 1;

prenos = find(prenos & singframenos);
postnos = find(postnos & singframenos);

preno = [];
for no = prenos
  if SET(no).RoiN > 1 && any(strcmp({SET(no).Roi.Name},'Blood'))
    preno = [preno no];
  end
end
postno = [];
for no = postnos
  if SET(no).RoiN > 1 && any(strcmp({SET(no).Roi.Name},'Blood'))
    postno = [postno no];
  end
end
if isempty(preno) || isempty(postno)
  myfailed(['Could not find single frame pre- and post-Gd image stacks' ...
    ' containing multiple ROI''s including one labelled ''Blood''.']);
  return
elseif numel(preno) > 1 || numel(postno) > 1
  mywarning(['Found more than one stack of properly prepared ' ...
    'pre- or post-Gd images. Taking first (arbitrary decision).']);
  preno = preno(1);
  postno = postno(1);
end

while true
  hem = inputdlg('Hematocrit','Input value');
  if isempty(hem)
    return
  else
    hem = str2double(hem{1});
    if hem >= 0 && hem <= 1
      break
    end
  end
end

premyo = [];
postmyo = [];
for preroi = SET(preno).Roi
  roiindinpost = find(strcmp(preroi.Name,{SET(postno).Roi.Name}));
  if numel(roiindinpost) > 1
    mywarning(['Found multiple ROI''s with same label in post-Gd ' ...
      'image stack. Taking first (arbitrary decision).']);
    roiindinpost = roiindinpost(1);
  end
  if ~isempty(roiindinpost)
    postroi = SET(postno).Roi(roiindinpost);
    premask = segment('createmask',...
      [SET(preno).XSize SET(preno).YSize],...
      preroi.Y,preroi.X);
    postmask = segment('createmask',...
      [SET(postno).XSize SET(postno).YSize],...
      postroi.Y,postroi.X);
    preim = calcfunctions('calctruedata',SET(preno).IM(:,:,1,preroi.Z),preno);
    postim = calcfunctions('calctruedata',SET(postno).IM(:,:,1,postroi.Z),postno);
    if strcmp(preroi.Name,'Blood')
      preblood = mean(preim(premask));
      postblood = mean(postim(postmask));
    else
      rois = [rois preroi.Name];
      premyo = [premyo mean(preim(premask))];
      postmyo = [postmyo mean(postim(postmask))];
    end
  end
end

extracellularvolume = (1-hem).*(1./postmyo-1./premyo)/(1/postblood-1/preblood);

if nargout == 0
  stri = '';
  for i = 1:numel(rois)
    stri = [stri sprintf('%s: %0.4f\n',rois{i},extracellularvolume(i))];
  end
  mymsgbox(stri,'ECV Report');
  segment('cell2clipboard',[rois' num2cell(extracellularvolume')]);
end

%----------------------------------------------------
function opengui(preno,postno,hem,premaptf,postmaptf,multislice,lvsegslicepre,lvsegslicepost)
%----------------------------------------------------
%Open the ECV registration GUI. Performs automatic rigid registration of
%pre and post Gd images

%written by Helen Fransson 2014-06-12, updated 2015-06-23

global DATA SET

for k = 1:length(preno)
  %create image and myocardial mask
  preT1maptemp = SET(preno(k)).IM(:,:,premaptf(k),:);
  postT1maptemp = SET(postno(k)).IM(:,:,postmaptf(k),:);
%   preendomask = segment('createmask',...
%     [SET(preno(k)).XSize SET(preno(k)).YSize],...
%     SET(preno(k)).EndoY(:,premaptf(k)),SET(preno(k)).EndoX(:,premaptf(k)));
%   preepimask = segment('createmask',...
%     [SET(preno(k)).XSize SET(preno(k)).YSize],...
%     SET(preno(k)).EpiY(:,premaptf(k)),SET(preno(k)).EpiX(:,premaptf(k)));
%   premyomasktemp = preepimask-preendomask;
  
  %calculate blood intensity
  bloodroispost = find(strcmp({SET(postno(k)).Roi.Name},'Blood'));
  bloodroispre = find(strcmp({SET(preno(k)).Roi.Name},'Blood'));
  %post blood
  if multislice
    if numel(bloodroispost) > 1
      for roiloop = 1:length(bloodroispost)
        postroi = SET(postno(k)).Roi(bloodroispost(roiloop));
        roiz(roiloop) = postroi.Z;
      end
      for sliceloop = 1:length(lvsegslicepre)
        slicerois = find(roiz == lvsegslicepre(sliceloop));
        roiim = [];
        if isempty(slicerois)
          %no roi in this slice, take mean in rois from all other sices          
          for roiloop = 1:length(bloodroispost)
            postroi = SET(postno(k)).Roi(bloodroispost(roiloop));
            postmask = segment('createmask',...
              [SET(postno(k)).XSize SET(postno(k)).YSize],...
              postroi.Y(:,postmaptf(k)),postroi.X(:,postmaptf(k)));
            postim = calcfunctions('calctruedata',postT1maptemp(:,:,1,postroi.Z),postno(k));
            roiim = [roiim;postim(postmask)];
          end          
        else
          %take mean of the rois in this slice
          for roiloop = 1:length(slicerois)
            postroi = SET(postno(k)).Roi(bloodroispost(slicerois(roiloop)));
            postmask = segment('createmask',...
              [SET(postno(k)).XSize SET(postno(k)).YSize],...
              postroi.Y(:,postmaptf(k)),postroi.X(:,postmaptf(k)));
            postim = calcfunctions('calctruedata',postT1maptemp(:,:,1,postroi.Z),postno(k));
            roiim = [roiim;postim(postmask)];
          end
        end
        postblood(lvsegslicepre(sliceloop)) = mean(roiim);
      end
    elseif numel(bloodroispost) == 1
      postroi = SET(postno(k)).Roi(bloodroispost);
      postmask = segment('createmask',...
        [SET(postno(k)).XSize SET(postno(k)).YSize],...
        postroi.Y(:,postmaptf(k)),postroi.X(:,postmaptf(k)));
      postim = calcfunctions('calctruedata',postT1maptemp(:,:,1,postroi.Z),postno(k));
      postblood = repmat(mean(postim(postmask)),[1 length(lvsegslicepre)]);
    else
      myfailed('No ROI labeled ''Blood'' found in post T1 map');
      return;
    end
    %pre blood
    if numel(bloodroispre) > 1
      for roiloop = 1:length(bloodroispre)
        preroi = SET(preno(k)).Roi(bloodroispre(roiloop));
        roiz(roiloop) = preroi.Z;
      end
      for sliceloop = 1:length(lvsegslicepre)
        slicerois = find(roiz == lvsegslicepre(sliceloop));
        roiim = [];
        if isempty(slicerois)
          %no roi in this slice, take mean in rois from all other sices          
          for roiloop = 1:length(bloodroispre)
            preroi = SET(preno(k)).Roi(bloodroispre(roiloop));
            premask = segment('createmask',...
              [SET(preno(k)).XSize SET(preno(k)).YSize],...
              preroi.Y(:,premaptf(k)),preroi.X(:,premaptf(k)));
            preim = calcfunctions('calctruedata',preT1maptemp(:,:,1,preroi.Z),preno(k));
            roiim = [roiim;preim(premask)];
          end          
        else
          %take mean of the rois in this slice
          for roiloop = 1:length(slicerois)
            preroi = SET(preno(k)).Roi(bloodroispre(slicerois(roiloop)));
            premask = segment('createmask',...
              [SET(preno(k)).XSize SET(preno(k)).YSize],...
              preroi.Y(:,premaptf(k)),preroi.X(:,premaptf(k)));
            preim = calcfunctions('calctruedata',preT1maptemp(:,:,1,preroi.Z),preno(k));
            roiim = [roiim;preim(premask)];
          end
        end
        preblood(lvsegslicepre(sliceloop)) = mean(roiim);
      end
    elseif numel(bloodroispre) == 1
      preroi = SET(preno(k)).Roi(bloodroispre);
      premask = segment('createmask',...
        [SET(preno(k)).XSize SET(preno(k)).YSize],...
        preroi.Y(:,premaptf(k)),preroi.X(:,premaptf(k)));
      preim = calcfunctions('calctruedata',preT1maptemp(:,:,1,preroi.Z),preno(k));
      preblood = repmat(mean(preim(premask)),[1 length(lvsegslicepre)]);
    else
      myfailed('No ROI labeled ''Blood'' found in pre T1 map');
      return;
    end
      
  else
    if numel(bloodroispost) > 1
      disp('Found multiple ROI''s labeled ''Blood''. Taking first (arbitrary decision).');
      bloodroipost = bloodroispost(1);
    else
      bloodroipost = bloodroispost;
    end
    if numel(bloodroispre) > 1
      disp('Found multiple ROI''s labeled ''Blood''. Taking first (arbitrary decision).');
      bloodroipre = bloodroispre(1);
    else
      bloodroipre = bloodroispre;
    end
    postroi = SET(postno(k)).Roi(bloodroipost);
    postmask = segment('createmask',...
      [SET(postno(k)).XSize SET(postno(k)).YSize],...
      postroi.Y(:,postmaptf(k)),postroi.X(:,postmaptf(k)));
    preroi = SET(preno(k)).Roi(bloodroipre);
    premask = segment('createmask',...
      [SET(preno(k)).XSize SET(preno(k)).YSize],...
      preroi.Y(:,premaptf(k)),preroi.X(:,premaptf(k)));
    preim = calcfunctions('calctruedata',preT1maptemp(:,:,1,preroi.Z),preno(k));
    postim = calcfunctions('calctruedata',postT1maptemp(:,:,1,postroi.Z),postno(k));
    preblood(k) = mean(preim(premask));
    postblood(k) = mean(postim(postmask));
  end
      
%   for preroi = SET(preno(k)).Roi
%     roiindinpost = find(strcmp(preroi.Name,{SET(postno(k)).Roi.Name}));
%     if numel(roiindinpost) > 1
%       mywarning(['Found multiple ROI''s with same label in post-Gd ' ...
%         'image stack. Taking first (arbitrary decision).']);
%       roiindinpost = roiindinpost(1);
%     end
%     if ~isempty(roiindinpost) && strcmp(preroi.Name,'Blood')
%       postroi = SET(postno(k)).Roi(roiindinpost);
%       premask = segment('createmask',...
%         [SET(preno(k)).XSize SET(preno(k)).YSize],...
%         preroi.Y(:,premaptf(k)),preroi.X(:,premaptf(k)));
%       postmask = segment('createmask',...
%         [SET(postno(k)).XSize SET(postno(k)).YSize],...
%         postroi.Y(:,postmaptf(k)),postroi.X(:,postmaptf(k)));
%       preim = calcfunctions('calctruedata',preT1maptemp(:,:,1,preroi.Z),preno(k));
%       postim = calcfunctions('calctruedata',postT1maptemp(:,:,1,postroi.Z),postno(k));
%       preblood = mean(preim(premask));
%       postblood = mean(postim(postmask));
%     end
%   end
  
  %calculate crop size (equal for all images)
  [cropcenterpretemp,cropcenterposttemp,cropradiustemp] = calccropsize(preno(k),postno(k),premaptf(k),postmaptf(k));
  if multislice
    cropcenterpre = squeeze(cropcenterpretemp);
    cropcenterpost = squeeze(cropcenterposttemp);
    cropradius = squeeze(cropradiustemp)';
  else
    cropcenterpre(:,k) = cropcenterpretemp;
    cropcenterpost(:,k) = cropcenterposttemp;
    cropradius(:,k) = cropradiustemp;
  end
  imagesizepre(:,k) = [SET(preno(k)).XSize SET(preno(k)).YSize];
  imagesizepost(:,k) = [SET(postno(k)).XSize SET(postno(k)).YSize];
  
end
cropradiusshared = calcsharedcropsize(cropcenterpre,cropcenterpost,cropradius,imagesizepre,imagesizepost);

for k = 1:length(preno)  
  %crop images and LV segmentation and ROIs
  if multislice
    for slice=1:length(lvsegslicepre)
      [preT1mapcroptemp(:,:,1,lvsegslicepre(slice)),postT1mapcroptemp(:,:,1,lvsegslicepre(slice)), ...
        endoxtemp(:,1,lvsegslicepre(slice)),endoytemp(:,1,lvsegslicepre(slice)), ...
        epixtemp(:,1,lvsegslicepre(slice)),epiytemp(:,1,lvsegslicepre(slice)), ...
        roitemp,roimasktemp{slice},roiareatemp{slice}] = ...
        cropmaps(preno(k),postno(k),premaptf(k),postmaptf(k), ...
        cropcenterpre(:,lvsegslicepre(slice)),cropcenterpost(:,lvsegslicepre(slice)),cropradiusshared,lvsegslicepre(slice));
      roistemp{slice} = roitemp;
    end
    preT1mapcrop = preT1mapcroptemp;
    postT1mapcrop = postT1mapcroptemp;
    endox{k} = endoxtemp;
    endoy{k} = endoytemp;
    epix{k} = epixtemp;
    epiy{k} = epiytemp;
    rois{k} = roistemp;
    roimask{k} = roimasktemp;
    roiarea{k} = roiareatemp;
  else
%   preT1mapcrop = []; postT1mapcrop = [];
    [preT1mapcrop,postT1mapcrop,endox{k},endoy{k},epix{k},epiy{k},rois{k},roimask{k},roiarea{k}] = ...
      cropmaps(preno(k),postno(k),premaptf(k),postmaptf(k),cropcenterpre(:,k),cropcenterpost(:,k),cropradiusshared,1);
  end
  
  %calculates true intensity data
  preT1map{k} = calcfunctions('calctruedata',preT1mapcrop,preno(k));
  postT1map{k} = calcfunctions('calctruedata',postT1mapcrop,postno(k));
  
%   %align images
%   setstr = SET(preno(1));
%   setstr.ImageType = [SET(preno(1)).ImageType ' Aligned'];
%   setstr.IM = [];
%   setstr.IM(:,:,1) = preT1map{k};
%   setstr.IM(:,:,2) = postT1map{k};
%   setstr.TSize = 2;
%   setnew = alignSlides(setstr,1,2);
%   postT1mapaligned{k} = setnew.IM(:,:,2);
  postT1mapaligned{k} = postT1map{k};
  
  %calculate ECV map
  if multislice %&& length(postblood)>1
    postT1mapslice = postT1map{k};
    preT1mapslice = preT1map{k};
    postT1mapalignedslice = postT1mapaligned{k};
    for slice = 1:length(lvsegslicepre)
      ecvmaptemp(:,:,1,lvsegslicepre(slice)) = (1-hem).*(1./postT1mapslice(:,:,1,lvsegslicepre(slice))-1./preT1mapslice(:,:,1,lvsegslicepre(slice)))/ ...
        (1/postblood(lvsegslicepre(slice))-1/preblood(lvsegslicepre(slice)));
      ecvmapalignedtemp(:,:,1,lvsegslicepre(slice)) = (1-hem).*(1./postT1mapalignedslice(:,:,1,lvsegslicepre(slice))-1./preT1mapslice(:,:,1,lvsegslicepre(slice)))/ ...
        (1/postblood(lvsegslicepre(slice))-1/preblood(lvsegslicepre(slice)));
    end
    ecvmap{k} = ecvmaptemp;
    ecvmapaligned{k} = ecvmapalignedtemp;
  else
    ecvmap{k} = (1-hem).*(1./postT1map{k}-1./preT1map{k})/(1/postblood(k)-1/preblood(k));
    ecvmapaligned{k} = (1-hem).*(1./postT1mapaligned{k}-1./preT1map{k})/(1/postblood(k)-1/preblood(k));
  end
%   ecvmapmasked = ecvmap{k}.*premyomask{k};
  %calculate ECV values for each ROI
  preimroiall = preT1map{k};
  postimroiall = postT1map{k};
  postimroialignedall = postT1mapaligned{k};
  roimaskz = roimask{k};
  roiecvtemp = [];
  roiecvalignedtemp = [];
  
  if multislice
    allrois = rois{k};
    for sliceloop = 1:length(allrois)
      roislices = allrois{sliceloop};
      for roiloop = 1:length(roislices)
        roi = roislices(roiloop);
        slice = roi.Z;
        preimroi = preimroiall(:,:,1,slice);
        postimroi = postimroiall(:,:,1,slice);
        postimroialigned = postimroialignedall(:,:,1,slice);
        roimaskzslice = roimaskz{sliceloop};
        preroimyo = preimroi(roimaskzslice{roiloop});
        postroimyo = postimroi(roimaskzslice{roiloop});
        postroimyoaligned = postimroialigned(roimaskzslice{roiloop});
        temproiecvtemp = (1-hem).*(1./postroimyo-1./preroimyo)/(1/postblood(slice)-1/preblood(slice));
        roiecvtemp{sliceloop}{roiloop} = max(0,min(1,temproiecvtemp(not(isinf(temproiecvtemp)))));  %exclude Inf values, ECV within 0-100
        temproiecvalignedtemp = (1-hem).*(1./postroimyoaligned-1./preroimyo)/(1/postblood(slice)-1/preblood(slice));
        roiecvalignedtemp{sliceloop}{roiloop} = max(0,min(1,temproiecvalignedtemp(not(isinf(temproiecvalignedtemp))))); %exclude Inf values, ECV within 0-100
      end
    end
  else
    j = 1;
    for roi = rois{k}
      preimroi = preimroiall;
      postimroi = postimroiall;
      postimroialigned = postimroialignedall;
      preroimyo = preimroi(roimaskz{j});
      postroimyo = postimroi(roimaskz{j});
      postroimyoaligned = postimroialigned(roimaskz{j});
      temproiecvtemp = (1-hem).*(1./postroimyo-1./preroimyo)/(1/postblood(k)-1/preblood(k));
      roiecvtemp{j} = max(0,min(1,temproiecvtemp(not(isinf(temproiecvtemp))))); %ECV within 0-100
      temproiecvalignedtemp = (1-hem).*(1./postroimyoaligned-1./preroimyo)/(1/postblood(k)-1/preblood(k));
      roiecvalignedtemp{j} = max(0,min(1,temproiecvalignedtemp(not(isinf(temproiecvalignedtemp))))); %ECV within 0-100
      j = j+1;
    end
  end
  roiecv{k} = roiecvtemp; 
  roiecvaligned{k} = roiecvalignedtemp;
end

%%%%%%%%%%%%%%%% Open GUI %%%%%%%%%%%%%%%%%
if isopengui('ecvregistration.fig');
  gui = DATA.GUI.ECVRegistration;
  figure(gui.fig);
else
  DATA.GUI.ECVRegistration = mygui('ecvregistration.fig');
  gui = DATA.GUI.ECVRegistration;
  myadjust(gui.fig,DATA.GUI.Segment);
end
gui.preno = preno;
gui.postno = postno;
gui.preT1map = preT1map;
gui.postT1map = postT1map;
gui.postT1mapaligned = postT1mapaligned;
gui.postT1mapplot = gui.postT1mapaligned;
gui.ecvmap = ecvmap;
gui.ecvmapaligned = ecvmapaligned;
gui.ecvmapplot = gui.ecvmapaligned;
gui.endox = endox;
gui.endoy = endoy;
gui.epix = epix;
gui.epiy = epiy;
gui.hem = hem;
gui.rois = rois;
% gui.roimask = roimask;
gui.roiarea = roiarea;
gui.roiecv = roiecv;
gui.roiecvaligned = roiecvaligned;
gui.multislice = multislice;

%for the reporter the following is stored into the SET
for i = 1:length(preno)
  SET(preno(i)).ECV.preno = 1;
  SET(postno(i)).ECV.preno = 0;
  SET(preno(i)).ECV.postno = 0;
  SET(postno(i)).ECV.postno = 1;
  SET(postno(i)).ECV.hem = hem;
  SET(preno(i)).ECV.hem = hem;
  SET(postno(i)).ECV.roiecv = roiecv;
  SET(preno(i)).ECV.roiecv = roiecv;
  SET(preno(i)).ECV.roiarea = roiarea;
  SET(postno(i)).ECV.roiarea = roiarea;
  SET(postno(i)).ECV.multislice = multislice;
  SET(preno(i)).ECV.multislice = multislice;
  SET(postno(i)).ECV.rois = rois;
  SET(preno(i)).ECV.rois = rois;
end

gui.premaptf = premaptf;
gui.postmaptf = postmaptf;
gui.lvsegslicepre = lvsegslicepre;
gui.lvsegslicepost = lvsegslicepost;

%set the view max T1 map values sliders and editboxes
gui.viewmaxpre = 2000;
gui.viewmaxpost = 2000;
set(gui.handles.viewmaxpreslider,'value',gui.viewmaxpre);
set(gui.handles.viewmaxpreedit,'value',gui.viewmaxpre);
set(gui.handles.viewmaxpreedit,'String',gui.viewmaxpre);
set(gui.handles.viewmaxpostslider,'value',gui.viewmaxpost);
set(gui.handles.viewmaxpostedit,'value',gui.viewmaxpost);
set(gui.handles.viewmaxpostedit,'String',gui.viewmaxpost);

%plot T1 maps and ECV map
set(gui.handles.applyalignmentradiobutton,'value',1);
set(gui.handles.disablealignmentradiobutton,'value',0);
plotimages(multislice,lvsegslicepre,lvsegslicepost); 
plotecv(multislice,lvsegslicepre);
%plot LV segmentation and ROIs from the pre T1 map in alla images
set(gui.handles.showsegmentationradiobutton,'value',1);
set(gui.handles.hidesegmentationradiobutton,'value',0);
plotLVsegmentation(multislice,lvsegslicepre);
plotROI(multislice,lvsegslicepre);


%-------------------------------------------
function z = reshape2layout(im,rows,cols,sz)
%-------------------------------------------
%Convert a 3D array to an layout:ed image with cols, and rows

z = repmat(im(1),rows*sz(1),cols*sz(2));
loop = 1;
for slice = 1:sz(3)
  c = 1+mod(loop-1,cols);
  r = ceil(loop/cols);
  z((1+(r-1)*sz(1)):(r*sz(1)),...
    (1+(c-1)*sz(2)):(c*sz(2))) = im(:,:,slice);
  loop = loop+1;
end

%-------------------------------------------------------------------------
function [cropcenterpre,cropcenterpost,cropradius] = calccropsize(preno,postno,premaptf,postmaptf)
%-------------------------------------------------------------------------
%calculate crop size for the maps based on the LV segmentation
global SET

%find proper crop size from preT1map
preLVcenterX = mean(SET(preno).EpiX(:,premaptf,:));
preLVcenterY = mean(SET(preno).EpiY(:,premaptf,:));
preLVwidth = max([max(SET(preno).EpiX(:,premaptf,:))-min(SET(preno).EpiX(:,premaptf,:)) ...
  max(SET(preno).EpiY(:,premaptf,:))-min(SET(preno).EpiY(:,premaptf,:))]);
prexmin = max(1,round(preLVcenterX-1.2*preLVwidth));
prexmax = min(SET(preno).XSize,round(preLVcenterX+1.2*preLVwidth));
preymin = max(1,round(preLVcenterY-1.2*preLVwidth));
preymax = min(SET(preno).YSize,round(preLVcenterY+1.2*preLVwidth));
prexsize = prexmax-prexmin;
preysize = preymax-preymin;
%find proper crop size from postT1map
postLVcenterX = mean(SET(postno).EpiX(:,postmaptf,:));
postLVcenterY = mean(SET(postno).EpiY(:,postmaptf,:));
postLVwidth = max([max(SET(postno).EpiX(:,postmaptf,:))-min(SET(postno).EpiX(:,postmaptf,:)) ...
  max(SET(postno).EpiY(:,postmaptf,:))-min(SET(postno).EpiY(:,postmaptf,:))]);
postxmin = max(1,round(postLVcenterX-1.2*postLVwidth));
postxmax = min(SET(postno).XSize,round(postLVcenterX+1.2*postLVwidth));
postymin = max(1,round(postLVcenterY-1.2*postLVwidth));
postymax = min(SET(postno).YSize,round(postLVcenterY+1.2*postLVwidth));
postxsize = postxmax-postxmin;
postysize = postymax-postymin;

%find final crop size, the same for preT1map and postT1map
minxsize = min(prexsize,postxsize);
minysize = min(preysize,postysize);
if any(minxsize-prexsize<0)
  reducex = prexsize-minxsize;
  prexmin = prexmin+ceil(reducex/2);
  prexmax = prexmax-floor(reducex/2);
end
if any(minxsize-postxsize < 0)
  reducex = postxsize-minxsize;
  postxmin = postxmin+ceil(reducex/2);
  postxmax = postxmax-floor(reducex/2);
end
if any(minysize-preysize < 0)
  reducey = preysize-minysize;
  preymin = preymin+ceil(reducey/2);
  preymax = preymax-floor(reducey/2);
end
if any(minysize-postysize < 0)
  reducey = postysize-minysize;
  postymin = postymin+ceil(reducey/2);
  postymax = postymax-floor(reducey/2);
end

%exclude slices with no LV segmentation
prexmin(isnan(preLVcenterX)) = NaN;
prexmax(isnan(preLVcenterX)) = NaN;
preymin(isnan(preLVcenterX)) = NaN;
preymax(isnan(preLVcenterX)) = NaN;
postxmin(isnan(preLVcenterX)) = NaN;
postxmax(isnan(preLVcenterX)) = NaN;
postymin(isnan(preLVcenterX)) = NaN;
postymax(isnan(preLVcenterX)) = NaN;

% cropcenterpre = round([mean([prexmin prexmax]) mean([preymin preymax])]);
% cropcenterpost = round([mean([postxmin postxmax]) mean([postymin postymax])]);

cropcenterpre = round([preLVcenterX preLVcenterY]);
cropcenterpost = round([postLVcenterX postLVcenterY]);
cropradius = round((prexmax-prexmin)/2);

%-------------------------------------------------------------------------
function maxcropradius = calcsharedcropsize(cropcenterpre,cropcenterpost,cropradius,imagesizepre,imagesizepost)
%---------------------------------------------------------------------
%calculate shared crop size for the maps based on proposed crop sizes for
%each image

maxcropradius = max(cropradius);
xmin = min([cropcenterpre(1,:)-maxcropradius cropcenterpost(1,:)-maxcropradius]);
ymin = min([cropcenterpre(2,:)-maxcropradius cropcenterpost(2,:)-maxcropradius]);
xmaxpre = max(cropcenterpre(1,:)+maxcropradius); 
xmaxpost = max(cropcenterpost(1,:)+maxcropradius);
ymaxpre = max(cropcenterpre(2,:)+maxcropradius); 
ymaxpost = max(cropcenterpost(2,:)+maxcropradius);
if xmin<1 || ymin<1 || xmaxpre>min(imagesizepre(1,:)) || xmaxpost>min(imagesizepost(1,:)) || ...
    ymaxpre>min(imagesizepre(2,:)) || ymaxpost>min(imagesizepost(2,:))
  %reduce the radius so the crop is inside all images
  %find largest exceeding value
  exceeding = [xmin-1 ymin-1 min(imagesizepre(1,:))-xmaxpre min(imagesizepost(1,:))-xmaxpost min(imagesizepre(2,:))-ymaxpre min(imagesizepost(2,:))-ymaxpost];
  minval = min(exceeding);
  maxcropradius = maxcropradius+minval;
  disp('Coorecting common crop size');
end


%----------------------------------------------------------------
function [preT1map,postT1map,endox,endoy,epix,epiy,rois,roimask,roiarea] = ...
  cropmaps(preno,postno,premaptf,postmaptf,cropcenterpre,cropcenterpost,cropradius,slice)
%----------------------------------------------------------------
%crop the T1 maps based on the LV segmentation, so the two T1 maps have the
%same size afterwards

global SET

%crop images
preT1map = SET(preno).IM(cropcenterpre(1)-cropradius:cropcenterpre(1)+cropradius, ...
  cropcenterpre(2)-cropradius:cropcenterpre(2)+cropradius,premaptf,slice);
postT1map = SET(postno).IM(cropcenterpost(1)-cropradius:cropcenterpost(1)+cropradius, ...
  cropcenterpost(2)-cropradius:cropcenterpost(2)+cropradius,postmaptf,slice);
%crop LV segmentation
endox = SET(preno).EndoX(:,premaptf,slice)-(cropcenterpre(1)-cropradius)+1;
endoy = SET(preno).EndoY(:,premaptf,slice)-(cropcenterpre(2)-cropradius)+1;
epix = SET(preno).EpiX(:,premaptf,slice)-(cropcenterpre(1)-cropradius)+1;
epiy = SET(preno).EpiY(:,premaptf,slice)-(cropcenterpre(2)-cropradius)+1;
if SET(preno).RoiN > 1 %ROIs in addition to blood pool ROI
  rois = [];
  j = 1;
  for preroi = SET(preno).Roi
    if ~strcmp(preroi.Name,'Blood')
      if preroi.Z == slice
        temproi = preroi;
        temproi.X = temproi.X-(cropcenterpre(1)-cropradius)+1;
        temproi.Y = temproi.Y-(cropcenterpre(2)-cropradius)+1;
        rois = [rois temproi];
        roimask{j} = segment('createmask',...
          [size(preT1map,1) size(preT1map,2)],...
          preroi.Y(:,premaptf)-(cropcenterpre(2)-cropradius)+1,preroi.X(:,premaptf)-(cropcenterpre(1)-cropradius)+1);
        roiarea(j) = temproi.Area(:,premaptf);
        j = j+1;
      end
    end
  end
else
  rois = [];
  roimask = [];
  roiarea = [];
end

%----------------------------------
function setviewmax_Callback(button) %#ok<DEFNU>
%----------------------------------
%read in translation and rotation

global DATA

gui = DATA.GUI.ECVRegistration;

switch button
  case 'viewmaxpreslider'  %left-right translation
    viewmaxpre = mygetslider(gui.handles.viewmaxpreslider);
    set(gui.handles.viewmaxpreedit,'String',round(viewmaxpre));
    set(gui.handles.viewmaxpreedit,'value',viewmaxpre);
    gui.viewmaxpre = viewmaxpre;
  case 'viewmaxpreedit'
    viewmaxpre = min(max(str2num(mygetedit(gui.handles.viewmaxpreedit)),get(gui.handles.viewmaxpreslider,'min')),get(gui.handles.viewmaxpreslider,'max')); %#ok<ST2NM>
    set(gui.handles.viewmaxpreslider,'value',viewmaxpre);
    set(gui.handles.viewmaxpreedit,'String',round(viewmaxpre));
    set(gui.handles.viewmaxpreedit,'value',viewmaxpre);
    gui.viewmaxpre = viewmaxpre;
  case 'viewmaxpostslider'  %left-right translation
    viewmaxpost = mygetslider(gui.handles.viewmaxpostslider);
    set(gui.handles.viewmaxpostedit,'String',round(viewmaxpost));
    set(gui.handles.viewmaxpostedit,'value',viewmaxpost);
    gui.viewmaxpost = viewmaxpost;
  case 'viewmaxpostedit'
    viewmaxpost = min(max(str2num(mygetedit(gui.handles.viewmaxpostedit)),get(gui.handles.viewmaxpostslider,'min')),get(gui.handles.viewmaxpostslider,'max')); %#ok<ST2NM>
    set(gui.handles.viewmaxpostslider,'value',viewmaxpost);
    set(gui.handles.viewmaxpostedit,'String',round(viewmaxpost));
    set(gui.handles.viewmaxpostedit,'value',viewmaxpost);
    gui.viewmaxpost = viewmaxpost;
end
updateimages;

%-----------------------------
function plotimages(multislice,lvsegslicepre,lvsegslicepost)
%------------------------------
%Update plots of pre T1map and post T1map

global DATA

gui = DATA.GUI.ECVRegistration;

%set colormap
graymap = colormap(gray(256));


if length(gui.preno) == 1  && not(multislice) %plot only one image slice
  preT1mapimage = gui.preT1map{1};
  postT1mapalignedimage = gui.postT1mapplot{1};
  %set view max value of T1 maps
  preT1mapimage = min(preT1mapimage,2000);
  postT1mapalignedimage = min(postT1mapalignedimage,2000);
  %image value between 0 and 1
  preT1mapimage  = preT1mapimage *1/max(preT1mapimage (:));
  postT1mapalignedimage  = postT1mapalignedimage *1/max(postT1mapalignedimage (:));
  %set colormap
  viewimpre = spect.spectperfusionsegmentation('remapuint8',preT1mapimage,graymap);
  viewimpost = spect.spectperfusionsegmentation('remapuint8',postT1mapalignedimage,graymap);
  %plot images
  gui.handles.preimage = imagesc(viewimpre,'parent',gui.handles.preaxes);
  gui.handles.postimage = imagesc(viewimpost,'parent',gui.handles.postaxes);  
  axis(gui.handles.preaxes,'off','image')
  axis(gui.handles.postaxes,'off','image')
  
else
  if multislice %plot multiple slices
    nbrrows = length(lvsegslicepre);
    preT1map = gui.preT1map{1};
    postT1map = gui.postT1mapplot{1};
    for k=1:nbrrows
      preT1mapimage(:,:,k) = squeeze(preT1map(:,:,:,lvsegslicepre(k)));
      postT1mapalignedimage(:,:,k) = squeeze(postT1map(:,:,:,lvsegslicepost(k)));
    end
  else %if length(gui.preno) <= 3 %plot 2-3 image slices
    nbrrows = length(gui.preno);
    for k=1:nbrrows
      preT1mapimage(:,:,k) = gui.preT1map{k};
      postT1mapalignedimage(:,:,k) = gui.postT1mapplot{k};
    end
  end
  nbrcols = 1;
  szpre = size(preT1mapimage);
  viewimpre = reshape2layout(preT1mapimage,nbrrows,nbrcols,[szpre(1) szpre(2) nbrrows]);
  viewimpost = reshape2layout(postT1mapalignedimage,nbrrows,nbrcols,[szpre(1) szpre(2) nbrrows]);
  
  %set view max value of T1 maps
  viewimpre = min(viewimpre,gui.viewmaxpre);
  viewimpost = min(viewimpost,gui.viewmaxpost);
  
  %image value between 0 and 1
  viewimpre = viewimpre*1/max(viewimpre(:));
  viewimpost = viewimpost*1/max(viewimpost(:));

  %create visulaization images
  viewimpremap = spect.spectperfusionsegmentation('remapuint8',viewimpre,graymap);
  viewimpostmap = spect.spectperfusionsegmentation('remapuint8',viewimpost,graymap);
    
  %plot images
  gui.handles.preimage = imagesc(viewimpremap,'parent',gui.handles.preaxes);
  gui.handles.postimage = imagesc(viewimpostmap,'parent',gui.handles.postaxes);
  axis(gui.handles.preaxes,'off','image');
  axis(gui.handles.postaxes,'off','image');
    
end

%plot colorbars
precolorbartemp = linspace(0,1,255);
precolorbar = spect.spectperfusionsegmentation('remapuint8',flipud(precolorbartemp'),graymap);
gui.handles.precolorbarimage = imagesc(precolorbar,'parent',gui.handles.precolorbaraxes);
set(gui.handles.precolorbaraxes,'ytick',linspace(1,255,11),'YTickLabel',linspace(gui.viewmaxpre,0,11),'YAxisLocation','right','xtick',[]);
postcolorbartemp = linspace(0,1,255);
postcolorbar = spect.spectperfusionsegmentation('remapuint8',flipud(postcolorbartemp'),graymap);
gui.handles.postcolorbarimage = imagesc(postcolorbar,'parent',gui.handles.postcolorbaraxes);
set(gui.handles.postcolorbaraxes,'ytick',linspace(1,255,11),'YTickLabel',linspace(gui.viewmaxpost,0,11),'YAxisLocation','right','xtick',[]);


%---------------
function plotecv(multislice,lvsegslicepre)
%---------------
%Update plot of ECV map

global DATA

gui = DATA.GUI.ECVRegistration;

%set colormap
% hotmap = colormap(jet(256)); 
load('colormapecv');
hotmap = flipud(colormapecv)./256;

%plot ECV map

if length(gui.preno) == 1 && not(multislice)
  nbrrows = 1;
  viewimecv = spect.spectperfusionsegmentation('remapuint8',gui.ecvmapplot{1},hotmap);
  gui.handles.ecvimage = imagesc(viewimecv,'parent',gui.handles.ecvaxes);
  axis(gui.handles.ecvaxes,'off','image');
  set(gui.handles.ecvaxes,'clim',[0 1]);      
else
  if multislice
    nbrrows = length(lvsegslicepre);
    rows = lvsegslicepre;
    ecvmapplot = gui.ecvmapplot{1};
    for k=1:nbrrows
      ecvmapimage(:,:,k) = squeeze(ecvmapplot(:,:,:,lvsegslicepre(k)));
    end
  else%if length(gui.preno) <= 3
    nbrrows = length(gui.preno);
    rows = gui.preno;
    for k=1:nbrrows
      ecvmapimage(:,:,k) = gui.ecvmapplot{k};
    end
  end
  nbrcols = 1;
  szpre = size(ecvmapimage);
  gui.viewimecv = reshape2layout(ecvmapimage,nbrrows,nbrcols,[szpre(1) szpre(2) nbrrows]);
  viewimecvmap = spect.spectperfusionsegmentation('remapuint8',gui.viewimecv,hotmap);
  gui.handles.ecvimage = imagesc(viewimecvmap,'parent',gui.handles.ecvaxes);
  axis(gui.handles.ecvaxes,'off','image')
  set(gui.handles.ecvaxes,'clim',[0 1]);
end
  
%plot ecv colorbar
ecvcolorbartemp = linspace(0,1,255);
ecvcolorbar = spect.spectperfusionsegmentation('remapuint8',flipud(ecvcolorbartemp'),hotmap);
gui.handles.ecvcolorbarimage = imagesc(ecvcolorbar,'parent',gui.handles.ecvcolorbaraxes);
set(gui.handles.ecvcolorbaraxes,'ytick',linspace(1,255,11),'YTickLabel',linspace(100,0,11),'YAxisLocation','right','xtick',[]);

%print stack/slice numbers in the lower left corner of the ecv image
if nbrrows > 1
  for k = 1:nbrrows
    x = 0.005;
    y = 1-k/nbrrows+0.005;
    gui.textNO(k) = text(x,y,num2str(rows(k)),'FontSize',10,'Color','w','Units','normalized','VerticalAlignment','bottom','Parent',gui.handles.ecvaxes);
  end
end


%--------------------------
function plotLVsegmentation(multislice,lvsegslicepre)
%--------------------------
%plot LV segmentation in all images

global DATA

gui = DATA.GUI.ECVRegistration;

if length(gui.preno) == 1 && not(multislice)
  gui.handles.endoyView = gui.endoy{1};
  gui.handles.endoxView = gui.endox{1};
  gui.handles.epiyView = gui.epiy{1};
  gui.handles.epixView = gui.epix{1};
  
else
  %reshape LV segmentation to view
  zvisualization = 0;
  sz = size(gui.preT1map{1});
  sz = sz(1:2);
  nbrcols = 1;
  if multislice
    endox = gui.endox{1};
    endoy = gui.endoy{1};
    epix = gui.epix{1};
    epiy = gui.epiy{1};
    nbrrows = length(lvsegslicepre);
    for zloop = 1:nbrrows
      zvisualization = zvisualization+1;  %zloop-paneldivision(panelloop)+1;
      [xofs,yofs] = spect.spectperfusionsegmentation('calcoffset',zvisualization,nbrcols,sz);
      %Endocontour
      gui.handles.endoxView((1+(zvisualization-1)*(DATA.NumPoints+1)):(zvisualization*(DATA.NumPoints+1)-1)) = ...
        endox(:,1,lvsegslicepre(zloop))+xofs;
      gui.handles.endoyView((1+(zvisualization-1)*(DATA.NumPoints+1)):(zvisualization*(DATA.NumPoints+1)-1)) = ...
        endoy(:,1,lvsegslicepre(zloop))+yofs;
      %Epicontour
      gui.handles.epixView((1+(zvisualization-1)*(DATA.NumPoints+1)):(zvisualization*(DATA.NumPoints+1)-1)) = ...
        epix(:,1,lvsegslicepre(zloop))+xofs;
      gui.handles.epiyView((1+(zvisualization-1)*(DATA.NumPoints+1)):(zvisualization*(DATA.NumPoints+1)-1)) = ...
        epiy(:,1,lvsegslicepre(zloop))+yofs;
    end
  else  
    nbrrows = length(gui.preno);
    for zloop = 1:nbrrows
      zvisualization = zvisualization+1;  %zloop-paneldivision(panelloop)+1;
      [xofs,yofs] = spect.spectperfusionsegmentation('calcoffset',zvisualization,nbrcols,sz);
      %Endocontour
      gui.handles.endoxView((1+(zvisualization-1)*(DATA.NumPoints+1)):(zvisualization*(DATA.NumPoints+1)-1)) = ...
        gui.endox{zloop}+xofs;
      gui.handles.endoyView((1+(zvisualization-1)*(DATA.NumPoints+1)):(zvisualization*(DATA.NumPoints+1)-1)) = ...
        gui.endoy{zloop}+yofs;
      %Epicontour
      gui.handles.epixView((1+(zvisualization-1)*(DATA.NumPoints+1)):(zvisualization*(DATA.NumPoints+1)-1)) = ...
        gui.epix{zloop}+xofs;
      gui.handles.epiyView((1+(zvisualization-1)*(DATA.NumPoints+1)):(zvisualization*(DATA.NumPoints+1)-1)) = ...
        gui.epiy{zloop}+yofs;
    end
  end    
  gui.handles.endoxView(gui.handles.endoxView==0) = NaN;
  gui.handles.endoyView(gui.handles.endoyView==0) = NaN;
  gui.handles.epixView(gui.handles.epixView==0) = NaN;
  gui.handles.epiyView(gui.handles.epiyView==0) = NaN;
end

%plot LV segmentation
hold(gui.handles.preaxes,'on');
gui.handles.preendo = plot(gui.handles.preaxes,gui.handles.endoyView,gui.handles.endoxView,'r-');
gui.handles.preepi = plot(gui.handles.preaxes,gui.handles.epiyView,gui.handles.epixView,'g-');
hold(gui.handles.preaxes,'off');
hold(gui.handles.postaxes,'on');
gui.handles.postendo = plot(gui.handles.postaxes,gui.handles.endoyView,gui.handles.endoxView,'r-');
gui.handles.postepi = plot(gui.handles.postaxes,gui.handles.epiyView,gui.handles.epixView,'g-');
hold(gui.handles.postaxes,'off');
hold(gui.handles.ecvaxes,'on');
gui.handles.ecvendo = plot(gui.handles.ecvaxes,gui.handles.endoyView,gui.handles.endoxView,'w-');
gui.handles.ecvepi = plot(gui.handles.ecvaxes,gui.handles.epiyView,gui.handles.epixView,'w-');
hold(gui.handles.ecvaxes,'off');
   

%---------------
function plotROI(multislice,lvsegslicepre)
%---------------
%plot ROIs from the pre T1 map in alla images

global DATA

gui = DATA.GUI.ECVRegistration;

%present ECV result
AxesTables = [];
AxesTables.result = axestable(gui.handles.reportaxes);
AxesTables.result.backgroundcolor = [0.94 0.94 0.94]; %[0 0 0];%
AxesTables.result.fontcolor = [0 0 0]; % [1 1 1];%
AxesTables.result.fontsize = 10;
AxesTables.result.ystep = 22;
AxesTables.result.addTable(sprintf('ECV Result.  Hematocrit: %0.2f',gui.hem),5,4,[0.18 0.25 0.19 0.19 0.19]);
AxesTables.result.addKey('title','Label',[],{'Area[mm2]','Mean','Min','Max'});

zvisualization = 1;
sz = size(gui.preT1map{1});
sz = sz(1:2);
nbrcols = 1;
if multislice
  nbrrows = length(lvsegslicepre);
  roiareatemp = gui.roiarea{1};
  roiecvalignedtemp = gui.roiecvaligned{1};
  allroistemp = gui.rois{1};
else
  nbrrows = length(gui.preno);
end

for zloop = 1:nbrrows
  %ecv result string
  roiloop = 1;
  if multislice
    tabletext = sprintf('Slice %0.0f\n',lvsegslicepre(zloop));
    roiarea = roiareatemp{zloop};
    roiecv = roiecvalignedtemp{zloop};
    prerois = allroistemp{zloop};
    tfindex = 1;
  else
    tabletext = sprintf('Image stack %0.0f\n',zloop);
    roiarea = gui.roiarea{zloop};
    roiecv = gui.roiecvaligned{zloop};
    prerois = gui.rois{zloop};
    tfindex = zloop;
  end
  AxesTables.result.addTable(tabletext,5,4,[0.2 0.23 0.19 0.19 0.19]);
%   AxesTables.result.addKey('title2',tabletext,[],{[],[],[],[]});
  %ROI contour
  [xofs,yofs] = spect.spectperfusionsegmentation('calcoffset',zvisualization,nbrcols,sz);
  zvisualization = zvisualization+1;
  roixView = [];
  roiyView = [];
  roilabelX = [];
  roilabelY = [];
  for preroi = prerois
    roixView((1+(roiloop-1)*(DATA.NumPoints+1)):(roiloop*(DATA.NumPoints+1)-1)) = ...
      preroi.X(:,gui.premaptf(tfindex))+xofs;
    roiyView((1+(roiloop-1)*(DATA.NumPoints+1)):(roiloop*(DATA.NumPoints+1)-1)) = ...
      preroi.Y(:,gui.premaptf(tfindex))+yofs;
    roilabelX(roiloop) = mynanmean(preroi.X(:,gui.premaptf(tfindex))+xofs);
    roilabelY(roiloop) = mynanmean(preroi.Y(:,gui.premaptf(tfindex))+yofs);
    AxesTables.result.addKey(preroi.Name,preroi.Name,[],{dprintf('%0.2f',roiarea(roiloop)),round(100*mean(roiecv{roiloop})),round(100*min(roiecv{roiloop})),round(100*max(roiecv{roiloop}))});
    roiloop = roiloop+1;
  end
  
  %plot roi
  roixView(roixView==0) = NaN;
  roiyView(roiyView==0) = NaN;
  hold(gui.handles.preaxes,'on');
  gui.handles.preroiplot(zloop) = plot(gui.handles.preaxes,roiyView,roixView,'b-');
  gui.handles.preroilabel{zloop} = text(roilabelY,roilabelX,{prerois.Name}, ...
    'Parent',gui.handles.preaxes,'Color','b','FontSize',8,'HorizontalAlignment','center');
  hold(gui.handles.preaxes,'off');
  hold(gui.handles.postaxes,'on');
  gui.handles.postroiplot(zloop) = plot(gui.handles.postaxes,roiyView,roixView,'b-');
  hold(gui.handles.postaxes,'off');
  hold(gui.handles.ecvaxes,'on');
  gui.handles.ecvroiplot(zloop) = plot(gui.handles.ecvaxes,roiyView,roixView,'k-');
  hold(gui.handles.ecvaxes,'off');
  %store the roi plottig
  gui.handles.roixView{zloop} = roixView;
  gui.handles.roiyView{zloop} = roiyView;
  gui.handles.roilabels{zloop} = {prerois.Name};
end

% add table to axes
AxesTables.result.draw();
gui.AxesTables.result = AxesTables.result;


%--------------------
function updateimages
%--------------------
%Update plots of pre T1map and post T1map

global DATA

gui = DATA.GUI.ECVRegistration;

%set colormap
graymap = colormap(gray(256));

%plot only one image slice
if length(gui.preno) == 1 && not(gui.multislice) 
  preT1mapimage = gui.preT1map{1};
  postT1mapplotimage = gui.postT1mapplot{1};
  %set view max value of T1 maps
  preT1mapimage = min(preT1mapimage,gui.viewmaxpre);
  postT1mapplotimage = min(postT1mapplotimage,gui.viewmaxpost);
  %image value between 0 and 1
  preT1mapimage  = preT1mapimage *1/max(preT1mapimage (:));
  postT1mapplotimage  = postT1mapplotimage *1/max(postT1mapplotimage (:));
  %create visulaization images
  viewimpre = spect.spectperfusionsegmentation('remapuint8',preT1mapimage,graymap);
  viewimpost = spect.spectperfusionsegmentation('remapuint8',postT1mapplotimage,graymap);
  %plot images
  set(gui.handles.preimage,'cdata',viewimpre);
  set(gui.handles.postimage,'cdata',viewimpost);  
   
%plot multiple image slices  
else  
  if gui.multislice %plot multiple slices
    nbrrows = length(gui.lvsegslicepre);
    preT1map = gui.preT1map{1};
    postT1map = gui.postT1mapplot{1};
    for k=1:nbrrows
      preT1mapimage(:,:,k) = squeeze(preT1map(:,:,:,gui.lvsegslicepre(k)));
      postT1mapalignedimage(:,:,k) = squeeze(postT1map(:,:,:,gui.lvsegslicepost(k)));
    end
  else %if length(gui.preno) <= 3 %plot 2-3 image slices
    nbrrows = length(gui.preno);
    for k=1:nbrrows
      preT1mapimage(:,:,k) = gui.preT1map{k};
      postT1mapalignedimage(:,:,k) = gui.postT1mapplot{k};
    end
  end
  nbrcols = 1;
  szpre = size(preT1mapimage);
  viewimpre = reshape2layout(preT1mapimage,nbrrows,nbrcols,[szpre(1) szpre(2) nbrrows]);
  viewimpost = reshape2layout(postT1mapalignedimage,nbrrows,nbrcols,[szpre(1) szpre(2) nbrrows]);
  
  %set view max value of T1 maps
  viewimpre = min(viewimpre,gui.viewmaxpre);
  viewimpost = min(viewimpost,gui.viewmaxpost);
  
  %image value between 0 and 1
  viewimpre = viewimpre*1/max(viewimpre(:));
  viewimpost = viewimpost*1/max(viewimpost(:));

  %create visulaization images
  viewimpremap = spect.spectperfusionsegmentation('remapuint8',viewimpre,graymap);
  viewimpostmap = spect.spectperfusionsegmentation('remapuint8',viewimpost,graymap);
      
  %plot images
  set(gui.handles.preimage,'cdata',viewimpremap);
  set(gui.handles.postimage,'cdata',viewimpostmap);
 
end
  
%plot colorbars
set(gui.handles.precolorbaraxes,'YTickLabel',round(linspace(gui.viewmaxpre,0,11)));
set(gui.handles.postcolorbaraxes,'YTickLabel',round(linspace(gui.viewmaxpost,0,11)));


%-----------------
function updateecv
%-----------------
%Update plot of ECV map

global DATA

gui = DATA.GUI.ECVRegistration;

%set colormap
hotmap = colormap(jet(256)); 

%plot ECV map

if length(gui.preno) == 1 && not(gui.multislice)
  viewimecvmap = spect.spectperfusionsegmentation('remapuint8',gui.ecvmapplot{1},hotmap);
  set(gui.handles.ecvimage,'cdata',viewimecvmap);
      
else
  if gui.multislice
    nbrrows = length(gui.lvsegslicepre);
    rows = gui.lvsegslicepre;
    ecvmapplot = gui.ecvmapplot{1};
    for k=1:nbrrows
      ecvmapimage(:,:,k) = squeeze(ecvmapplot(:,:,:,gui.lvsegslicepre(k)));
    end
  else%if length(gui.preno) <= 3
    nbrrows = length(gui.preno);
    rows = gui.preno;
    for k=1:nbrrows
      ecvmapimage(:,:,k) = gui.ecvmapplot{k};
    end
  end
  nbrcols = 1;
  szpre = size(ecvmapimage);
  gui.viewimecv = reshape2layout(ecvmapimage,nbrrows,nbrcols,[szpre(1) szpre(2) nbrrows]);
  viewimecvmap = spect.spectperfusionsegmentation('remapuint8',gui.viewimecv,hotmap);
  set(gui.handles.ecvimage,'cdata',viewimecvmap);
end


%-----------------------
function updateroiresult
%-----------------------
%update the ROI ECV result based on alignment on or alignment off

global DATA

gui = DATA.GUI.ECVRegistration;
AxesTables.result = gui.AxesTables.result;

if gui.multislice
  nbrrows = length(gui.lvsegslicepre);
  roiareatemp = gui.roiarea{1};
  if get(gui.handles.applyalignmentradiobutton,'Value')
    roiecvalignedtemp = gui.roiecvaligned{1};
  else
    roiecvalignedtemp = gui.roiecv{1};
  end
  allroistemp = gui.rois{1};
else
  nbrrows = length(gui.preno);
end

for zloop = 1:nbrrows
  %ecv result string
  roiloop = 1;
  if gui.multislice
%     string = [string dprintf('Slice %0.0f\n',gui.lvsegslicepre(zloop))];
    roiarea = roiareatemp{zloop};
    roiecv = roiecvalignedtemp{zloop};
    prerois = allroistemp{zloop};
  else
%     string = [string dprintf('Image stack %0.0f\n',zloop)];
    roiarea = gui.roiarea{zloop};
    if get(gui.handles.applyalignmentradiobutton,'Value')
      roiecv = gui.roiecvaligned{zloop};
    else
      roiecv = gui.roiecv{zloop};
    end
    prerois = gui.rois{zloop};
  end
  for preroi = prerois
    updatestruct = {dprintf('%0.2f',roiarea(roiloop)),round(100*mean(roiecv{roiloop})),round(100*min(roiecv{roiloop})),round(100*max(roiecv{roiloop}))};
    AxesTables.result.updateKey(preroi.Name,updatestruct,true);
    roiloop = roiloop+1;
  end
end

% add table to axes
AxesTables.result.draw();
gui.AxesTables.result = AxesTables.result;


%----------------------------------
function showsegmentations_Callback %#ok<DEFNU>
%----------------------------------
%Show ROI and LV segmentation in all images

global DATA
  
gui = DATA.GUI.ECVRegistration;

set(gui.handles.showsegmentationradiobutton,'value',1);
set(gui.handles.hidesegmentationradiobutton,'value',0);
%LV segmentation
set(gui.handles.preendo,'xdata',gui.handles.endoyView);
set(gui.handles.preendo,'ydata',gui.handles.endoxView);
set(gui.handles.preepi,'xdata',gui.handles.epiyView);
set(gui.handles.preepi,'ydata',gui.handles.epixView);
set(gui.handles.postendo,'xdata',gui.handles.endoyView);
set(gui.handles.postendo,'ydata',gui.handles.endoxView);
set(gui.handles.postepi,'xdata',gui.handles.epiyView);
set(gui.handles.postepi,'ydata',gui.handles.epixView);
set(gui.handles.ecvendo,'xdata',gui.handles.endoyView);
set(gui.handles.ecvendo,'ydata',gui.handles.endoxView);
set(gui.handles.ecvepi,'xdata',gui.handles.epiyView);
set(gui.handles.ecvepi,'ydata',gui.handles.epixView);
%ROI segmentation
if gui.multislice
  slices = gui.lvsegslicepre';
else
  slices = 1:length(gui.preT1map);
end
for zloop = slices
  set(gui.handles.preroiplot(zloop),'xdata',gui.handles.roiyView{zloop});
  set(gui.handles.preroiplot(zloop),'ydata',gui.handles.roixView{zloop});
  set(gui.handles.postroiplot(zloop),'xdata',gui.handles.roiyView{zloop});
  set(gui.handles.postroiplot(zloop),'ydata',gui.handles.roixView{zloop});
  set(gui.handles.ecvroiplot(zloop),'xdata',gui.handles.roiyView{zloop});
  set(gui.handles.ecvroiplot(zloop),'ydata',gui.handles.roixView{zloop});
  roilabel = gui.handles.roilabels{zloop};
  set(gui.handles.preroilabel{zloop},{'String'},roilabel');
end


%----------------------------------
function hidesegmentations_Callback %#ok<DEFNU>
%----------------------------------
%Hide ROI and LV segmentation in all images

global DATA
  
gui = DATA.GUI.ECVRegistration;

set(gui.handles.showsegmentationradiobutton,'value',0);
set(gui.handles.hidesegmentationradiobutton,'value',1);
%LV segmentation
set(gui.handles.preendo,'xdata',[]);
set(gui.handles.preendo,'ydata',[]);
set(gui.handles.preepi,'xdata',[]);
set(gui.handles.preepi,'ydata',[]);
set(gui.handles.postendo,'xdata',[]);
set(gui.handles.postepi,'ydata',[]);
set(gui.handles.postepi,'xdata',[]);
set(gui.handles.postendo,'ydata',[]);
set(gui.handles.ecvendo,'xdata',[]);
set(gui.handles.ecvendo,'ydata',[]);
set(gui.handles.ecvepi,'xdata',[]);
set(gui.handles.ecvepi,'ydata',[]);
%ROI segmentation
if gui.multislice
  slices = gui.lvsegslicepre';
else
  slices = 1:length(gui.preT1map);
end
for zloop = slices
  set(gui.handles.preroiplot(zloop),'xdata',[]);
  set(gui.handles.preroiplot(zloop),'ydata',[]);
  set(gui.handles.postroiplot(zloop),'xdata',[]);
  set(gui.handles.postroiplot(zloop),'ydata',[]);
  set(gui.handles.ecvroiplot(zloop),'xdata',[]);
  set(gui.handles.ecvroiplot(zloop),'ydata',[]);
  set(gui.handles.preroilabel{zloop},'String',[]);
end


%-----------------------------
function alignmentoff_Callback %#ok<DEFNU>
%-----------------------------
%Turn off the alignment of the post T1map

global DATA

gui = DATA.GUI.ECVRegistration;

set(gui.handles.applyalignmentradiobutton,'value',0);
set(gui.handles.disablealignmentradiobutton,'value',1);

gui.postT1mapplot = gui.postT1map;
gui.ecvmapplot = gui.ecvmap;
updateimages;
updateecv;
updateroiresult;


%----------------------------
function alignmenton_Callback %#ok<DEFNU>
%----------------------------
%Turn off the alignment of the post T1map

global DATA

gui = DATA.GUI.ECVRegistration;

set(gui.handles.applyalignmentradiobutton,'value',1);
set(gui.handles.disablealignmentradiobutton,'value',0);

gui.postT1mapplot = gui.postT1mapaligned;
gui.ecvmapplot = gui.ecvmapaligned;
updateimages;
updateecv;
updateroiresult;


%-----------------------
function export_Callback %#ok<DEFNU>
%-----------------------
%Exports the ECV result to clipboard

global DATA SET

gui = DATA.GUI.ECVRegistration;

% output = cell(19,2+strlen+rstlen);
output{1,1} = 'Patient name';
output{1,2} = SET(gui.preno(1)).PatientInfo.Name;
output{2,1} = 'Hematocrit';
output{2,2} = gui.hem;
output{3,1} = 'ECV ROI';
output{3,2} = 'ECV Area';
output{3,3} = 'ECV Mean';
output{3,4} = 'ECV Min';
output{3,5} = 'ECV Max';

%ecv ROI result
row = 4;
if gui.multislice
  allrois = gui.rois{1};
  allroiarea = gui.roiarea{1};
  if get(gui.handles.applyalignmentradiobutton,'Value')
    allroiecv = gui.roiecvaligned{1};
  else
    allroiecv = gui.roiecv{1};
  end
  for zloop = 1:length(allrois)
    output{row,1} = dprintf('Slice %s',get(gui.textNO(zloop),'String'));
    roiloop = 1;
    roiarea = allroiarea{zloop};
    roiecv = allroiecv{zloop};
    row = row+1;
    for preroi = allrois{zloop}
      %roi result string
      output{row,1} = preroi.Name;
      output{row,2} = roiarea(roiloop);
      output{row,3} = 100*mean(roiecv{roiloop});
      output{row,4} = 100*min(roiecv{roiloop});
      output{row,5} = 100*max(roiecv{roiloop});
      roiloop = roiloop+1;
      row = row+1;
    end
  end
else
  for zloop = 1:length(gui.preno)
    output{row,1} = dprintf('Image stack %0.0f',zloop);
    roiloop = 1;
    roiarea = gui.roiarea{zloop};
    if get(gui.handles.applyalignmentradiobutton,'Value')
      roiecv = gui.roiecvaligned{zloop};
    else
      roiecv = gui.roiecv{zloop};
    end
    row = row+1;
    for preroi = gui.rois{zloop}
      %roi result string
      output{row,1} = preroi.Name;
      output{row,2} = roiarea(roiloop);
      output{row,3} = 100*mean(roiecv{roiloop});
      output{row,4} = 100*min(roiecv{roiloop});
      output{row,5} = 100*max(roiecv{roiloop});
      roiloop = roiloop+1;
      row = row+1;
    end
  end
end
  
segment('cell2clipboard',output);


%-----------------------
function create_Callback %#ok<DEFNU>
%-----------------------
%Create an ECV map image stack and close the ECV GUI

global DATA SET NO

gui = DATA.GUI.ECVRegistration;

%create ECV map image stack
for zloop = 1:length(gui.preno)
  nbr = length(SET)+1;
  NO = nbr;
  no = nbr;
  SET(no) = SET(gui.preno(1));
  %fix the time resolution
  if SET(no).TSize>1
    DATA.Silent = 1;
    ind = 1;
    tools('removetimeframes',ind);
    DATA.Silent = 0;
    [SET(no).TSize,SET(no).OrgTSize,SET(no).StartSlice,SET(no).EndSlice, ...
      SET(no).StartAnalysis,SET(no).EndAnalysis,SET(no).CurrentTimeFrame] = deal(1);
    SET(no).TIncr = 0;
  end
  %fix image
  SET(no).IM = max(0,min(1,gui.ecvmap{zloop}));
  SET(no).IntensityScaling = 100;
  SET(no).IntensityOffset = 0;
  SET(no).OrgZSize = SET(no).ZSize;
  SET(no).OrgTSize = SET(no).TSize;
  SET(no).ImageType = 'ECV map';
  SET(no).SeriesDescription = 'ECV map';
  [SET(no).XSize,SET(no).OrgXSize] = deal(size(gui.ecvmap{zloop},2));
  [SET(no).YSize,SET(no).OrgYSzie] = deal(size(gui.ecvmap{zloop},1));
  [SET(no).CenterX,SET(no).Mmode.X] = deal(round(SET(no).XSize/2));
  [SET(no).CenterY,SET(no).Mmode.Y] = deal(round(SET(no).YSize/2));
  SET(no).Linked = no;
  annotationpoint('pointclearall'); %erase annotation points
  %fix LV segmentation
  SET(no).EndoX = gui.endox{zloop};
  SET(no).EndoY = gui.endoy{zloop};
  SET(no).EpiX = gui.epix{zloop};
  SET(no).EpiY = gui.epiy{zloop};
  %fix ROI segmentation  
  if gui.multislice
    allrois = gui.rois{zloop};    
    SET(no).Roi = [];
    for sliceloop = 1:length(allrois)
      rois = allrois{sliceloop};
      SET(no).Roi = [SET(no).Roi rois];
    end
  else
    SET(no).Roi = gui.rois{zloop};
  end
  SET(no).RoiN = length(SET(no).Roi);
  SET(no).RoiCurrent = min(SET(no).RoiCurrent,SET(no).RoiN);
  for rloop=1:SET(no).RoiN
    SET(no).Roi(rloop).X = SET(no).Roi(rloop).X(:,gui.premaptf);
    SET(no).Roi(rloop).Y = SET(no).Roi(rloop).Y(:,gui.premaptf);
    SET(no).Roi(rloop).T = 1;
  end
  
  tools('setcolormap_Callback','ecv',no);
  %erase MaR and scar data
  SET(no).MaR = [];
  SET(no).Scar = [];
end

segment('updatevolume'); %Included in next callback??
segment('viewrefreshall_Callback');
force = true;
segment('switchtoimagestack',no,force);
tools('setcolormap_Callback','ecv',no);
% drawfunctions('drawimageno');
% drawfunctions('drawsliceno');
% drawfunctions('drawall',length(DATA.ViewPanels));


%----------------------
function close_Callback
%----------------------
%close the ECV GUI

global DATA

try
  DATA.GUI.ECVRegistration = close(DATA.GUI.ECVRegistration);
catch %#ok<CTCH>
  close(gcf)
end
