function setstruct = loadfieldhelper(setstruct)
%This fcn fixes backward compability when a .mat file is loaded

global DATA SET 

if nargin < 1
  setstruct = SET;
end

if not(isfield(setstruct, 'SpectSpecialTag'))
  setstruct(1).SpectSpecialTag = [];
end

if not(isfield(setstruct,'Stress'))
  setstruct(1).Stress = [];
end

if not(isfield(setstruct,'Scar'))
  setstruct(1).Scar = [];
end

for no = 1:length(setstruct)
   if ~isempty(setstruct(no).Scar)
       
     %removed scar segmentation method SCJ (replaced by weighted)
     if isequal(setstruct(no).Scar.Mode,'new')
       setstruct(no).Scar.Mode = 'weighted';
     end

    %UseWeighting field
    if ~isfield(setstruct(no).Scar,'UseWeighting')
      if isequal(setstruct(no).Scar.Mode,'weighted')
        setstruct(no).Scar.UseWeighting = true;
      else
        setstruct(no).Scar.UseWeighting = false;
      end
    end
        
    %NoReflow field
    if ~isfield(setstruct(no).Scar,'NoReflow')
      try
        setstruct(no).Scar.NoReflow = false(size(setstruct(no).Scar.Auto));
      catch
      end
    end
    
    %MOPercentage field
    if ~isfield(setstruct(no).Scar,'MOPercentage')
      setstruct(no).Scar.MOPercentage = NaN;
    end
    
    %GreyZone field
    if ~isfield(setstruct(no).Scar,'GreyZone')
      setstruct(no).Scar.GreyZone = [];
      setstruct(no).Scar.GreyZone.map = [];
      setstruct(no).Scar.GreyZone.limit = 0.5;
    end
    
    if ~isfield(setstruct(no).Scar,'MR')
      setstruct(no).Scar.MR.Rotation=0;
      setstruct(no).Scar.MR.Artery=[];
    end

   end
   
end
          
if ~isfield(setstruct,'EndoInterpOngoing')
  for no = 1:length(setstruct)    
     setstruct(no).EndoInterpOngoing = false;
  end
end
if ~isfield(setstruct,'EpiInterpOngoing')
  for no = 1:length(setstruct)    
     setstruct(no).EpiInterpOngoing = false;
  end
end
if ~isfield(setstruct,'RVEndoInterpOngoing')
  for no = 1:length(setstruct)    
     setstruct(no).RVEndoInterpOngoing = false;
  end
end
if ~isfield(setstruct,'RVEpiInterpOngoing')
  for no = 1:length(setstruct)    
     setstruct(no).RVEpiInterpOngoing = false;
  end
end

if not(isfield(setstruct,'RoiN'))
  for no = 1:length(setstruct)
    roi('roireset',no);    
  end
end

if not(isfield(setstruct,'CenterX'))
  for no=1:length(setstruct)
    setstruct(no).CenterX = setstruct(no).XSize/2;
  end
end

if not(isfield(setstruct,'CenterY'))
  for no=1:length(setstruct)
    setstruct(no).CenterY = setstruct(no).YSize/2;
  end
end

if not(isfield(setstruct,'Roi'))&&isfield(setstruct,'RoiX')
  for no=1:length(setstruct)
    if setstruct(no).RoiN==0||isempty(setstruct(no).RoiX)
      setstruct(no).Roi = roi('roireset');
    else
      for loop=1:setstruct(no).RoiN
        setstruct(no).Roi(loop).X = setstruct(no).RoiX(:,:,loop);
        setstruct(no).Roi(loop).Y = setstruct(no).RoiY(:,:,loop);
        setstruct(no).Roi(loop).T = 1:setstruct(no).TSize;
        setstruct(no).Roi(loop).Z = setstruct(no).RoiZ(loop);
        setstruct(no).Roi(loop).Name = setstruct(no).RoiName{loop};
        setstruct(no).Roi(loop).Sign = setstruct(no).RoiSign(loop);
        setstruct(no).Roi(loop).LineSpec = setstruct(no).RoiLineSpec{loop};
        %update mean area std in case of cropped by border
        %[meanarea,area] = calcfunctions('calcroiarea',no,loop);
        %meanarea = NaN;
        area = NaN;
        setstruct(no).Roi(loop).Area = area;
        %[m,sd]=calcfunctions('calcroiintensity',no,loop);
        m = 0;
        sd = 0;
        setstruct(no).Roi(loop).Mean = m;
        setstruct(no).Roi(loop).StD = sd;
      end
    end
  end
  setstruct=rmfield(setstruct,'RoiX');
  setstruct=rmfield(setstruct,'RoiY');
  setstruct=rmfield(setstruct,'RoiZ');
  setstruct=rmfield(setstruct,'RoiName');
  setstruct=rmfield(setstruct,'RoiSign');
  setstruct=rmfield(setstruct,'RoiLineSpec');
elseif isfield(setstruct,'Roi')
  for no = 1:length(setstruct)
    if (numel(setstruct(no).Roi) == 1 && isequal(setstruct(no).Roi.Name,{})) || isempty(setstruct(no).Roi)
      setstruct(no).Roi = roi('roireset');
    end
  end
end
if isfield(setstruct,'Roi')
  for no = 1:length(setstruct)
    if not(isfield(setstruct(no).Roi,'Area'))
      %calculate mean area and std
      for loop=1:setstruct(no).RoiN
        [~,area] = calcfunctions('calcroiarea',no,loop);
        setstruct(no).Roi(loop).Area = area;
        [m,sd] = calcfunctions('calcroiintensity',no,loop);
        setstruct(no).Roi(loop).Mean = m;
        setstruct(no).Roi(loop).StD = sd;
      end
    end
  end
end
%RoiCurrent field check
for no = 1:length(setstruct)
  if isequal(setstruct(no).RoiCurrent,0) || isequal(setstruct(no).RoiN,0)
    setstruct(no).RoiCurrent = [];
  end
  if setstruct(no).RoiN>0 && isempty(setstruct(no).RoiCurrent)
    setstruct(no).RoiCurrent = 1;
  end
end

%Respiratory
if not(isfield(setstruct,'Respiratory'))
  for no = 1:length(setstruct)
    setstruct(no).Respiratory = [];
  end
end

%New General Pen
if not(isfield(setstruct,'GeneralPenX'))
  for no = 1:length(setstruct)
    setstruct(no).GeneralPenX = [];
    setstruct(no).GeneralPenY = [];
  end
end

%New General Pen Interp
if not(isfield(setstruct,'GeneralPenInterpX'))
  for no = 1:length(setstruct)
    setstruct(no).GeneralPenInterpX = [];
    setstruct(no).GeneralPenInterpY = [];
  end
end

%Start analysis field check
if not(isfield(setstruct,'StartAnalysis'))
  for no = 1:length(setstruct)
    setstruct(no).StartAnalysis = 1;
    setstruct(no).EndAnalysis = setstruct(no).TSize;
  end
end
for no = 1:length(setstruct)
  if isempty(setstruct(no).StartAnalysis)
    setstruct(no).StartAnalysis = 1;
    setstruct(no).EndAnalysis = setstruct(no).TSize;    
  end
end

%Resolution field check
if not(isfield(setstruct,'ResolutionX'))
  for no = 1:length(setstruct)
    setstruct(no).ResolutionX = setstruct(no).Resolution;
    setstruct(no).ResolutionY = setstruct(no).Resolution;
  end
end
for no = 1:length(setstruct)
  if isempty(setstruct(no).ResolutionX) && isfield(setstruct(no),'Resolution')
    setstruct(no).ResolutionX = setstruct(no).Resolution;
    setstruct(no).ResolutionY = setstruct(no).Resolution;    
  end  
end

%Endocenter field check
if not(isfield(setstruct,'EndoCenter'))
  for no = 1:length(setstruct)
    setstruct(no).EndoCenter = DATA.Pref.EndoCenter;
  end
end
for no = 1:length(setstruct)
  if isempty(setstruct(no).EndoCenter)
    setstruct(no).EndoCenter = DATA.Pref.EndoCenter;
  end
end

if not(isfield(setstruct,'NormalZoomState'))
  for no = 1:length(setstruct)
    setstruct(no).NormalZoomState = [];
  end   
end

if not(isfield(setstruct,'MontageZoomState'))
  for no = 1:length(setstruct)
    setstruct(no).MontageZoomState = [];
  end
end

if not(isfield(setstruct,'MontageRowZoomState'))
  for no = 1:length(setstruct)
    setstruct(no).MontageRowZoomState = [];
  end
end

if not(isfield(setstruct,'MontageFitZoomState'))
  for no = 1:length(setstruct)
    setstruct(no).MontageFitZoomState = [];
  end
end

if not(isfield(setstruct,'Rotated'))
  for no = 1:length(setstruct)
    setstruct(no).Rotated = false;
  end
end

if not(isfield(setstruct,'Flow'))
  for no = 1:length(setstruct)
    setstruct(no).Flow = [];
  end
end

for no = 1:length(setstruct)
  if not(isempty(setstruct(no).Flow))
    if not(isfield(setstruct(no).Flow,'PhaseCorrTimeResolved'))    
      setstruct(no).Flow.PhaseCorrTimeResolved = false;
    end
  end
end
% check if number of ROIs corresponds to the length of ROI stuct
for no = 1:length(setstruct)
  numrois = length(setstruct(no).Roi);
  if setstruct(no).RoiN ~= 0 && setstruct(no).RoiN ~= numrois
    for roiloop = 1:numrois
      % check
      if isempty(setstruct(no).Roi(roiloop).X) || isempty(setstruct(no).Roi(roiloop).Y)
        setstruct(no).Roi(roiloop) = [];
      end
    end
  end
end

%This addition Writes the flow area if it exists over the area parameter
for no = 1:length(setstruct)
  if not(isempty(setstruct(no).Flow))
    if isfield(setstruct(no).Flow,'Result')
        for roiloop = 1:length(setstruct(no).Flow.Result)
          if roiloop <= setstruct(no).RoiN
            if isfield(setstruct(no).Flow.Result(roiloop),'area')
              lengthflowresult = nnz(setstruct(no).Flow.Result(roiloop).area);
              if lengthflowresult ~= setstruct(no).TSize
                % correct results if their length without zero values is not the same as TSize                
                numtimeframes = setstruct(no).TSize;
                setstruct(no).Flow.Result(roiloop).area = rebuildflowresults(setstruct(no).Flow.Result(roiloop).area,numtimeframes);
                setstruct(no).Flow.Result(roiloop).velmean = rebuildflowresults(setstruct(no).Flow.Result(roiloop).velmean,numtimeframes);
                setstruct(no).Flow.Result(roiloop).velstd = rebuildflowresults(setstruct(no).Flow.Result(roiloop).velstd,numtimeframes);
                setstruct(no).Flow.Result(roiloop).velmax = rebuildflowresults(setstruct(no).Flow.Result(roiloop).velmax,numtimeframes);
                setstruct(no).Flow.Result(roiloop).velmin = rebuildflowresults(setstruct(no).Flow.Result(roiloop).velmin,numtimeframes);
                setstruct(no).Flow.Result(roiloop).kenergy = rebuildflowresults(setstruct(no).Flow.Result(roiloop).kenergy,numtimeframes);
                setstruct(no).Flow.Result(roiloop).netflow = rebuildflowresults(setstruct(no).Flow.Result(roiloop).netflow,numtimeframes);
                setstruct(no).Flow.Result(roiloop).posflow = rebuildflowresults(setstruct(no).Flow.Result(roiloop).posflow,numtimeframes);
                setstruct(no).Flow.Result(roiloop).negflow = rebuildflowresults(setstruct(no).Flow.Result(roiloop).negflow,numtimeframes);
              end
              setstruct(no).Roi(roiloop).Area = setstruct(no).Flow.Result(roiloop).area;
            end
          else
            % there is no corresponding ROI for the result -> delete result
            setstruct(no).Flow.Result(roiloop) = [];
          end
        end
    end
  end
end

if not(isfield(setstruct,'VENC'))
  for no = 1:length(setstruct)
    setstruct(no).VENC = 0;
  end
end

if not(isfield(setstruct,'OrigFileName'))
  for no = 1:length(setstruct)
    setstruct(no).OrigFileName = setstruct(no).FileName;
  end
end

% JU Measure now time dependent
if not(isfield(setstruct,'Measure'))
  for no = 1:length(setstruct)
    setstruct(no).Measure = [];
  end
else
  for no = 1:length(setstruct)
      if (isfield(setstruct(no).Measure,'X')&&not(isfield(setstruct(no).Measure,'T')))
          for loop=1:length(setstruct(no).Measure)
            setstruct(no).Measure(loop).T = NaN;
          end
      end
      if isfield(setstruct(no).Measure,'Z') 
        for loop=1:length(setstruct(no).Measure)
          if numel(setstruct(no).Measure(loop).Z) == 1
            setstruct(no).Measure(loop).Z = setstruct(no).Measure(loop).Z*[1;1];
          end
        end
      end
  end
end

if not(isfield(setstruct,'Report'))
  for no = 1:length(setstruct)
    setstruct(no).Report = [];
  end
end

for no = 1:length(setstruct)
  if ~isfield(setstruct(no),'RV')
    setstruct(no).RV = [];
  end
end

if not(isfield(setstruct,'LevelSet'))
  for no = 1:length(setstruct)
    setstruct(no).LevelSet = [];
  end
end

%Orgsize field check
if not(isfield(setstruct,'OrgXSize'))
  for no = 1:length(setstruct)
    setstruct(no).OrgXSize = 256;
    setstruct(no).OrgYSize = 256;
    setstruct(no).OrgTSize = setstruct(no).TSize;
    setstruct(no).OrgZSize = setstruct(no).ZSize;
  end
end

%Orgresolution field check
if not(isfield(setstruct,'OrgRes'))
  for no = 1:length(setstruct)
    setstruct(no).OrgRes = [setstruct(no).ResolutionX,...
        setstruct(no).ResolutionY,setstruct(no).SliceThickness + ...
        setstruct(no).SliceGap,setstruct(no).TIncr];
  end
end

for no = 1:length(setstruct)
    if isempty(setstruct(no).OrgRes)
        setstruct(no).OrgRes = [setstruct(no).ResolutionX,...
            setstruct(no).ResolutionY,setstruct(no).SliceThickness + ...
            setstruct(no).SliceGap,setstruct(no).TIncr];
    end
end

for no = 1:length(setstruct)
  if isempty(setstruct(no).OrgXSize)
    setstruct(no).OrgXSize = 256;
    setstruct(no).OrgYSize = 256;
    setstruct(no).OrgTSize = setstruct(no).TSize;
    setstruct(no).OrgZSize = setstruct(no).ZSize;
  end
end

if not(isfield(setstruct,'GEVENCSCALE'))
  for no = 1:length(setstruct)
    setstruct(no).GEVENCSCALE = 0;
  end
end

%Check if there are any autolongaxis, and ask to reset
for no = 1:length(setstruct)
  if ~isfield(setstruct(no), 'AutoLongaxis') || isempty(setstruct(no).AutoLongaxis) 
    setstruct(no).AutoLongaxis = false;
  end
end
  
autolongaxis = cat(1,setstruct.AutoLongaxis) & (cat(1,setstruct.TSize)>1);
autolongaxis = sum(autolongaxis);
if autolongaxis
  if true %yesno('Old longaxis compensation is active in the loaded file. This may compromise quantified values unless you know exactly what you are doing. Do you want to disable this? (Yes=recommended)');
    for no = 1:length(setstruct)
      setstruct(no).AutoLongaxis = false;
      setstruct(no).Longaxis = 1;
    end
  end
end

if not(isfield(setstruct,'Longaxis'))
  for no = 1:length(setstruct)
    setstruct(no).Longaxis = 1; %Remove this from default, EH: 2017-06-19, 1=>first choice in listbox => 0. Was 1 before as well...
  end
end

for no = 1:length(setstruct)
  if isempty(setstruct(no).Longaxis)
    setstruct(no).Longaxis = 1; %Remove this from default, EH: 2017-06-19, 1=>first choice in listbox => 0. Was 1 before as well...
  end
end

if not(isfield(setstruct,'AutoLongaxis'))
  for no = 1:length(setstruct)
    setstruct(no).AutoLongaxis = false; %Remove this from default, EH: 2017-06-19
  end
else
  for no = 1:length(setstruct)
    if isempty(setstruct(no).AutoLongaxis)
      setstruct(no).AutoLongaxis = false; %Remove this from default, EH: 2017-06-19
    end
  end
end


if not(isfield(setstruct,'Strain'))
  for no = 1:length(setstruct)
    setstruct(no).Strain = [];
  end
end

if not(isfield(setstruct,'StrainTagging'))
  for no = 1:length(setstruct)
    setstruct(no).StrainTagging = [];
  end
end

%RV field check
if not(isfield(setstruct,'RVEndoX'))
  for no=1:length(setstruct)
    setstruct(no).RVEndoX = [];
    setstruct(no).RVEndoY = [];
    setstruct(no).RVEpiX = [];
    setstruct(no).RVEpiY = [];
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).RVEndoX)
    setstruct(no).RVEndoX = [];
    setstruct(no).RVEndoY = [];
  end
  if isempty(setstruct(no).RVEpiX)
    setstruct(no).RVEpiX = [];
    setstruct(no).RVEpiY = [];
  end
end

if not(isfield(setstruct,'RVEndoPinX'))
  for no=1:length(setstruct)
    setstruct(no).RVEndoPinX = [];
    setstruct(no).RVEndoPinY = [];
    setstruct(no).RVEndoPinXView = [];
    setstruct(no).RVEndoPinYView = [];    
    setstruct(no).RVEpiPinX = [];
    setstruct(no).RVEpiPinY = [];
    setstruct(no).RVEpiPinXView = [];
    setstruct(no).RVEpiPinYView = [];        
  end
end

if not(isfield(setstruct,'EndoInterpX'))
  for no=1:length(setstruct)
    setstruct(no).EndoInterpX = [];
    setstruct(no).EndoInterpY = [];
    setstruct(no).EpiInterpX = [];
    setstruct(no).EpiInterpY = [];
    setstruct(no).RVEndoInterpX = [];
    setstruct(no).RVEndoInterpY = [];
    setstruct(no).RVEpiInterpX = [];
    setstruct(no).RVEpiInterpY = [];
  end
end

if not(isfield(setstruct,'RVV'))
  for no = 1:length(setstruct)
    setstruct(no).RVV = zeros(1,setstruct(no).TSize);
    setstruct(no).RVEPV = zeros(1,setstruct(no).TSize);
    setstruct(no).RVM = 0;    
    setstruct(no).RVEDV = 0;
    setstruct(no).RVESV = 0;
    setstruct(no).RVSV = 0;
    setstruct(no).RVEF = 0;    
  end
end

%--- Point check
if not(isfield(setstruct,'Point'))
  for no = 1:length(setstruct)
    setstruct(no).Point.X = [];
    setstruct(no).Point.Y = [];
    setstruct(no).Point.T = [];
    setstruct(no).Point.Z = [];    
    setstruct(no).Point.Label = {};
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).Point)
    setstruct(no).Point.X = [];
    setstruct(no).Point.Y = [];
    setstruct(no).Point.T = [];
    setstruct(no).Point.Z = [];    
    setstruct(no).Point.Label = {};
  end
end

%--- 3D Point check
if not(isfield(setstruct,'Point3D'))
  for no = 1:length(setstruct)
    setstruct(no).Point3D.X = [];
    setstruct(no).Point3D.Y = [];
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).Point3D)
    setstruct(no).Point3D.X = [];
    setstruct(no).Point3D.Y = [];
  end
end

%--- IntensityMapping
if not(isfield(setstruct,'IntensityMapping'))
  for no = 1:length(setstruct)
    setstruct(no).IntensityMapping.Brightness = 0.5;
    setstruct(no).IntensityMapping.Contrast = 1;
    setstruct(no).IntensityMapping.Compression = [];
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).IntensityMapping)
    setstruct(no).IntensityMapping.Brightness = 0.5;
    setstruct(no).IntensityMapping.Contrast = 1;
    setstruct(no).IntensityMapping.Compression = [];
  end
end

for no=1:length(setstruct)
  if ~isfield(setstruct(no),'Colormap') %||isempty(setstruct(no).Colormap)
    %remove the empty condition once [] defaults to gray.
%    setstruct(no).Colormap = gray(256);
    setstruct(no).Colormap = [];
  end
end

for loop = 1:length(setstruct)
  if isempty(setstruct(no).PatientInfo)
    setstruct(no).PatientInfo.Name = setstruct(no).FileName;
    setstruct(no).PatientInfo.ID = '';
    setstruct(no).PatientInfo.BirthDate = '';
    setstruct(no).PatientInfo.Age = '';
    setstruct(no).PatientInfo.AcquisitionDate = '';
    setstruct(no).PatientInfo.BSA = 0;
    setstruct(no).PatientInfo.Weight = 0;
    setstruct(no).PatientInfo.Length = 0;
  end
end

for no = 1:length(setstruct)
  if ~isfield(setstruct(no).PatientInfo,'BSA')
    setstruct(no).PatientInfo.BSA = 0;
  end
end

for no = 1:length(setstruct)
  if setstruct(no).PatientInfo.BSA==0
    setstruct(no).PatientInfo.BSA = calcfunctions('calcbsa',...
      setstruct(no).PatientInfo.Weight,...
      setstruct(no).PatientInfo.Length);
  end
end

for no = 1:length(setstruct)
  if isequal(setstruct(no).Flow,0)
    setstruct(no).Flow = [];
  end
  if not(isempty(setstruct(no).Flow))
    if not(isfield(setstruct(no).Flow,'PhaseX'))
      setstruct(no).Flow.PhaseX = [];
      setstruct(no).Flow.PhaseY = [];
    end
    if not(isfield(setstruct(no).Flow,'Angio'))
      setstruct(no).Flow.Angio = [];
      setstruct(no).Flow.VelMag = [];
    end
    if not(isfield(setstruct(no).Flow,'PhaseNo'))
      setstruct(no).Flow.PhaseNo = [];
    end   
    if not(isfield(setstruct(no).Flow,'Result'))
      if isfield(setstruct(no).Flow,'nettotvol')
        try
          setstruct(no).Flow.Result.nettotvol = setstruct(no).Flow.nettotvol;
          setstruct(no).Flow.Result.velmean = setstruct(no).Flow.velmean;
          setstruct(no).Flow.Result.velstd = setstruct(no).Flow.velstd;
          setstruct(no).Flow.Result.velmax = setstruct(no).Flow.velmax;
          setstruct(no).Flow.Result.velmin = setstruct(no).Flow.velmin;
          setstruct(no).Flow.Result.kenergy = setstruct(no).Flow.kenergy;
          setstruct(no).Flow.Result.area = setstruct(no).Flow.area;
          setstruct(no).Flow.Result.netflow = setstruct(no).Flow.netflow;
          setstruct(no).Flow.Result.posflow = setstruct(no).Flow.posflow;
          setstruct(no).Flow.Result.negflow = setstruct(no).Flow.negflow;
          setstruct(no).Flow.Result.diameter = setstruct(no).Flow.diameter;
          setstruct(no).Flow.Result.netforwardvol = setstruct(no).Flow.netforwardvol;
          setstruct(no).Flow.Result.netbackwardvol = setstruct(no).Flow.netbackwardvol;
          setstruct(no).Flow.Result.regfrac = setstruct(no).Flow.regfrac;
          setstruct(no).Flow.Result.sv = setstruct(no).Flow.sv;
        catch
        end
      else
        setstruct(no).Flow.Result = [];
      end
    end        
    if isfield(setstruct(no).Flow,'parameter')
      if not(isfield(setstruct(no).Flow.parameter,'expandoutward'))
        setstruct(no).Flow.parameter.expandoutward = 0; %mm 
        %This will probably be changed when we know that this is a good idea.
      end
    end
  end
end

%--- ImagePosition
if not(isfield(setstruct,'ImagePosition'))
  for no = 1:length(setstruct)
    setstruct(no).ImagePosition = [0 0 0];
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).ImagePosition)
    setstruct(no).ImagePosition = [0 0 0];
  end
end

%--- ImageOrientation
if not(isfield(setstruct,'ImageOrientation'))
  for no = 1:length(setstruct)
    setstruct(no).ImageOrientation = [1 0 0 0 1 0];
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).ImageOrientation)
    setstruct(no).ImageOrientation = [1 0 0 0 1 0];
  end
end

%--- RotatedImagePosition
if not(isfield(setstruct,'RotatedImagePosition'))
  for no = 1:length(setstruct)
    setstruct(no).RotatedImagePosition = setstruct(no).ImagePosition;
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).RotatedImagePosition)
    setstruct(no).RotatedImagePosition = setstruct(no).ImagePosition;
  end
end

%--- RotatedImageOrientation
if not(isfield(setstruct,'RotatedImageOrientation'))
  for no = 1:length(setstruct)
    setstruct(no).RotatedImageOrientation = setstruct(no).ImageOrientation;
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).RotatedImageOrientation)
    setstruct(no).RotatedImageOrientation = setstruct(no).ImageOrientation;
  end
end

%--- EchoTime
if not(isfield(setstruct,'EchoTime'))
  for no = 1:length(setstruct)
    setstruct(no).EchoTime = 0;
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).EchoTime)
    setstruct(no).EchoTime = 0;
  end
end

%--- RepetitionTime
if not(isfield(setstruct,'RepetitionTime'))
  for no = 1:length(setstruct)
    setstruct(no).RepetitionTime = 0;
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).RepetitionTime)
    setstruct(no).RepetitionTime = 0;
  end
end

%--- InversionTime
if not(isfield(setstruct,'InversionTime'))
  for no = 1:length(setstruct)
    setstruct(no).InversionTime = 0;
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).InversionTime)
    setstruct(no).InversionTime = 0;
  end
end

%--- FlipAngle
if not(isfield(setstruct,'FlipAngle'))
  for no = 1:length(setstruct)
    setstruct(no).FlipAngle = 0;
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).FlipAngle)
    setstruct(no).FlipAngle = 0;
  end
end

%--- NumberOfAverages
if not(isfield(setstruct,'NumberOfAverages'))
  for no = 1:length(setstruct)
    setstruct(no).NumberOfAverages = 0;
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).NumberOfAverages)
    setstruct(no).NumberOfAverages = 0;
  end
end

%--- Scanner
if not(isfield(setstruct,'Scanner'))
  for no = 1:length(setstruct)
    setstruct(no).Scanner = '';
  end
end
for no=1:length(setstruct)
  if isempty(setstruct(no).Scanner)
    setstruct(no).Scanner = '';
  end
end

%--- View
if not(isfield(setstruct,'View'))
  for no = 1:length(setstruct)
    setstruct(no).View = [];
  end
%--- ViewPanelsMatrix subfield
elseif not(isempty(setstruct(1).View)) && not(isfield(setstruct(1).View,'ViewPanelsMatrix'))
  panels = setstruct(1).View.ViewPanels;
  setstruct(1).View.ViewPanelsMatrix = cell(1,numel(panels));
  for i = 1:numel(panels)
    if panels(i) ~= 0
      [rows,cols] = calcfunctions('calcrowscols',panels(i),setstruct(panels(i)).ZSize);
    else
      [rows,cols] = deal(1);
    end
    setstruct(1).View.ViewPanelsMatrix{i} = [rows cols];
  end
elseif not(isempty(setstruct(1).View)) && not(isfield(setstruct(1).View,'LVNO'))
  setstruct(1).View.LVNO = [];
  setstruct(1).View.RVNO = [];
  setstruct(1).View.FlowNO = [];
  setstruct(1).View.FlowROI = [];
end
if not(isempty(setstruct(1).View)) && not(isfield(setstruct(1).View,'CurrentTheme'))
  setstruct(1).View.CurrentTheme = 'lv';
end

%--- RotationCenter
if not(isfield(setstruct,'RotationCenter'))
  setstruct(1).RotationCenter = [];
end
for loop = 1:length(setstruct)
  if setstruct(loop).Rotated
    if isempty(setstruct(loop).RotationCenter)
      setstruct(loop).RotationCenter = setstruct(loop).OrgYSize/2-setstruct(loop).YMin;
    end
  end
end

if not(isfield(setstruct,'SequenceName'))
  setstruct(1).SequenceName = '';
end
for loop = 1:length(setstruct)
  if isempty(setstruct(loop).SequenceName)
    setstruct(loop).SequenceName = '';
  end
end

if not(isfield(setstruct,'SeriesDescription'))
  setstruct(1).SeriesDescription = '';
end
for loop = 1:length(setstruct)
  if isempty(setstruct(loop).SeriesDescription)
    setstruct(loop).SeriesDescription = '';
  end
end

if not(isfield(setstruct,'DICOMImageType'))
  setstruct(1).DICOMImageType = '';
end
for loop = 1:length(setstruct)
  if isempty(setstruct(loop).DICOMImageType)
    setstruct(loop).DICOMImageType = '';
  end
end

if not(isfield(setstruct,'Fusion'))
  setstruct(1).Fusion = [];
end

%EH,JT: Fix for backwards compability of .mat files with
%2D or 3D flow.
setstruct(1).ProgramVersion = DATA.ProgramVersion;
for loop=1:length(setstruct)
  setstruct(loop).ProgramVersion = setstruct(1).ProgramVersion;
end

% Backwards compability for old pin structs
for loop=1:length(setstruct)
    if not(isempty(setstruct(loop).EndoPinX)) && all(all(cellfun('isempty',setstruct(loop).EndoPinX)))
        setstruct(loop).EndoPinX = [];
        setstruct(loop).EndoPinY = [];
    end
    if not(isempty(setstruct(loop).EpiPinX)) && all(all(cellfun('isempty',setstruct(loop).EpiPinX)))
        setstruct(loop).EpiPinX = [];
        setstruct(loop).EpiPinY = [];
    end
    if not(isempty(setstruct(loop).RVEndoPinX)) && all(all(cellfun('isempty',setstruct(loop).RVEndoPinX)))
        setstruct(loop).RVEndoPinX = [];
        setstruct(loop).RVEndoPinY = [];
    end
    if not(isempty(setstruct(loop).RVEpiPinX)) && all(all(cellfun('isempty',setstruct(loop).RVEpiPinX)))
        setstruct(loop).RVEpiPinX = [];
        setstruct(loop).RVEpiPinY = [];
    end
end

if not(isfield(setstruct,'AcquisitionTime'))
  for no = 1:length(setstruct)
    setstruct(no).AcquisitionTime = 0;
  end
end

if not(isfield(setstruct,'SeriesNumber'))
  for no = 1:length(setstruct)
    setstruct(no).SeriesNumber = '';
  end
end

if not(isfield(setstruct,'PapillaryIM'))
  for no=1:length(setstruct)
    setstruct(no).PapillaryIM=[];
  end
end

if isfield(setstruct,'StudyUid')
  setstruct = rmfield(setstruct,'StudyUid');
end
%Found case where Endo interp isn't full cell make sure that if there is
%some segmentation the cell is (Tsize,Zsize) i.e pad
%Found case where Endo interp is cell with empty matrix in it performing safety check for this.
for no=1:length(setstruct)
  if iscell(setstruct(no).EndoInterpX)
    if all(cellfun(@isempty,setstruct(no).EndoInterpX))
      setstruct(no).EndoInterpX=[];
      setstruct(no).EndoInterpY=[];
    else
      correcttsize = (size(setstruct(no).EndoInterpX,1) == setstruct(no).TSize);
      correctzsize = (size(setstruct(no).EndoInterpX,2) == setstruct(no).ZSize);
      if not(correcttsize) && correctzsize
        try
          pad = cell(setstruct(no).TSize-size(setstruct(no).EndoInterpX,1),setstruct(no).ZSize);
          setstruct(no).EndoInterpX=[setstruct(no).EndoInterpX;pad];
          setstruct(no).EndoInterpY=[setstruct(no).EndoInterpY;pad];
        catch
          disp('Issue in reading RV interpolation points. Deleting the points.');
          setstruct(no).EndoInterpX = cell(setstruct(no).TSize,setstruct(no).ZSize);
          setstruct(no).EndoInterpY = setstruct(no).EndoInterpX;
        end
      elseif not(correctzsize)
        disp('Issue in reading RV interpolation points. Deleting the points.');
        setstruct(no).EndoInterpX = cell(setstruct(no).TSize,setstruct(no).ZSize);
        setstruct(no).EndoInterpY = setstruct(no).EndoInterpX;
      end
    end
  end
  
  if iscell(setstruct(no).EpiInterpX)
    if all(cellfun(@isempty,setstruct(no).EpiInterpX))
      setstruct(no).EpiInterpX=[];
      setstruct(no).EpiInterpY=[];
    else
      correcttsize = (size(setstruct(no).EpiInterpX,1) == setstruct(no).TSize);
      correctzsize = (size(setstruct(no).EpiInterpX,2) == setstruct(no).ZSize);
      if not(correcttsize) && correctzsize
        try
          pad = cell(setstruct(no).TSize-size(setstruct(no).EpiInterpX,1),setstruct(no).ZSize);
          setstruct(no).EpiInterpX=[setstruct(no).EpiInterpX;pad];
          setstruct(no).EpiInterpY=[setstruct(no).EpiInterpY;pad];
        catch
          disp('Issue in reading RV interpolation points. Deleting the points.');
          setstruct(no).EpiInterpX = cell(setstruct(no).TSize,setstruct(no).ZSize);
          setstruct(no).EpiInterpY = setstruct(no).EpiInterpX;
        end
      elseif not(correctzsize)
        disp('Issue in reading RV interpolation points. Deleting the points.');
        setstruct(no).EpiInterpX = cell(setstruct(no).TSize,setstruct(no).ZSize);
        setstruct(no).EpiInterpY = setstruct(no).EpiInterpX;
      end
    end
  end
  
  if iscell(setstruct(no).RVEndoInterpX)
    if all(cellfun(@isempty,setstruct(no).RVEndoInterpX))
      setstruct(no).RVEndoInterpX=[];
      setstruct(no).RVEndoInterpY=[];
    else
      correcttsize = (size(setstruct(no).RVEndoInterpX,1) == setstruct(no).TSize);
      correctzsize = (size(setstruct(no).RVEndoInterpX,2) == setstruct(no).ZSize);
      if not(correcttsize) && correctzsize
        try
          pad = cell(setstruct(no).TSize-size(setstruct(no).RVEndoInterpX,1),setstruct(no).ZSize);
          setstruct(no).RVEndoInterpX=[setstruct(no).RVEndoInterpX;pad];
          setstruct(no).RVEndoInterpY=[setstruct(no).RVEndoInterpY;pad];
        catch
          disp('Issue in reading RV interpolation points. Deleting the points.');
          setstruct(no).RVEndoInterpX = cell(setstruct(no).TSize,setstruct(no).ZSize);
          setstruct(no).RVEndoInterpY = setstruct(no).RVEndoInterpX;
        end
      elseif not(correctzsize)
        disp('Issue in reading RV interpolation points. Deleting the points.');
        setstruct(no).RVEndoInterpX = cell(setstruct(no).TSize,setstruct(no).ZSize);
        setstruct(no).RVEndoInterpY = setstruct(no).RVEndoInterpX;
      end
    end
  end
  
  if iscell(setstruct(no).RVEpiInterpX)
    if all(cellfun(@isempty,setstruct(no).RVEpiInterpX))
      setstruct(no).RVEpiInterpX=[];
      setstruct(no).RVEpiInterpY=[];
    else
      correcttsize = (size(setstruct(no).RVEpiInterpX,1) == setstruct(no).TSize);
      correctzsize = (size(setstruct(no).RVEpiInterpX,2) == setstruct(no).ZSize);
      if not(correcttsize) && correctzsize
        try
          pad = cell(setstruct(no).TSize-size(setstruct(no).RVEpiInterpX,1),setstruct(no).ZSize);
          setstruct(no).RVEpiInterpX=[setstruct(no).RVEpiInterpX;pad];
          setstruct(no).RVEpiInterpY=[setstruct(no).RVEpiInterpY;pad];
        catch
          disp('Issue in reading RV interpolation points. Deleting the points.');
          setstruct(no).RVEpiInterpX = cell(setstruct(no).TSize,setstruct(no).ZSize);
          setstruct(no).RVEpiInterpY = setstruct(no).RVEpiInterpX;
        end
      elseif not(correctzsize)
        disp('Issue in reading RV interpolation points. Deleting the points.');
        setstruct(no).RVEpiInterpX = cell(setstruct(no).TSize,setstruct(no).ZSize);
        setstruct(no).RVEpiInterpY = setstruct(no).RVEpiInterpX;
      end
    end
  end
end

%Check if valid interpx
for no = 1:length(setstruct)
  if ~isempty(setstruct(no).EndoInterpX)
    for tloop=1:setstruct(no).TSize
      for sloop=1:setstruct(no).ZSize
        if size(setstruct(no).EndoInterpX{tloop,sloop},2)>size(setstruct(no).EndoInterpX{tloop,sloop},1)
          setstruct(no).EndoInterpX{tloop,sloop} = setstruct(no).EndoInterpX{tloop,sloop}(:);
          setstruct(no).EndoInterpY{tloop,sloop} = setstruct(no).EndoInterpY{tloop,sloop}(:);          
          disp('Warning interpolation point size error. Now fixed.');          
        end
      end
    end
  end
end

for no = 1:length(setstruct)
  if ~isempty(setstruct(no).EpiInterpX)
    for tloop=1:setstruct(no).TSize
      for sloop=1:setstruct(no).ZSize
        if size(setstruct(no).EpiInterpX{tloop,sloop},2)>size(setstruct(no).EpiInterpX{tloop,sloop},1)
          setstruct(no).EpiInterpX{tloop,sloop} = setstruct(no).EpiInterpX{tloop,sloop}(:);
          setstruct(no).EpiInterpY{tloop,sloop} = setstruct(no).EpiInterpY{tloop,sloop}(:);          
          disp('Warning interpolation point size error. Now fixed.');          
        end
      end
    end
  end
end

for no = 1:length(setstruct)
  if ~isempty(setstruct(no).RVEndoInterpX)
    for tloop=1:setstruct(no).TSize
      for sloop=1:setstruct(no).ZSize
        if size(setstruct(no).RVEndoInterpX{tloop,sloop},2)>size(setstruct(no).RVEndoInterpX{tloop,sloop},1)
          setstruct(no).RVEndoInterpX{tloop,sloop} = setstruct(no).RVEndoInterpX{tloop,sloop}(:);
          setstruct(no).RVEndoInterpY{tloop,sloop} = setstruct(no).RVEndoInterpY{tloop,sloop}(:);          
          disp('Warning interpolation point size error. Now fixed.');          
        end
      end
    end
  end
end

for no = 1:length(setstruct)
  if ~isempty(setstruct(no).RVEpiInterpX)
    for tloop=1:setstruct(no).TSize
      for sloop=1:setstruct(no).ZSize
        if size(setstruct(no).RVEpiInterpX{tloop,sloop},2)>size(setstruct(no).RVEpiInterpX{tloop,sloop},1)
          setstruct(no).RVEpiInterpX{tloop,sloop} = setstruct(no).RVEpiInterpX{tloop,sloop}(:);
          setstruct(no).RVEpiInterpY{tloop,sloop} = setstruct(no).RVEpiInterpY{tloop,sloop}(:);          
          disp('Warning interpolation point size error. Now fixed.');          
        end
      end
    end
  end
end

%set image view plane based on old image type definition
if ~isfield(setstruct,'ImageViewPlane')
  for no = 1:length(setstruct) 
    oldImageType = setstruct(no).ImageType;
    [type,viewplane] = segment('imagedescription');
    settype = 0;
    setviewplane = 0;
    if isfield(setstruct,'ImageType')
      for typeloop = 1:length(type)
        if findstr(lower(oldImageType),lower(type{typeloop})) %#ok<FSTR>
          setstruct(no).ImageType = type{typeloop};
          settype = 1;
        end
      end
      %First use Einars auto detect algoritm to see if it finds the viewplane 
      outlabel=openfile('autodetectviewplane',setstruct(no).ImageOrientation);
      if isempty(outlabel)
        for viewplaneloop = 1:length(viewplane)
          if findstr(lower(oldImageType),lower(viewplane{viewplaneloop})) %#ok<FSTR>
            setstruct(no).ImageViewPlane = viewplane{viewplaneloop};
            setviewplane = 1;
          end
        end
      else
        setstruct(no).ImageViewPlane = outlabel;
        setviewplane = 1;
      end
    end
    if settype == 0  %couldn't find image type
      setstruct(no).ImageType = 'General';
    end
    if setviewplane == 0  %couldn't find image view plane
      setstruct(no).ImageViewPlane = 'Unspecified';
    end
  end
end

if not(isfield(setstruct,'MaR'))||isfield(setstruct,'MaRIM')% if field MaRIM exist an old way of implementing MaR has been used
  for no=1:length(setstruct)
    setstruct(no).MaR=[];
    if isfield(setstruct,'MaRIM')
      if not(isempty(setstruct(no).MaRIM))
        mar('initdefault',no);
        mar('createmyocardmask',no);
        setstruct(no).MaR.Manual=setstruct(no).MaRIM;
        mar('update',no,true);%true for silent
      end
    end
  end
end
if isfield(setstruct,'MaRIM')
  setstruct=rmfield(setstruct,'MaRIM');
end
for no=1:length(setstruct)
  if not(isempty(setstruct(no).MaR)) && not(isfield(setstruct(no).MaR, 'MPS'))
    setstruct(no).MaR.MPS.DefectThreshold=55;
    setstruct(no).MaR.MPS.DefectPercentile=100;
    setstruct(no).MaR.MPS.DefectMinVolume=10;
    setstruct(no).MaR.MPS.DefectRegion='roi';
    setstruct(no).MaR.MPS.UsingModel=1;
    setstruct(no).MaR.MPS.ModelType='mercator';
    setstruct(no).MaR.MPS.Dominant=[];
    setstruct(no).MaR.MPS.TPD=[];
    setstruct(no).MaR.MPS.TPDLAD=[];
    setstruct(no).MaR.MPS.TPDLCx=[];
    setstruct(no).MaR.MPS.TPDRCA=[];
    setstruct(no).MaR.MPS.SeverityIndex=[];
    setstruct(no).MaR.MPS.NadirIndex=[];
  end
end

for no=1:length(setstruct)
  if isfield(setstruct(no).MaR,'MR')
    if not(isfield(setstruct(no).MaR.MR,'Beta'))
      setstruct(no).MaR.MR.Beta = 0.2;
    end
    if not(isfield(setstruct(no).MaR.MR,'Radius'))
      setstruct(no).MaR.MR.Radius = 2*mean([setstruct(no).ResolutionX setstruct(no).ResolutionY]);
    end
    if not(isfield(setstruct(no).MaR.MR,'Artery'))
      setstruct(no).MaR.MR.Artery=[];
    end
    if not(isfield(setstruct(no).MaR.MR,'Rotation'))
      setstruct(no).MaR.MR.Rotation=[];
    end
  end
end

%Perfusion
if not(isfield(setstruct,'Perfusion'))
  for no = 1:length(setstruct)
    setstruct(no).Perfusion = [];
%     perfusion('initdefault',no);
%     perfusion('createmyocardmask',no);
  end
end

%Perfusion scoring
if not(isfield(setstruct,'PerfusionScoring'))
  for no = 1:length(setstruct)
    setstruct(no).PerfusionScoring = [];
  end
end

%T2*
if not(isfield(setstruct,'T2'))
  for no = 1:length(setstruct)
    setstruct(no).T2 = [];
%     perfusion('initdefault',no);
%     perfusion('createmyocardmask',no);
  end
end

%check to remove zeros that could have been inserted instead of NaNs when
%resampling in versions prior to r1021
if not(isempty(setstruct(no).EndoX))
  removezeros = (setstruct(no).EndoX==0) & (setstruct(no).EndoY==0);
  setstruct(no).EndoX(removezeros)=NaN;
  setstruct(no).EndoY(removezeros)=NaN;
end
if not(isempty(setstruct(no).EpiX))
  removezeros = (setstruct(no).EpiX==0) & (setstruct(no).EpiY==0);
  setstruct(no).EpiX(removezeros)=NaN;
  setstruct(no).EpiY(removezeros)=NaN;
end
if not(isempty(setstruct(no).RVEndoX))
  removezeros = (setstruct(no).RVEndoX==0) & (setstruct(no).RVEndoY==0);
  setstruct(no).RVEndoX(removezeros)=NaN;
  setstruct(no).RVEndoY(removezeros)=NaN;
end
if not(isempty(setstruct(no).RVEpiX))
  removezeros = (setstruct(no).RVEpiX==0) & (setstruct(no).RVEpiY==0);
  setstruct(no).RVEpiX(removezeros)=NaN;
  setstruct(no).RVEpiY(removezeros)=NaN;
end

%check if field LongName used for measure report tables is present
for no = 1:length(setstruct)
  if not(isempty(setstruct(no).Measure)) && not(isfield(setstruct(no).Measure,'LongName'))
    for sloop = 1:length(setstruct(no).Measure)
      setstruct(no).Measure(sloop).LongName = setstruct(no).Measure(sloop).Name;
    end
  end
end

%check if linking structure is present
for no = 1:length(setstruct)
  didnotlink = false;
  if not(isfield(setstruct(no),'Linked'))
    didnotlink = true;
    setstruct(no).Children = [];
    setstruct(no).Parent = [];
    setstruct(no).Linked = no;
  elseif isempty(setstruct(no).Linked)
    didnotlink = true;
    setstruct(no).Linked = no;
  end
  if ~isempty(setstruct(no).Flow) && didnotlink
    if no == setstruct(no).Flow.MagnitudeNo
      nos = [...
        setstruct(no).Flow.PhaseNo ...
        setstruct(no).Flow.PhaseX ...
        setstruct(no).Flow.PhaseY ...
        setstruct(no).Flow.Angio ...
        setstruct(no).Flow.VelMag ...
        ];
      setstruct(no).Children = unique([setstruct(no).Children nos]);
    else
      setstruct(no).Parent = setstruct(no).Flow.MagnitudeNo;
    end
    nos = [...
      setstruct(no).Flow.MagnitudeNo ...
      setstruct(no).Flow.PhaseNo ...
      setstruct(no).Flow.PhaseX ...
      setstruct(no).Flow.PhaseY ...
      setstruct(no).Flow.Angio ...
      setstruct(no).Flow.VelMag ...
      ];
    newlinkies = unique([setstruct(no).Linked nos]);
    if ~isempty(setstruct(newlinkies(1)).Parent) %make sure parent is at index = 1
      newlinkies = [setstruct(newlinkies(1)).Parent newlinkies(newlinkies ~= setstruct(newlinkies(1)).Parent)];
    end
    setstruct(no).Linked = newlinkies;
  end
end

%add field for overlays
if ~isfield(setstruct(no),'Overlay')
  setstruct(no).Overlay = [];
end

%--- AccessionNumber
if not(isfield(setstruct,'AccessionNumber'))
  for no = 1:length(setstruct)
    setstruct(no).AccessionNumber = '';
  end
end

%--- StudyID
if not(isfield(setstruct,'StudyID'))
  for no = 1:length(setstruct)
    setstruct(no).StudyID = '';
  end
end

%--- StudyUID
if not(isfield(setstruct,'StudyUID'))
  for no = 1:length(setstruct)
    setstruct(no).StudyUID = '';
  end
end

%--- Intersection
if not(isfield(setstruct,'Intersection'))
  for no = 1:length(setstruct)
    setstruct(no).Intersection = [];
  end
end

%Bug check for EST and EDT if its zero
for no = 1:length(setstruct)
  if setstruct(no).EST==0
    setstruct(no).EST=1;
  end
  if setstruct(no).EDT==0
    setstruct(no).EDT=1;  
  end
end

%--- TimeVector
for no = 1:length(setstruct) 
  if not(isfield(setstruct,'TimeVector')) || isempty(setstruct(no).TimeVector) || all(setstruct(no).TimeVector==0) 
    if isnan(setstruct(no).TIncr)
      setstruct(no).TIncr=1/(SET(no).TSize-1); %Normalised time
    end
    tvec = 0:setstruct(no).TIncr:(setstruct(no).TSize-1) * setstruct(no).TIncr;
    if isempty(tvec)
      tvec = 0;
    end
    setstruct(no).TimeVector = tvec;
  end
end

for no = 1:length(setstruct) 
  if (isempty(setstruct(no).CurrentSlice))
    setstruct(no).CurrentSlice = 1;
    setstruct(no).EndSlice = 1;
    setstruct(no).StartSlice = 1;
  end
end

%--- Comment
if not(isfield(setstruct,'Comment'))
  setstruct(1).Comment = [];
end

%--- AtrialScar
if not(isfield(setstruct,'AtrialScar'))
  setstruct(1).AtrialScar = [];
end

%--- PatientInfo -> Institution
for no = 1:length(setstruct)
  if ~isfield(setstruct(no).PatientInfo,'Institution')
    setstruct(no).PatientInfo.Institution = '';
  end
end

%--- PatientInfo -> Remove obsolete HeartRate field (later)
% for no = 1:length(setstruct)
%   if isfield(setstruct(no).PatientInfo,'HeartRate')
%     hr = setstruct(no).PatientInfo.HeartRate; %#ok<NASGU>
%     setstruct(no).PatientInfo = rmfield(setstruct(no).PatientInfo,'HeartRate');
%   end;
% end

%--- LastUserInfo 
for no = 1:length(setstruct)
  if ~isfield(setstruct(no),'LastUserInfo')
    setstruct(no).LastUserInfo = '';
  end
end


%--- Extra views for SPECT
if ~isfield(setstruct,'HLA')
  [setstruct.SAX3] = deal([]);
  [setstruct.HLA] = deal([]);
  [setstruct.VLA] = deal([]);
end
if ~isfield(setstruct,'GLA')
  [setstruct.GLA] = deal([]);
end

%--- Mmode struct
if isfield(setstruct,'Mmodex')
  [setstruct.Mmode] = deal([]);
  for no = 1:numel(setstruct)
    mmode = struct( ...
      'X', setstruct(no).Mmodex, ...
      'Y', setstruct(no).Mmodey, ...
      'Lx', setstruct(no).Mmodelx, ...
      'Ly', setstruct(no).Mmodely, ...
      'M1', setstruct(no).Mmodem1, ...
      'M2', setstruct(no).Mmodem2, ...
      'T1', 1, ...
      'T2', setstruct(no).TSize);
    setstruct(no).Mmode = mmode;
  end
  setstruct = rmfield(setstruct,{'Mmodex','Mmodey','Mmodelx','Mmodely', ...
    'Mmodem1','Mmodem2'});
end

if ~isfield(setstruct,'CT')
  [setstruct.CT] = deal([]);
end

if ~isfield(setstruct,'RVPFR')
  [setstruct.RVPFR] = deal(0);
  [setstruct.RVPER] = deal(0);
  [setstruct.RVPFRT] = deal(1);
  [setstruct.RVPERT] = deal(1);
end

if ~isfield(setstruct,'PapillaryThreshold')
  [setstruct.PapillaryThreshold] = deal(0);
end

if ~isfield(setstruct,'ECV')
  [setstruct.ECV] = deal([]);
end

if ~isfield(setstruct,'Developer')
  [setstruct.Developer] = deal([]);
end

if ~isfield(setstruct,'Line3D')
  %Used for 3D lines in Segment 3DP
  for loop = 1:length(setstruct)
    setstruct(loop).Line3D = [];
    setstruct(loop).Line3D.X = {};
    setstruct(loop).Line3D.Y = {};
    setstruct(loop).Line3D.Z = {};
    setstruct(loop).Line3D.Points = {};
    setstruct(loop).Line3D.Closed = [];
    setstruct(loop).Line3D.Sigma = [];
  end
end

if ~isfield(setstruct,'PVLoop')
  %Used for the PV-loop module
  setstruct(1).PVLoop = []; %empty is enough to initialise as
end

%%%%%%%%%%%%%
%%% Check validity of fields %%%
%%%%%%%%%%%%%

%Check if valid interpx
for no = 1:length(setstruct)
    if ~isempty(setstruct(no).EndoInterpX) && not(isempty(setstruct(no).EndoX))
        for tloop=1:setstruct(no).TSize
            for sloop=1:setstruct(no).ZSize
                if not(isempty(setstruct(no).EndoInterpX{tloop,sloop})) && not(isnan(setstruct(no).EndoX(1,tloop,sloop)))
                    ipx=setstruct(no).EndoInterpX{tloop,sloop};
                    ipy=setstruct(no).EndoInterpY{tloop,sloop};
                    contx = setstruct(no).EndoX(:,tloop,sloop);
                    conty = setstruct(no).EndoY(:,tloop,sloop);
                    
                    ipxrep=repmat(ipx',[length(contx) 1]);
                    contxrep=repmat(contx,[1 length(ipx)]);
                    ipyrep=repmat(ipy',[length(conty) 1]);
                    contyrep=repmat(conty,[1 length(ipy)]);
                    pindist2cont = (ipxrep-contxrep).^2+(ipyrep-contyrep).^2;
                    [~,mindistindex] = min(pindist2cont);
                    [~,sortindex] = sort(mindistindex);
                    ipx=ipx(sortindex);
                    ipy=ipy(sortindex);
                    
                    setstruct(no).EndoInterpX{tloop,sloop}=ipx;
                    setstruct(no).EndoInterpY{tloop,sloop}=ipy;
                end
            end
        end
    end
end

for no = 1:length(setstruct)
    if ~isempty(setstruct(no).EpiInterpX)  && not(isempty(setstruct(no).EpiX))
        for tloop=1:setstruct(no).TSize
            for sloop=1:setstruct(no).ZSize
                if not(isempty(setstruct(no).EpiInterpX{tloop,sloop})) && not(isnan(setstruct(no).EpiX(1,tloop,sloop)))
                    ipx=setstruct(no).EpiInterpX{tloop,sloop};
                    ipy=setstruct(no).EpiInterpY{tloop,sloop};
                    contx = setstruct(no).EpiX(:,tloop,sloop);
                    conty = setstruct(no).EpiY(:,tloop,sloop);
                    
                    ipxrep=repmat(ipx',[length(contx) 1]);
                    contxrep=repmat(contx,[1 length(ipx)]);
                    ipyrep=repmat(ipy',[length(conty) 1]);
                    contyrep=repmat(conty,[1 length(ipy)]);
                    pindist2cont = (ipxrep-contxrep).^2+(ipyrep-contyrep).^2;
                    [~,mindistindex] = min(pindist2cont);
                    [~,sortindex] = sort(mindistindex);
                    ipx=ipx(sortindex);
                    ipy=ipy(sortindex);
                    
                    setstruct(no).EpiInterpX{tloop,sloop}=ipx;
                    setstruct(no).EpiInterpY{tloop,sloop}=ipy;
                end
            end
        end
    end
end

for no = 1:length(setstruct)
    if ~isempty(setstruct(no).RVEndoInterpX)  && not(isempty(setstruct(no).RVEndoX))
        for tloop=1:setstruct(no).TSize
            for sloop=1:setstruct(no).ZSize
                if not(isempty(setstruct(no).RVEndoInterpX{tloop,sloop})) && not(isnan(setstruct(no).RVEndoX(1,tloop,sloop)))
                    ipx=setstruct(no).RVEndoInterpX{tloop,sloop};
                    ipy=setstruct(no).RVEndoInterpY{tloop,sloop};
                    contx = setstruct(no).RVEndoX(:,tloop,sloop);
                    conty = setstruct(no).RVEndoY(:,tloop,sloop);
                    
                    ipxrep=repmat(ipx',[length(contx) 1]);
                    contxrep=repmat(contx,[1 length(ipx)]);
                    ipyrep=repmat(ipy',[length(conty) 1]);
                    contyrep=repmat(conty,[1 length(ipy)]);
                    pindist2cont = (ipxrep-contxrep).^2+(ipyrep-contyrep).^2;
                    [~,mindistindex] = min(pindist2cont);
                    [~,sortindex] = sort(mindistindex);
                    ipx=ipx(sortindex);
                    ipy=ipy(sortindex);
                    
                    setstruct(no).RVEndoInterpX{tloop,sloop}=ipx;
                    setstruct(no).RVEndoInterpY{tloop,sloop}=ipy;
                end
            end
        end
    end
end

for no = 1:length(setstruct)
    if ~isempty(setstruct(no).RVEpiInterpX)  && not(isempty(setstruct(no).RVEpiX))
        for tloop=1:setstruct(no).TSize
            for sloop=1:setstruct(no).ZSize
                if not(isempty(setstruct(no).RVEpiInterpX{tloop,sloop})) && not(isnan(setstruct(no).RVEpiX(1,tloop,sloop)))
                    ipx=setstruct(no).RVEpiInterpX{tloop,sloop};
                    ipy=setstruct(no).RVEpiInterpY{tloop,sloop};
                    contx = setstruct(no).RVEpiX(:,tloop,sloop);
                    conty = setstruct(no).RVEpiY(:,tloop,sloop);
                    
                    ipxrep=repmat(ipx',[length(contx) 1]);
                    contxrep=repmat(contx,[1 length(ipx)]);
                    ipyrep=repmat(ipy',[length(conty) 1]);
                    contyrep=repmat(conty,[1 length(ipy)]);
                    pindist2cont = (ipxrep-contxrep).^2+(ipyrep-contyrep).^2;
                    [~,mindistindex] =min(pindist2cont);
                    [~,sortindex] =sort(mindistindex);
                    ipx=ipx(sortindex);
                    ipy=ipy(sortindex);
                    
                    setstruct(no).RVEpiInterpX{tloop,sloop}=ipx;
                    setstruct(no).RVEpiInterpY{tloop,sloop}=ipy;
                end
            end
        end
    end
end

%%%%%%%%%%%%%
%%% Obsolete fields
%%%%%%%%%%%%%

%if isfield(setstruct,'Colormap')
%  setstruct = rmfield(setstruct,'Colormap');
%end;

if isfield(setstruct,'CurrentTool')
  setstruct = rmfield(setstruct,'CurrentTool');
end

if isfield(setstruct,'Model')
  setstruct = rmfield(setstruct,'Model');
end

if isfield(setstruct,'ViewMode')
  setstruct = rmfield(setstruct,'ViewMode');
end

if isfield(setstruct,'rows')
  setstruct = rmfield(setstruct,{'rows','cols'});
end

if isfield(setstruct,'RoiArea')
  setstruct = rmfield(setstruct,'RoiArea');
end

if nargin < 1
  SET = setstruct;
end

%-------------------------------------------------------------------
function fieldvalues = rebuildflowresults(fieldvalues,numtimeframes)
%-------------------------------------------------------------------
% rewrites flow results so that their length correspond to the number of
% time frames and existing values are placed in the corrsponding places and
% arraz indexes without values are set to NaN

[~,ind] = find(fieldvalues);
newvalues = nan(1,numtimeframes);
newvalues(ind) = fieldvalues(ind);
fieldvalues = newvalues;
