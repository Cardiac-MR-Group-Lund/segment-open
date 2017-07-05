function varargout = plugin_summarize(fcn,varargin)
%--------------------------------------------------
%Example code to do batch processing of multiple files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This section is standard template %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0 
  myfailed('Requires more than one input argument.');
  return;
end;

switch fcn
	case 'getname'
		varargout = cell(1,1);
    
    %%%% Start your change heres %%%%
		varargout{1} = 'Summarize'; %What appears in the menu in Segment.

		%Segment with versions >1.636 calls with two input arguments where
		%the second input argument is the handle to the menu item.

		%Register submenus
		uimenu(varargin{1},'Label','Multiple mat files','Callback','plugin_summarize(''multiple'')');
		uimenu(varargin{1},'Label','One file','Callback','plugin_summarize(''one'')');

    %%%% End of changes section %%%%    
    
    %Set the main "summarize" menu to not perform a callback.
		set(varargin{1},'Callback','');
  case 'getdependencies'
    %Here: List all depending files. This is required if your plugin should
    %be possible to compile to the stand-alone version of Segment.
    varargout = cell(1,4);

    %M-files
    varargout{1} = {};

    %Fig-files
    varargout{2} = {};

    %Mat-files
    varargout{3} = {};

    %Mex-files
    varargout{4} = {};
    
  otherwise
    macro_helper(fcn,varargin{:});
		[varargout{1:nargout}] = feval(fcn,varargin{:}); % FEVAL switchyard
end;

%------------------------------
function [outdata,ind] = header
%------------------------------
outdata{1, 1} = 'FileName';
outdata{1, 2} = 'PatientName';
outdata{1, 3} = 'PatientID';
outdata{1, 4} = 'AcquisitionDate';
outdata{1, 5} = 'Age';
outdata{1, 6} = 'Length[cm]';
outdata{1, 7} = 'Weight[kg]';
outdata{1, 8} = 'Sex';
outdata{1, 9} = 'BSA[m2]';
outdata{1,10} = 'HeartRate[bpm]';
outdata{1,11} = 'R-R[ms]';
outdata{1,12} = 'LVM[ml]';
outdata{1,13} = 'LVM[g]';
outdata{1,14} = 'LVMI[g/m2]';
outdata{1,15} = 'EDV[ml]';
outdata{1,16} = 'EDVI[ml/m2]';
outdata{1,17} = 'ESV[ml]';
outdata{1,18} = 'ESVI[ml/m2]';
outdata{1,19} = 'SV[ml]';
outdata{1,20} = 'SVI[ml/m2]';
outdata{1,21} = 'EF[%]';
outdata{1,22} = 'CO[l/min]';
outdata{1,23} = 'CI[l/min]';
outdata{1,24} = 'PFR[ml/s]';
outdata{1,25} = 'PER[ml/s]';
outdata{1,26} = 'RVM[ml]';
outdata{1,27} = 'RVM[g]';
outdata{1,28} = 'RVMI[g/m2]';
outdata{1,29} = 'RVEDV[ml]';
outdata{1,30} = 'RVEDVI[ml/m2]';
outdata{1,31} = 'RVESV[ml]';
outdata{1,32} = 'RVESVI[ml/m2]';
outdata{1,33} = 'RVSV[ml]';
outdata{1,34} = 'RVSVI[ml/m2]';
outdata{1,35} = 'RVEF[%]';
outdata{1,36} = 'LVM_DE[ml]';
outdata{1,37} = 'Scar[%]';
outdata{1,38} = 'Scar[ml](fromLVMDE)';
outdata{1,39} = 'TotExt[%]';
outdata{1,40} = 'MeanTransmurality[%]';
outdata{1,41} = 'MaxTransmurality[%]';
outdata{1,42} = 'TotExt[%](Weighted)';
outdata{1,43} = 'MeanTransmurality[%](Weighted)';
outdata{1,44} = 'MaxTransmurality[%](Weighted)';
outdata{1,45} = 'MO region[%]';
outdata{1,46} = 'MO[%]';

ind = [9 14 16 18 20 23 28 30 32 34];

%----------------
function multiple %#ok<DEFNU>
%----------------
%Example function to summarize many .mat files and export
%to clipboard.

global DATA SET NO

%Select path
pathname = DATA.Pref.datapath;
pathname = myuigetdir(pathname,'Select a folder with .mat files');
if isequal(pathname,0)
  myfailed('Aborted.');
  return;
end;

%Find files to process
files2load = dir([pathname filesep '*.mat']);
numfiles = length(files2load);

if numfiles==0
  myfailed('Found no files to summarize.');
  return;
end;

includenormalized = yesno('Do you want to include BSA normalized values?');

%Create output matrix
outdata = cell(numfiles+1,37); %+1 since header, 37 since header size
[temp,indforbsa] = header;

%--- Write header
for loop=1:length(temp)
  outdata{1,loop} = temp{loop};
end;

%Loop over all files
h = mywaitbarstart(numfiles,'Please wait, loading and summarizing files.',1);
for fileloop=1:numfiles

  %--- Load file
  DATA.Silent = true; %Turn on "silent" mode to avoid to much update on screen when loading etc.
  
  disp(sprintf('Loading %s.',files2load(fileloop).name)); %#ok<DSPS>
  
  %Set filename
  outdata{fileloop+1,1} = files2load(fileloop).name;
  
  SET = []; %Make sure a fresh start
  
  %---- try 
  try
    SET=[];
    %Load
    load([pathname filesep files2load(fileloop).name],'-mat');
    
    %Assign
    SET = setstruct;
    clear setstruct;
    
    %Call to intialize all variables correcly after loaded data.
    openfile('setupstacksfrommat',1);
    segment('renderstacksfrommat');
    
    %Call one to get info in one line.
    temp = one(false); %no header;
    for loop=2:size(temp,2);
      outdata{fileloop+1,loop} = temp{1,loop};
    end;
  
    %Do special with flow and multiple
    
    %--- Find correct image stacks
    [~,~,flowno] = findfunctions('findno');
    
    roipos = 41;
    if ~isempty(flowno)
      for loop = 1:length(flowno)
        NO = flowno(loop);
        
        if SET(flowno(loop)).RoiN>0
          %Calculate flow, call to get total flow.
          reportflow;
          tots = reportflow('gettotal');
          reportflow('close_Callback');
          
          for rloop=1:length(tots)
            outdata{1,roipos} = 'ROI_name';
            outdata{1,roipos+1} = 'tot-flow';
            outdata{fileloop+1,roipos} = SET(NO).Roi(rloop).Name;
            outdata{fileloop+1,roipos+1} = tots(rloop);
            roipos = roipos+2;
          end;
        end;
      end;
    end;
    
    measurementpos = roipos;
    for sloop=1:length(SET)
      for mloop=1:length(SET(sloop).Measure)
        outdata{1,measurementpos} = 'Measure_name';
        outdata{1,measurementpos+1} = '[mm]';
        outdata{fileloop+1,measurementpos} = SET(sloop).Measure(mloop).Name;
        outdata{fileloop+1,measurementpos+1} = SET(sloop).Measure(mloop).Length;
        measurementpos = measurementpos+2;
      end;
    end;
    
   catch %#ok<CTCH>
     %--- Some thing went wrong
     outdata{fileloop+1,2} = 'FAILED.';
   end
    h = mywaitbarupdate(h);
end; %loop over files
mywaitbarclose(h);

%If not wanted remove non normalized values.
ind = true(1,size(outdata,2));
if ~includenormalized
  ind(indforbsa) = false;
end;
outdata = outdata(:,ind);

%--- Output to a string
segment('cell2clipboard',outdata);

%Make sure starting with something fresh.
segment('filecloseall_Callback',true);

%Stop the silent mode.
DATA.Silent = false;

%---------------------------------
function varargout = one(doheader)
%---------------------------------
%Example code to just summarize the currently loaded image stack
%If called with an extra argument then also add the header.

global SET NO

[no,scarno,flowno,varargout] = findfunctions('findno');

if isnan(no)
  no = NO;
end;

outdata = cell(1,20);

if nargin==0
  doheader = true;
end;

if doheader
  %Write Header
  outdata = header;
  row = 2;
else
  row = 1;
end;

if nargout>0
  varargout = cell(1,nargout);
end;

%Write data
outdata{row, 1} = SET(no).FileName;
outdata{row, 2} = SET(no).PatientInfo.Name;
outdata{row, 3} = SET(no).PatientInfo.ID;
outdata{row, 4} = SET(no).PatientInfo.AcquisitionDate;
outdata{row, 5} = SET(no).PatientInfo.Age;
outdata{row, 6} = SET(no).PatientInfo.Length;
outdata{row, 7} = SET(no).PatientInfo.Weight;
outdata{row, 8} = SET(no).PatientInfo.Sex;
outdata{row, 9} = SET(no).PatientInfo.BSA;
outdata{row,10} = SET(no).HeartRate;
outdata{row,11} = 1000*SET(no).TSize*SET(no).TIncr;

lvm = 0.5*(SET(no).LVM(SET(no).EDT)+SET(no).LVM(SET(no).EST));
try
  rvm = 0.5*(SET(no).RVM(SET(no).EDT)+SET(no).RVM(SET(no).EST));
catch %#ok<CTCH>
  rvm = NaN;
end;

if SET(no).PatientInfo.BSA==0
  bsa_1 = NaN;
else
  bsa_1 = 1/SET(no).BSA;
end;

outdata{row,12} = lvm;
outdata{row,13} = lvm*1.05;
outdata{row,14} = bsa_1*lvm*1.05;
outdata{row,15} = SET(no).EDV;
outdata{row,16} = bsa_1*SET(no).EDV;
outdata{row,17} = SET(no).ESV;
outdata{row,18} = bsa_1*SET(no).ESV;
outdata{row,19} = SET(no).SV;
outdata{row,20} = bsa_1*SET(no).SV;
outdata{row,21} = 100*SET(no).EF;
outdata{row,22}= SET(no).HeartRate*SET(no).SV/1000;
outdata{row,23}= bsa_1*SET(no).HeartRate*SET(no).SV/1000;

%PEF PFR
outdata{row,24} = SET(no).PFR;
outdata{row,25} = SET(no).PER;

%rmv...
outdata{row,26} = rvm;
outdata{row,27} = rvm*1.05;
outdata{row,28} = bsa_1*rvm*1.05;
outdata{row,29} = SET(no).RVEDV;
outdata{row,30} = bsa_1*SET(no).RVEDV;
outdata{row,31} = SET(no).RVESV;
outdata{row,32} = bsa_1*SET(no).RVESV;
outdata{row,33} = SET(no).RVSV;
outdata{row,34} = bsa_1*SET(no).RVSV;
outdata{row,35} = 100*SET(no).RVEF;

%--- Check if multiple scardata.
if length(scarno)>1
  %Find best scar data to take. Take those with endo and scar data
  scar2use = false(size(scarno));
  for sloop=1:length(scarno)
    scar2use(sloop) = existfunctions('existendo',scarno(sloop)) && not(isempty(SET(scarno(sloop)).Scar));
  end;
  scarno = scarno(scar2use);
  if length(scarno)>1
    mywarning('Detected multiple scar data. Taking first (arbitrary decision)');
    scarno = scarno(1); %Take first just arbitrary.
  end;
end;

if ~isempty(scarno)
  lvmde = SET(scarno).LVM(SET(scarno).EDT);
  outdata{row,36} = lvmde;
  if ~isempty(SET(scarno).Scar) %Moved this check here to get LVM out.
    outdata{row,37} = SET(scarno).Scar.Percentage;
    outdata{row,38} = lvmde*SET(scarno).Scar.Percentage/100;

    %Total extent
    tempstart = SET(scarno).StartSlice;
    tempend = SET(scarno).EndSlice;
    SET(scarno).StartSlice = 1;
    SET(scarno).EndSlice = SET(scarno).ZSize;
    [~,maxtransmurality,meantrans4infarct,totextent] = ...
      viability('calctransmuralityline',24,scarno);
    SET(scarno).StartSlice = tempstart;
    SET(scarno).EndSlice = tempend;
    outdata{row,39} = totextent;
    outdata{row,40} = meantrans4infarct;
    outdata{row,41} = max(maxtransmurality(:));
    
    if SET(scarno).Scar.UseWeighting
      %Total extent for weighted
      tempstart = SET(scarno).StartSlice;
      tempend = SET(scarno).EndSlice;
      SET(scarno).StartSlice = 1;
      SET(scarno).EndSlice = SET(scarno).ZSize;
      [~,~,maxtransmurality,meantrans4infarct,totextent] = ...
        viability('calctransmuralityweighted',24,scarno);
      SET(scarno).StartSlice = tempstart;
      SET(scarno).EndSlice = tempend;
      outdata{row,42} = totextent;
      outdata{row,43} = meantrans4infarct;
      outdata{row,44} = max(maxtransmurality(:));
    else
      outdata{row,42} = NaN;
      outdata{row,43} = NaN;
      outdata{row,44} = NaN;
    end

    %MO
    volscale = SET(scarno).ResolutionX*SET(scarno).ResolutionY*(SET(scarno).SliceThickness+SET(scarno).SliceGap)/1e3;
    outdata{row,45} = sum(SET(scarno).Scar.NoReflow(:))*volscale; %Region that is marked as where MO might reside.
    
    %For backwards compability.
    if isnan(SET(scarno).Scar.MOPercentage)
      viability('viabilitycalcvolume',scarno);
    end;
    
    outdata{row,46} = SET(scarno).Scar.MOPercentage;
  end;
else
  outdata{row,36} = NaN;
  outdata{row,37} = NaN;
  outdata{row,38} = NaN;
  outdata{row,39} = NaN;
  outdata{row,40} = NaN;
  outdata{row,41} = NaN;
  outdata{row,42} = NaN;
  outdata{row,43} = NaN;
  outdata{row,44} = NaN;
  outdata{row,45} = NaN;
  outdata{row,46} = NaN;
end;


%Flow
roipos = 47;
if ~isempty(flowno)
  for loop=1:length(flowno)
    NO=flowno(loop);
    
    if SET(flowno(loop)).RoiN>0
      %Calculate flow, call to get total flow.
      reportflow;
      tots = reportflow('gettotal');
      reportflow('close_Callback');
      
      for rloop=1:length(tots)
        if doheader
          outdata{1,roipos} = sprintf('%s[tot ml]',SET(NO).Roi(rloop).Name);
        end;
        outdata{row,roipos} = tots(rloop);
        roipos = roipos+1;
      end;
    end;
    
  end;
end;

%Measurements
measurementpos = roipos;
for sloop=1:length(SET)
  for mloop=1:length(SET(sloop).Measure)
    if doheader
      outdata{1,measurementpos} = [SET(sloop).Measure(mloop).Name '[mm]'];
    end;
    outdata{row,measurementpos} = SET(sloop).Measure(mloop).Length;
    measurementpos = measurementpos+1;
  end;
end;
  
%If called with no output arguments then copy to clipboard.
if nargout==0
  segment('cell2clipboard',outdata);
else
  varargout{1} = outdata;
end;