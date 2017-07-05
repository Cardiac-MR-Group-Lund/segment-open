function matlab2ensightgold_multiple_parts(pathname,filename,timesteps,datastruct)
%MATLAB2ENSIGHTGOLD_MULTIPLE_PARTS   Saves Matlab data to a part in a EnsightGold fileformat.
%   Use of several 3D parts is supported, and is it
%   possibly to use several variables for each part. A limitation
%   is that the geometry is not allowed to be time-resolved.
%
%   MATLAB2ENSIGHTGOLD(PathName,FileName,TimeSteps,DataStruct)
%
%   - PathName       Pathname where files are saved. Set to '.' if
%                    current path is to be used.
%
%   - FileName       Filename of the ensight file (should
%                    not have suffix, i.e. .case).
%
%   - TimeSteps      vector of timesteps. zero if not timeresolved.
%
%   - DataStruct     A structure (array) where each element of the structure describes
%                    one part with potentially many variables. The struct needs to have
%                    the fields Data and VarName. Other fields are
%                    optional.
%
%   The fields are:
%
%   - Data           An array or a cell array of data where each cell
%                    contains data for one variable. Each element in must
%                    have the same size. For example a time-resolved 3D
%                    vector volume is given as x,y,z,t,n where n=1..3.
%
%   - VarName        A string or a cell array of strings the names for the
%                    variables in Data.
%
%   - [Size]         A vector describing the size of the data. For scalar
%                    variables this is the same as size(Data{1}), but
%                    differs for vector parts.
%
%   - [Desc]         string describing the part.
%
%   - [Origin]       three element vector of origin of data (x,y,z)
%
%   - [Delta]        three element vector describing spacing of data (x,y,z).
%
%   - [DirDim]       3x3 array describing the direction cosines of the data
%                    (column vectors)
%
%   Note1: If a variable called iblank is supplied as the first variable then
%   this variable is used to create an iblank array. This variable is
%   assumed to be one's or zeros. One is interior nodes and zeros are
%   exterior nodes and are blanked out. If the data is x,y,z,t then
%   iblank just need to be x,y,z.
%
%   Note2: NaN's are converted to a value twice as high as the highest
%   variable and reported as undefined to ensight.
%   
%   Note3: Known limitations are only 2D,2D+T,3D,3D+T parts are supported.

%  Version history
%  Created by Einar Brandt 2000-06-27
%  - EB 2000-07-04 Added support for iblank and changing geometry.
%  - EB 2000-07-03 Debugged for time-resolved parts
%  - PD 2006-09-07 *Changed from 'block' to 'curvelinear' and added support for multiple parts.
%                  *Temporarily removed support for iblank.
%                  *Changed the name of the function.
%  - PD 2006-10-09 *Added support for iblank
%  - EH 2006-11-28 Removed dependence of tensorarrays.
%                  Changed calling syntax
%                  Removed support for changing geometry
%                  Removed support for append allowed.
%                  Removed support for tensors.

%%% Error checking
if nargin < 4,
  error('Expected four input arguments.');
end;

if ~isa(filename,'char')
  error('FileName needs to be a string .');
end;

if ~isa(pathname,'char')
  error('PathName needs to be a string .');
end;

%Remove ending filesep if present
if ~isempty(pathname)
  if isequal(pathname(end),filesep)
    pathname = pathname(1:(end-1));
  end;
end;

if not(isa(datastruct,'struct'))
  error('Expected a struct as datastruct.');
end;

numparts = length(datastruct);
numsteps = length(timesteps); %Number of time steps

%--- General error checking for the datastruct
if ~isfield(datastruct,'Data')
  error('Data is a required field.');
end;

%Ensure that data is contained in a cell array
for partno=1:numparts
  if ~isa(datastruct(partno).Data,'cell')
    datastruct(partno).Data = {datastruct(partno).Data};
  end;
end;

%--- Fill in potentially missing fields
if ~isfield(datastruct,'Desc')
  for partno=1:numparts
    datastruct(partno).Desc = sprintf('Part description for part %d',partno);
  end;
end;

if ~isfield(datastruct,'DirDim')
  for partno=1:numparts
    datastruct(partno).DirDim = eye(3,3);
  end;
end;

if ~isfield(datastruct,'Origin')
  for partno=1:numparts
    datastruct(partno).Origin = zeros(1,3);
  end;
end;

if ~isfield(datastruct,'Delta')
  for partno=1:numparts
    datastruct(partno).Delta = ones(1,3);
  end;
end;

if ~isfield(datastruct,'Size')
  for partno=1:numparts
    datastruct(partno).Size = size(datastruct(partno).Data{1});
  end;
end;

if ~isfield(datastruct,'VarName')
  for partno=1:numparts
    datastruct(partno).VarName = cell(1,length(datastruct(partno).Data));
    for vloop=1:length(datastruct(partno).VarName)
      datastruct(partno).VarName{vloop} = sprintf('part%02dvar%02d',partno,vloop);
    end;
  end;
end;

%Ensure it is cell array
for partno=1:numparts
  if ~isa(datastruct(partno).VarName,'cell')
    datastruct(partno).VarName = {datastruct(partno).VarName};
  end;
end;

%Make place for field IsVector
for partno = 1:numparts
  datastruct(partno).IsVector = false(1,length(datastruct(partno).Data));
end;

%--- Check if valid values and check if iblanked and vectors
iblanked = false(1,numparts);
for partno = 1:numparts

  if ~isa(datastruct(partno).Desc,'char') %partdesc
    error(sprintf('Partdescription for part %d expected to be a string.',partno)); %#ok<*SPERR>
  end;

  if length(datastruct(partno).Desc)>78 %partdesc
    warning(sprintf('Too long partdescription for part %d, maximum 78.',partno)); %#ok<*SPWRN>
    datastruct(partno).Desc = datastruct(partno).Desc(1:78);
  end;

  if isempty(datastruct(partno).Size)
    datastruct(partno).Size = size(datastruct(partno).Data{1});
  end;
  
  if isempty(datastruct(partno).VarName)
    datastruct(partno).VarName = cell(1,length(datastruct(partno).Data));
    for vloop=1:length(datastruct(partno).VarName)
      datastruct(partno).VarName{vloop} = sprintf('part%02dvar%02d',partno,vloop);
    end;
  end;

  datastruct(partno).Origin = datastruct(partno).Origin(:);
  if length(datastruct(partno).Origin)~=3
    error(sprintf('Expected a three element vector as origin in part %d.',partno));
  end;

  datastruct(partno).Delta = datastruct(partno).Delta(:);
  if length(datastruct(partno).Delta)~=3
    error('Expected a three element vector as delta.');
  end;

  %Check and store if iblanked
  isiblanked(partno) = strcmp(datastruct(partno).VarName{1}, 'iblank');

  if isiblanked(partno)
    disp(sprintf('Using iblank array in part %d.', partno));
    if length(datastruct(partno).Data)<2
      error(sprintf('Part:%d. When iblanked the length of Data cell-array needs to be two or more.',...
        partno));
    end;
  end;

  %Check if match with temporal dimensions match
  if numsteps>1
    temp = datastruct(partno).Size;
    if not(isequal(temp(end),numsteps))
      error(sprintf(...
        'Number of timesteps and size in temporal dimension does not agree for part %d.',partno));
    end;
  end;

  if isempty(datastruct(partno).Data)
    error(sprintf('Empty data for part %d.',partno));
  end;

  %--- Check that all variables have the same size
  for vloop=1:length(datastruct(partno).Data)

    tempsize = size(datastruct(partno).Data{vloop});
    if length(tempsize)<4
      tempsize = [tempsize 1]; %#ok<AGROW>
    end;
    if length(tempsize)<4
      tempsize = [tempsize 1]; %#ok<AGROW>
    end;
    
    if length(tempsize)>length(datastruct(partno).Size)
      datastruct(partno).IsVector(vloop) = true;
    end;

    if not(datastruct(partno).IsVector(vloop))
      if not(isequal(datastruct(partno).Size,tempsize))
        error(sprintf('Size for variable %d in part %d does not match',vloop,partno));
      end;
    else
      %vector
      tempsize = tempsize(1:(end-1));
      if not(isequal(datastruct(partno).Size,tempsize))
        error(sprintf('Size for variable %d in part %d does not match',vloop,partno));
      end;
    end;
  end;

end %End of error checking of sizes

% ---- Create Case file

%Create and open the file
disp('Creating case file.');
casefilefid = fopen([pathname filesep filename '.case'],'w');
if isequal(casefilefid,-1)
  error(sprintf('Could not open %s for output.',[pathname filesep filename ...
    '.case']));
end;

%Write FORMAT
fprintf(casefilefid,'FORMAT\n\n');
fprintf(casefilefid,'type:\tensight gold\n\n');

%Write GEOMETRY section
if numsteps>1
  wildcard = '****';
else
  wildcard = '';
end;
fprintf(casefilefid,'GEOMETRY\n\n');
fprintf(casefilefid,'model:\t\t1\t%s.geo\n\n', filename);

%Write VARIABLE section
fprintf(casefilefid,'VARIABLE\n\n');
for partno = 1:numparts
  for vloop=1:length(datastruct(partno).Data)
    if ~(iblanked(partno)&(vloop==1))
      name = datastruct(partno).VarName{vloop};
      if ~datastruct(partno).IsVector(vloop)
        %scalar
        fprintf(casefilefid,'scalar per node:\t1\t%s\t%spart%03d_%s%s\n',name,filename,partno,name,wildcard);
      else
        %vector
        fprintf(casefilefid,'vector per node:\t1\t%s\t%spart%03d_%s%s\n',name,filename,partno,name,wildcard);
      end;
    end;
  end;
end
fprintf(casefilefid,'\n');

%Write TIME section
fprintf(casefilefid,'TIME\n');
fprintf(casefilefid,'time set:\t\t1\tModel\n'); %default 1
fprintf(casefilefid,'number of steps:\t%d\n',numsteps);
fprintf(casefilefid,'filename start number:\t1\n'); %default 1
fprintf(casefilefid,'filename increment:\t1\n'); %default 1
fprintf(casefilefid,'time values:\t%12.5e ',timesteps(1));

if numsteps>1
  for loop=2:numsteps
    fprintf(casefilefid,'%12.5e ',timesteps(loop)); %default 1
    if mod(loop,4)==0, fprintf(casefilefid,'\n'); end;
  end;
  fprintf(casefilefid,'\n');
end;

%Close file
fclose(casefilefid);

%%% ---- Create geometry file

newname = sprintf('%s%s%s.geo',pathname,filesep,filename);
geofilefid = fopen(newname,'w'); %  geofilefid = fopen(newname,'wb');
if isequal(geofilefid,-1)
  error(sprintf('Could not open %s for output.',newname));
end;

% Print Header
fprintf(geofilefid,'%-80s','C Binary');
fprintf(geofilefid,'%-80s','3D and 2D parts');
fprintf(geofilefid,'%-80s',sprintf('Exported by %s',mfilename));
fprintf(geofilefid,'%-80s','node id off');
fprintf(geofilefid,'%-80s','element id off');

for partno = 1:numparts

  dirdim = double(datastruct(partno).DirDim);
  dd1 = dirdim(:,1);
  dd2 = dirdim(:,2);
  dd3 = dirdim(:,3);

  voxspace = double(datastruct(partno).Delta); % voxel spacing
  origin = double(datastruct(partno).Origin);

  %  Transfer coordinates to patient reference system
  varsize = datastruct(partno).Size;
  disp(sprintf('Transfering coordinates in part %d to patient reference system', partno));
  if length(varsize)>2
    %3D
    [i j k] = ndgrid(0:(varsize(1)-1), 0:(varsize(2)-1), 0:(varsize(3)-1));
  else
    %2D
    [i j k] = ndgrid(0:(varsize(1)-1), 0:(varsize(2)-1), -.5:.5);
    voxspace(3) = 0.001;
    varsize(3) = 1;
    disp('   Note: The slice thickness of this part has been set to 1 mm.')
  end

  L = i*voxspace(1)*dd1(1) + j*voxspace(2)*dd2(1) + k*voxspace(3)*dd3(1) + origin(1);
  P = i*voxspace(1)*dd1(2) + j*voxspace(2)*dd2(2) + k*voxspace(3)*dd3(2) + origin(2);
  H = i*voxspace(1)*dd1(3) + j*voxspace(2)*dd2(3) + k*voxspace(3)*dd3(3) + origin(3);
  clear i j k
  L = L(:);
  P = P(:);
  H = H(:);

  %Start with the part
  fprintf(geofilefid,'%-80s', 'part');

  %Write part number
  fwrite(geofilefid, sprintf('%d',partno), 'int32');

  %Write description
  fprintf(geofilefid,'%-80s', datastruct(partno).Desc);

  %Write mark to be curvelinear grid.
  if iblanked(partno)
    fprintf(geofilefid,'%-80s', 'block iblanked');
  else
    fprintf(geofilefid,'%-80s', 'block');
  end;

  %Write size
  fwrite(geofilefid,varsize(1:3),'int32');

  %Write grid
  fwrite(geofilefid,[L; P; H],'float32');

  %Write iblank
  if iblanked(partno)
    temp = datastruct(partno).Data{1};
    fwrite(geofilefid,temp(:),'int32');
    clear temp; %Free some memory
  end;
end

% Close geo file
fclose(geofilefid);

%  end; %Timeframe loop
clear L P H %Free some memory

%%% ---- Create variable files
for partno = 1:numparts
  for vloop=1:length(datastruct(partno).Data)
    if ~(iblanked(partno)&(vloop==1))	%Do not export variable if this is iblank
      temp = datastruct(partno).Data{vloop};
      disp(sprintf('Finding maximum value for part %d variable %s',partno,datastruct(partno).VarName{vloop}));
      %Find maximum to get what is NaN's
      maxvalue = max(temp(:));
      maxvalue = maxvalue*2; %Take a larger value
      
      %Replace NaN's with maxvalue
      temp(isnan(temp))=maxvalue;
      name = datastruct(partno).VarName{vloop};

      %--- Loop over timeframes
      for timeframe=1:numsteps
        %Create and open the file
        newname = [pathname filesep filename sprintf('part%03d',partno) '_' name ];
        if numsteps>1
          newname = [newname sprintf('%04d',timeframe)];
        end;
        disp(sprintf('Creating variable file [%s]',newname));

        varfilefid = fopen(newname,'w');
        if isequal(varfilefid,-1)
          error(sprintf('Could not open %s for output.',newname));
        end;

        %Write description
        fprintf(varfilefid,'%-80s', datastruct(partno).Desc);

        %Write part
        fprintf(varfilefid,'%-80s', 'part');

        %Write part number
        fwrite(varfilefid,sprintf('%d',partno),'int32');

        %Write block data
        fprintf(varfilefid,'%-80s', 'block undef');

        %Write undef value
        fwrite(varfilefid,maxvalue,'float32');  %Not a number

        if ~datastruct(partno).IsVector(vloop)
          % Data is scalar valued
          if numsteps>1
            if ndims(temp)==3
              fwrite(varfilefid,temp(:,:,timeframe),'float32');
            else
              fwrite(varfilefid,temp(:,:,:,timeframe),'float32');              
            end;
          else
            fwrite(varfilefid,temp,'float32');
          end;
        else
          % Data is vector-valued

          % Preparations for changing vectors to the patient (Left,
          % Postierior, Head) coordinate system.
          dirdim = double(datastruct(partno).DirDim);
          dd1 = dirdim(:,1); % first index is along this direction in LPH
          dd2 = dirdim(:,2); % second index 
          dd3 = dirdim(:,3); % third index
          
          %first x-coords, y-coords, z-coords
          if numsteps>1
            % Field is time-dependent
            if ndims(temp)==4
              %2D vector part
              
              %%%%% H�R SKALL DE ROTERAS %%%%%
              fwrite(varfilefid,temp(:,:,timeframe,1),'float32');
              fwrite(varfilefid,temp(:,:,timeframe,2),'float32');
              fwrite(varfilefid,temp(:,:,timeframe,3),'float32');              
            else
              %3D vector part
                            
              % v{L,P,H}: velocities in Left, Posterior and
              % Head-direction.
              v1 = temp(:,:,:,timeframe,1);
              v2 = temp(:,:,:,timeframe,2);
              v3 = temp(:,:,:,timeframe,3);
              vL = v1*dd1(1) + v2*dd2(1) + v3*dd3(1);
              vP = v1*dd1(2) + v2*dd2(2) + v3*dd3(2);
              vH = v1*dd1(3) + v2*dd2(3) + v3*dd3(3);
              
              fwrite(varfilefid,vL,'float32');
              fwrite(varfilefid,vP,'float32');
              fwrite(varfilefid,vH,'float32');
            end;
          else
            % Field is constant in time, output the same field for every
            % timestep.
            if ndims(temp)==3
              %2D vector part
              
              %%%%% H�R SKALL DE ROTERAS %%%%%              
              fwrite(varfilefid,temp(:,:,1),'float32');
              fwrite(varfilefid,temp(:,:,2),'float32');
              fwrite(varfilefid,temp(:,:,3),'float32');              
            else
              %3D vector part
              
              % Note (JT): this is untested, but should work
              v1 = temp(:,:,:,1);
              v2 = temp(:,:,:,2);
              v3 = temp(:,:,:,3);
              vL = v1*dd1(1) + v2*dd2(1) + v3*dd3(1);
              vP = v1*dd1(2) + v2*dd2(2) + v3*dd3(2);
              vH = v1*dd1(3) + v2*dd2(3) + v3*dd3(3);
              
              fwrite(varfilefid,vL,'float32');
              fwrite(varfilefid,vP,'float32');
              fwrite(varfilefid,vH,'float32');
            end;
          end;
        end;

        %Close file
        fclose(varfilefid);

      end; %loop over timeframes
    end; %End if iblank statment
  end; %loop over variables
end; %loop over parts
