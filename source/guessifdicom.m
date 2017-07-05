function isdicom = guessifdicom(filename,mode)
%Guess if a file is a dicomfile or not.
%ISDICOM = GUESSIFDICOM(FILENAME,MODE)
%
%guess if dicom file based on three different modes 
%1) heuristics (default) 
%2) file extension==.dcm or .ima
%3) take all files(not recommended)

%Einar Heiberg

if nargin == 1
  mode = 1;
end

isdicom = false;
switch mode
  case 1
    rawfilename = filename;
    [~,filename,ext] = fileparts(filename); 
    ext = lower(ext);
    if isequal(ext,'.dcm') || isequal(ext,'.segdicom') || isequal(ext,'.dicom')
      isdicom = true; %It is a .dcm file
      return;
    end;
    
    %Starts with Modality and contains a lot of .
    if length(filename)>4
      if ...
          ~isempty(strfind(filename,'MR.')) || ...
          ~isempty(strfind(filename,'MRe.')) || ...          
          ~isempty(strfind(filename,'SC.')) || ...                    
          ~isempty(strfind(filename,'CT.')) || ...
          ~isempty(strfind(filename,'CTe.')) || ...                    
          ~isempty(strfind(filename,'PT.')) || ...
          ~isempty(strfind(filename,'NM.')) || ...
          ~isempty(strfind(filename,'OT.')) || ...
          ~isempty(strfind(filename,'XR.'))
        if sum(filename=='.')>=6
          isdicom = true;
          return;
        end;        
      end;
      if isequal(filename(1:3),'im-')
        isdicom = true;
      end;
    end;
        
    
    %Files with only numbers and no extension
    if isequal(filename,removechars(filename)) && isempty(ext)
      isdicom = true;
      return;
    end;
    
    % It's not a dicom file if the file ends with '.???'
    if regexp(rawfilename, '\....$') 
      if ~isequal(lower(ext),'.ima')
        isdicom = false;
        return;
      end
    end;
    
    % if we're still not sure try and load the file
    %disp(rawfilename);
    
    try
      segdicomtags.readfiles({rawfilename});
      isdicom = true;
    catch e
      if not(strcmp(e.identifier, 'SEGMENT:ERROR'))
        rethrow(e);
      end
    end
    
  case 2
    [~,~,ext] = fileparts(filename); 
    ext = lower(ext);
    if isequal(ext,'.dcm') || isequal(ext,'.ima')
      isdicom = true; %It is a .dcm file
    end;
    
  case 3
    if not(...
        isequal(filename,'thumbs.cache') || ...
        isequal(filename,'folders.cache') || ...
        isequal(filename,'dicom.cache'))
      isdicom = true;
    end
  otherwise
    error('Expected mode to be 1..3');    
end

