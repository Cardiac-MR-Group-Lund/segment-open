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
    end
    
    %Starts with Modality and contains a lot of .
    if length(filename)>4
      if ...
          contains(filename,'MR.') || ...
          contains(filename,'MRe.') || ...          
          contains(filename,'SC.') || ...                    
          contains(filename,'CT.') || ...
          contains(filename,'CTe.') || ...                    
          contains(filename,'PT.') || ...
          contains(filename,'NM.') || ...
          contains(filename,'OT.') || ...
          contains(filename,'XR.')
        if sum(filename=='.')>=6
          isdicom = true;
          return;
        end
      end
      if  contains(filename,'PSg.') || ...
          contains(filename,'REG.') || ...
          contains(filename,'SRc.') || ...
          contains(filename,'SRd.') || ...
          contains(filename,'SRe.') || ...
          contains(filename,'SRt.')
        
        isdicom = false;
        return;
      end
      if isequal(filename(1:3),'im-')
        isdicom = true;
      end
    end
        
    
    %Files with only numbers and no extension
    if isequal(filename,removechars(filename)) && isempty(ext)
      isdicom = true;
      return;
    end
    
    % It's not a dicom file if the file ends with '.???'
    if regexp(rawfilename, '\....$') 
      if ~isequal(lower(ext),'.ima')
        isdicom = false;
        return;
      end
    end
    if isequal(lower(ext),'.cache')
      isdicom = false;
      return;
    end
%     if  filename=='.'
%       isdicom =false;
%       return;
%     end
    
    if  contains(filename,'VERSION')
      isdicom = false;
      return;
    end
    
    if not(...
        isequal(filename,'thumbs.cache') || ...
        isequal(filename,'folders.cache') || ...
        isequal(filename,'dicom.cache'))
 
      try        
        segdicomtags.readfiles({rawfilename});
        isdicom = true;
        return;
      catch
        try
          temp = dicominfo(rawfilename);
          if isfield(temp,'ImagingFrequency') %Fixed for Toshiba/Canon 3T MRI sending images directly from MRISystem
            if (round(temp.ImagingFrequency) == 123) && ( strcmpi(temp.Manufacturer,'toshiba')|| strcmpi(temp.Manufacturer,'canon'))
              if strcmp(temp.ImageType,'ORIGINAL\PRIMARY\OTHER')&&strcmp(temp.SequenceName,'PSMRA')
                isdicom = true;
              elseif strcmp(temp.ImageType,'ORIGINAL\PRIMARY\OTHER') || strcmp(temp.ImageType,'DERIVED\SECONDARY\SHIMMING')
                isdicom = false;
              else
                isdicom = true;
              end
            else
              isdicom = true; %For all other MR
            end
          else
            isdicom = true; %In the case that it is not MR but CT
          end
          
          return;
        catch
          %mywarning('The File is not DICOM');
          isdicom = false;
          return;
        end
        %       temp=dicominfo(rawfilename);
        %       if isfield(dicominfo(rawfilename),'Format')
        %         if strcmp(temp.Format,'DICOM')
        %           isdicom = true;
        %           return;
        %         end
        %       end
      end
      
      
      % if we're still not sure try and load the file
      %disp(rawfilename);
    end
    try
      segdicomtags.readfiles({rawfilename});
      isdicom = true;
      return
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
    end
    
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

