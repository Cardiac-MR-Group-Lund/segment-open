classdef thumbnails < handle
  properties(SetAccess = private, GetAccess = private)
    temppath
    files
    studyuid
    description
    patname
    wbfac
    didfail
    time
    sequencename
    venc
    docleanup = true;
  end
  
  methods(Access = public)
    function self = thumbnails(files, wbf)
      self.wbfac = wbf;
      wb = self.wbfac(4, 'Generating thumbnails...');
      self.temppath = tempname;
      mkdir(self.temppath);
      
      % Load the data
      loader = segloader();
      loader.adddicomfiles(files);
      
      wb.update();
      
      self.didfail = false;
      
      % May throw exception
      [~, data] = loader.render('', []);
      
      wb.update();
      
      % Get the uid, the patient name and the description
      self.studyuid = data(1).preview.StudyUID;
      self.patname = data(1).preview.PatientInfo.Name;
      if (~isfield(data(1).preview,'SequenceName') || strcmp(strtrim(data(1).preview.SequenceName), '')) && (~isfield(data(1).preview,'SeriesDescription') || strcmp(strtrim(data(1).preview.SeriesDescription), ''))
        self.description = 'No description';
        self.sequencename = 'No sequence name';
      elseif ~isfield(data(1).preview,'SequenceName') || strcmp(strtrim(data(1).preview.SequenceName), '')
        self.description = strtrim(data(1).preview.SeriesDescription);
        self.sequencename = 'No sequence name';
      else
        t = data(1).preview.AcquisitionTime;
        hours = floor(t/3600);
        minutes = floor((t-hours*3600)/60);
        seconds = floor(t-hours*3600-minutes*60);
        self.time = sprintf(...
          '%02d:%02d:%02d', ...
          hours, minutes, seconds);
        self.description = strtrim(data(1).preview.SeriesDescription);
        self.sequencename = strtrim(data(1).preview.SequenceName);
      end
      if size(data(1).IM, 5) ~= 1
        self.description = [self.description ' VENC data is present'];
      end
      
      self.files = {};
      picnum = 0;
      
      %Loop over rendered stacks
      for no = 1:numel(data)
        imsz = size(data(no).IM);
        if numel(imsz) < 5
          imsz(5) = 1;
          if numel(imsz) < 4
            imsz(4) = 1;
          end
        end
        
        % Get the thumbnails
        frames = imsz(3);
        slices = imsz(4);
        if frames>30
          %If more than 30 frames, then take middle, else take first. This is
          %then a good guess to get enddiastole.
          timeslice = round((imsz(3)+1)/2);
        else
          timeslice = 1;
        end;
        
        %Assign to wholepic.
        %If few timeframes, treat as slices
          if frames < 4 && slices == 1
            imsz = [imsz(1:2) 1 imsz(3) imsz(5)];
            wholepic = reshape(data(no).IM,imsz);
            frames = 1;
            slices = imsz(4);
          else
            wholepic = data(no).IM;
          end
          
        %--- Automatically crop the data if timeresolved.
        roisizeslice = 300; %mm
        roisizestack = 250; %mm
        if (slices>1) && (frames>15)
          xsize = roisizestack/data(no).preview.ResolutionX;
          ysize = roisizestack/data(no).preview.ResolutionY;
          wholepic = autocrop(wholepic,xsize,ysize);
        elseif (frames>15)
          xsize = roisizeslice/data(no).preview.ResolutionX;
          ysize = roisizeslice/data(no).preview.ResolutionY;
          wholepic = autocrop(wholepic,xsize,ysize);
        end;
        
        %Loop over slices
        for i=1:size(wholepic, 4)
          picnum = picnum + 1;
          
          %Extract data
          pic = wholepic(:, :, timeslice, i);
          
          %Autocontrast of the image
          pic = transfer.scaleim(pic);
          
          %Generate stationary thumbnails
          imwrite(pic, sprintf('%s%sthumb-%d-big.jpg', self.temppath, filesep, picnum));
          pic_small = imresize(pic, [128 round(128*size(pic,2)/size(pic,1))]);
          imwrite(pic_small, sprintf('%s%sthumb-%d-small.jpg', self.temppath, filesep, picnum));
          
          %Generate animations
          pic = transfer.scaleim(wholepic(:,:,:,i));
          pic = uint8(floor(pic*255));
          createanimgif(pic,sprintf('%s%sthumb-%d-anim.gif', self.temppath, filesep, picnum));
          
          %Send the files
          self.files{end+1} = {...
            sprintf('%s%sthumb-%d-big.jpg',self.temppath, filesep, picnum), ...
            sprintf('%s%sthumb-%d-small.jpg',self.temppath, filesep, picnum), ...
            sprintf('%s%sthumb-%d-anim.gif',self.temppath, filesep, picnum)};
        end
      end
      
      wb.update();
    end
    
    function r = getthumbs(self)
      r = self.files;
    end

    function r = getstudyuid(self)
      r = self.studyuid;
    end

    function r = getstudyname(self)
      r = self.patname;
    end

    function r = getdescription(self)
      r = self.description;
    end

    function delete(self)
      if self.docleanup
        transfer.deltree(self.temppath);
      end
    end
    
    function r = ok(self)
      r = ~self.didfail;
    end;
    
    function r = getsequence(self)
      r = self.sequencename;
    end
    
    function r = gettime(self)
      r = self.time;
    end
    
    function r = getvenc(self)
      r = self.venc;
    end
    
    function donotcleanup(self)
      self.docleanup = false;
    end
    
    
  end
  
end