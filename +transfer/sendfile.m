classdef sendfile < handle
  properties(SetAccess = private, GetAccess = private)
    data
    bytes_left
    chunksize
    sendfunc
  end
  
  methods(Access = public)
    function self = sendfile(sendfunc)
      self.data = struct('data', {}, 'pos', {}, 'fileid', {});
      self.chunksize = 2000000;
      self.bytes_left = self.chunksize;
      self.sendfunc = sendfunc;
    end
    
    function addfile(self, fileinfo)
      f = fopen(fileinfo.name, 'r');
      bytes = fread(f, inf, '*uint8');
      fclose(f);
      
      data_added = 0;
      while data_added ~= numel(bytes)
        data_added = data_added + ...
          self.add_data(bytes((1+data_added):end), data_added, fileinfo.id);
      end
    end
    
    function finalize(self)
      self.send();
    end
  end
  
  methods(Access = private)
    function send(self)
      if(numel(self.data) == 0)
        return
      end
      self.sendfunc(self.data);
      clear self.data;
      self.data = struct('data', {}, 'pos', {}, 'fileid', {});
      self.bytes_left = self.chunksize;
    end
    
    function data_added = add_data(self, raw_bytes, pos, id)
      % Find out how much data we want and b64encode it
      data_added = numel(raw_bytes);
      bytes = transfer.base64encode(raw_bytes(1:data_added));
      if(numel(bytes) > self.bytes_left)
        data_added = floor(self.bytes_left/1.3509)-100;
        bytes = transfer.base64encode(raw_bytes(1:data_added));
      end
      if(numel(bytes) > self.bytes_left)
        error('SEGMENT:PANIC', 'Couldn''t find good size');
      end
      
      self.data(end+1) = struct(...
        'data', bytes, ...
        'pos', pos, ...
        'fileid', id);
      self.bytes_left = self.bytes_left - numel(bytes);
      
      if(self.bytes_left < 10000)
        self.send();
      end
    end
  end

  methods(Static, Access = public)
  end
end