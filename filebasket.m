classdef filebasket < handle
  properties
    fid
    len
  end
  
  methods
    function self = filebasket(filename, mode)
      self.fid = fopen(filename, mode);
      if self.fid == -1
        error('SEGMENT:ERROR', 'Couldn''t open file');
      end
      fseek(self.fid, 0, 'eof');
      self.len = ftell(self.fid);
      fseek(self.fid, 0, 'bof');
    end
    
    function add(self, data)
      if strcmp(class(data), 'memorybasket')
        data = data.render();
      end
      fwrite(self.fid, data);
      self.len = self.len + length(data);
    end
    
    function r = pop(self, len)
      len = double(len);
      r = fread(self.fid, len, '*uint8')';
      if len > length(r)
        error('SEGMENT:ERROR', ...
          'Invalid stream: Not Enough bytes left');
      end
    end
    
    function r = bytes_left(self)
      r = self.len - ftell(self.fid);
    end
    
    function close(self)
      fclose(self.fid);
    end
  end
end