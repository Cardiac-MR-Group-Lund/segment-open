classdef memoryhandle < handle
    properties
        mem
        pos
    end
    
    methods
        function self = memoryhandle(size)
            self.mem = repmat(uint8(0), 1, size);
            self.pos = 1;
        end
        
        function write(self, data)
            self.mem(self.pos:(self.pos+length(data)-1)) = data;
            self.pos = self.pos + length(data);
        end
        
        function mem = read(self)
            mem = self.mem;
        end
    end
end
