classdef waitbar < handle
  properties(SetAccess = private, GetAccess = private)
    h
  end
  
  methods(Access = public)
    function self = waitbar(n, msg)
      if n==1
        self.h = mywaitbarstart(2, msg);
        mywaitbarupdate(self.h);
      else
        self.h = mywaitbarstart(n, msg);
      end
    end
    
    function update(self)
      self.h = mywaitbarupdate(self.h);
    end
    
    function delete(self)
      mywaitbarclose(self.h);
      flushlog;      
    end
  end
end