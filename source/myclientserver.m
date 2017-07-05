%MYCLIENTSERVER Class to create interface object to communicate with a 
%computational server that can perform tasks in background.
%
%B = MYCLIENTSERVER(NAME) Sets up a computational client with the name NAME.
%
%METHODS
%
%DISPLAY
%  Overloads display.
%
%SEND(B,FCN,CMD,ARG). Sends a non blocking message. 
%  FCN is function handle / function that is executed upon return of 
%  message. The function handle is called with:
%
%  1) Name of slave
%  2) Received status, 0 = sucess
%  3) Received command
%  4) Received arguments, often empty, depends on command
%
%  Note that there is NO buffering of commands, so you can not stack
%  multiple commands. If you do not give the slave time to read the message
%  before you send another message.
%
%SET(B,'PROPERTY',VALUE)
%   Set property values
%
%OK = PING(B,TIMEOUT). 
%  Pings client.
%
%OK = ALIVE(B). 
%  Test if server/slave is alive
%
%DELETE(B). 
%  Destructor for class. Sends kill command to server/slave.
%
%B = MYCLIENTSERVER instanciate object.
%
%[STATUS,RES,<ARGS>] = SENDBLOCK(B,CMD,ARG). 
%  Sends a blocking message.
%
%VALUE = GET(B,'PROPERTY')

%Communication protocol:
%-----------------------
%See decodemessage and encodemessage
%
%Definitions:
%------------
%- Client (i.e program that uses the slave).
%- Slave (i.e program that performs the tasks).
%
%See also SLAVE.

%Einar Heiberg.
%Modified by Jane Sjögren and Jonatan Wulcan.

classdef (InferiorClasses = {?timer}) myclientserver < handle %Inherits from handles to get persistent objects.

  properties (SetAccess = 'private',Hidden)
    Name = '';
    Path = '';
    Pause = 0.05;
    Verbose = true;
    PollTimer = [];
    AtReceive = [];
    LastMessageId = 0;
    LastReceivedId = 0;
  end;

  methods (Access = 'public')
    
    %-----------------------    
    function b = myclientserver(name)
      %Constructor method

      if nargin~=1
        myfailed('Expected one input argument.');
        return;
      end;
      
      %Clear spaces and bars from name
      name = removeforbiddenchars(name);

      b.Name = name;
      b.Path = getpreferencespath;
      b.Verbose = true;
      b.PollTimer = timer(...
        'Name','internalpolltimer',...
        'Period',0.05,...
        'StartDelay',0.05,...
        'executionmode','fixeddelay',...
        'TimerFcn',{@poll,b});
            
      %Ensure inbox is empty
      outfid = fopen(getinboxname(b),'w');
      if outfid>0
        fprintf(outfid,'\n');        
        fclose(outfid);
      else        
        pause(b.Pause);
        outfid = fopen(getinboxname(b),'w');
        if outfid>0
          fprintf(outfid,'\n');
          fclose(outfid);
        else
          error('Could not open slave''s inbox.');
        end;
      end;
      
      %Ensure outbox is empty
      infid = fopen(getoutboxname(b),'w');
      if infid>0
        fprintf(infid,'\n');
        fclose(infid);
      else        
        pause(b.Pause);
        infid = fopen(getoutboxname(b),'w');
        if infid>0
          fprintf(infid,'\n');
          fclose(infid);
        else
          error('Could not open slave''s outbox.');
        end;
      end;      
            
      start(b.PollTimer);
      
      %Start the server/slave
      system(sprintf('slave%s %s &',myexecutableext,name));
      
      pause(3);
    end

    %-----------------------
    function send(b,fcn,cmd,args)
      %SEND(B,FCN,CMD,ARGS) Send message to slave B.
      %FCN is function handle to be executed when answer returns
      %CMD is command to send
      %ARGS are optional arguments to send.
      
      %--- Error checking
      if nargin<3
        error('Expected at least three input arguments.');
      end;

      if nargin<4
        args = '';
      end;
           
      %Compose string.
      b.LastMessageId = b.LastMessageId + 1;
      stri = encodemessage(b.Name,0,b.LastMessageId,cmd,args);

      if b.Verbose
        tempstri = stri;
        tempstri(tempstri==sprintf('\n')) = ' ';
        disp(dprintf('Sending string: %s',tempstri)); %#ok<DSPS>
      end;

      %--- Setup struct for AtReceive
      if isempty(b.AtReceive)
        b.AtReceive = [];
        b.AtReceive.ID = b.LastMessageId;
        b.AtReceive.Fcn = fcn;
      else
        n = length(b.AtReceive);
        b.AtReceive(n+1).ID = b.LastMessageId;
        b.AtReceive(n+1).Fcn = fcn;
      end;
      
      %--- Send message
      outfid = fopen(getinboxname(b),'w');
      if outfid<1
        pause(b.Pause);
      else
        fprintf(outfid,'%s\n',stri);
        fclose(outfid);
        return; %to avoid second attempt
      end;
      
      outfid = fopen(getinboxname(b),'w');
      if outfid<1
        myfailed('Could not send message to slave.');
      else
        fprintf(outfid,'%s\n',stri);
        fclose(outfid);
        return; 
      end;            

    end; %end method send

    %-----------------------
    function ok = sendblock(b,cmd,arg,varargin)
      %Pings the slave, return true if response within timeout
      persistent exitnow
      
      ok = false;
      
      if nargin<4
        %This clause is the normal exectution pathway
        exitnow = false;
        
        send(b,{@b.sendblock 1},cmd,arg);
        
        while not(exitnow)
          pause(0.05);
        end;
        
        ok = exitnow;
        
      else
        %this clause is when called from poll.
        exitnow = true;
      end;
          
    end; %end method sendblock
    
    %-----------------------
    function display(b)
      %Display method 
      disp(dprintf('myclient with name %s',b.Name)); %#ok<DSPS>
    end; %end method display

    %-----------------------
    function ok = alive(b)
      %Check if slave is alive.
      ok = ping(b,2);
    end; %method alive
    
    %-----------------------    
    function ok = ping(b,timeout,varargin)
      %Pings the slave, return true if response within timeout
      persistent exitnow
      
      ok = false;
      
      switch nargin        
        case {1,2}
          %This clause is the normal exectution pathway
          starttime = now;
          exitnow = false;
          
          if nargin==1
            timeout = 10;
          end;

          send(b,{@b.ping 1},'ping','');

          while not(exitnow) && ((now-starttime)*24*3600)<timeout
            pause(0.05);
          end;

          ok = exitnow;
          
        otherwise
          %this clause is when called from poll.
          exitnow = true;
      end;
          
    end; %end method ping
  
    %-----------------------    
    function [varargout] = delete(b)
      %Delete method 
      
      try
        stop(b.PollTimer);
        send(b,0,'kill',''); %Try die in style
      catch me
        disp('Failed trying to delete myclientserver object.');
        mydispexception(me);
      end;

      %Delete the files.
      pause(0.5);
      try
        delete(getinboxname(b));
        delete(getoutboxname(b));
      catch %#ok<CTCH>
        disp('Could not delete connection.');
      end;
      
      %if used as b = delete(b);
      if nargout>0
        varargout = cell(1,nargout);
        varargout{1} = [];
      end;
      
    end; %end method delete
    
    %-----------------------    
    function out = get(b,parameter)
      %Get parameters
      
      if nargin<2
        error('Expected two input arguments.');
      end;
      
      switch lower(parameter)
        case 'name'
          out = b.Name;
        case 'verbose'
          out = b.Verbose;
        otherwise
          error(dprintf('Invalid parameter %s.',parameter)); %#ok<SPERR>
      end;
    end; %get method
      
    %-----------------------
    function set(b,varargin)
      %Set paramters 
      
      nargs = length(varargin);
      if ~isequal(round(nargs/2)*2,nargs)
        error('Expected even parameter pairs.');        
      end;
      
      npar = nargs/2;
      
      %Loop over parameter pairs.
      for loop=1:npar
        par = varargin{loop*2-1};
        value = varargin{loop*2};
        switch lower(par)
          case {'name'}
            error('Can not change name or server port. Please delete object instead and recreate.');            
          case 'verbose'
            if ~islogical(value)
              error('Expected logical value as verbose');
            end;
            if numel(value)~=1
              error('Expected scalar value as verbose');
            end;            
            b.Verbose = value(1);
          otherwise
            error(dprintf('Invalid parameter %s.',par)); %#ok<SPERR>
        end;
      end;
    end; %set method  

    %-----------------------    
    function poll(obj, event, b, varargin)        %#ok<INUSL>
      %This function is called by the timer object. It reads from the port
      %and if something has arrived it looks to which call it belongs.
      %When it calls the other function it adds 
      %1) Name of slave
      %2) Received status, 0 = sucess
      %3) Received command
      %4) Received arguments
      
      %Read string
      infid = fopen(getoutboxname(b),'r'); %get the slaves outbox to read from.
      if infid>0        
        
        %Try to read from the file
        fullstri = '';
        stri = '';
        try
          while ~feof(infid) && ~isnumeric(stri);
            stri = fgetl(infid);
            if ~isnumeric(stri) && ~isempty(stri)
              fullstri = [fullstri stri sprintf('\n')]; %#ok<AGROW>
            end;
          end;
          fclose(infid);
          if ~ischar(fullstri)
            fullstri = '';
          end;
        catch %#ok<CTCH>
          fullstri = '';
        end;
      else        
        fullstri = '';
      end;
   
      %Avoid []
      if ~ischar(fullstri)
        fullstri = '';
      end;
      
      %If the string is not empty then we received something.
      if ~isempty(fullstri)

        %Clear it. This is dangerous!!!
        %infid = fopen(getoutboxname(b),'w');
        %fprintf(infid,'\n');
        %fclose(infid);
        
        %Decode the message
        [rname,rs,rid,rcmd,rarg] = decodemessage(fullstri);
        
        if rid <= b.LastReceivedId
          %Already received this
          return;
        else
          b.LastReceivedId= rid;
        end;

        if b.Verbose
          disp(dprintf('Received:%s',fullstri)); %#ok<DSPS>
        end;
        
        %AtReceive is a list of registered commands in mybatchworkerobject
        if isempty(b.AtReceive)
          error('No commands registered.');
        end;
        
        %Loop over registered commands to try to find a match
        n = length(b.AtReceive);

        logind = true(1,n); %Keep track if matches, false = match..
        for loop=1:n
          if isequal(b.AtReceive(loop).ID,rid)
            logind(loop)=false;
            
            %Extract function
            fcn = b.AtReceive(loop).Fcn;
            
            %Execute function
            if iscell(fcn)                            
              feval(fcn{1},rname,rs,rcmd,rarg,fcn{2:end});              
            else
              feval(fcn,rname,rs,rcmd,rarg);
            end;
          end; %matching command received
        end;
        
        %Remove this command from the list
        b.AtReceive = b.AtReceive(logind);
      end;
      
    end; %End of method poll
    
    %------------------------------------
    function name = getoutboxname(b)
      %Returns the name of the slaves outbox      
      name = [b.Path filesep '.' b.Name '.outbox'];
    end;

    %------------------------------------
    function name = getinboxname(b)
      %Returns the name of the slaves inbox
      name = [b.Path filesep '.' b.Name '.inbox'];
    end;
 
  end; %End private methods
end