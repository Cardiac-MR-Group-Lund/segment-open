classdef gui < maingui
  
  methods
    
    %-------------------------------
    function g = gui(programversion)
    %-------------------------------
    %Constructor
    
    %Call superclass constructor
    g = g@maingui(programversion);
    
    end
    
    %----------------------
    function initLogFile(g)
    %----------------------
    %This overloaded method sets the name of the log-file.
    pathname = g.getpreferencespath;
    g.LogFile = [pathname filesep sprintf('segmenttransferlog_%s.log',datestr(now,'yyyymmddHHMMSS'))];
    fid = fopen(g.LogFile,'w');
    if isequal(fid,-1)
      myfailed(dprintf('Could not create .log file %s.',g.LogFile));
    else
      fclose(fid); %Close the file, to make it empty
      g.startlog(g.LogFile); %Start diary.
      fprintf(1,'Segment Transfer Version %s\n',g.ProgramVersion);
    end;
    end
    
    %---------------------------
    function loadguipositions(g)
    %---------------------------
    %Load prestored positions of GUI's.
    %Overload this with empty function, do not need to load gui positions
    end;
    
    %----------------------
    function initmaingui(g)
    %----------------------
    %Initiates main GUI. Overloads mainGUI method to use Segment-transfer GUI instead
    g.GUI.Segment=mygui(['+transfer' filesep 'segmenttransfer.fig']);%Load custom GUI
    g.fig = g.GUI.Segment.fig;
    g.Handles = killhandles(g.fig,'segment.fig');
    g.ProgramName = 'Segment-Transfer';
    
    gui = g.GUI.Segment;
    
    %Create the progressbar
    gui.handles.overallpatch = patch(...
      [0 0 0 0],...
      [0 1 1 0],...
      [1 1 1 1],...
      'parent',gui.handles.overallaxes);
    hold(gui.handles.overallaxes,'on');    hold(gui.handles.overallaxes,'on');
    gui.handles.overallline = plot(gui.handles.overallaxes,[0 0 1 1 0],[0 1 1 0 0],'k-');
    hold(gui.handles.overallaxes,'off');
    xlim(gui.handles.overallaxes,[0 1]);
    axis(gui.handles.overallaxes,'off');
    set(gui.handles.overalltext,'string','');

    set(gui.handles.overallline,'visible','off');
    set(gui.handles.overallpatch,'visible','off');
    set(gui.handles.overallaxes,'visible','off');

    %Create the detail waitbar
    gui.handles.detailpatch = patch(...
      [0 0 0 0],...
      [0 1 1 0],...
      [1 1 1 1],...
      'parent',gui.handles.detailaxes);
    hold(gui.handles.detailaxes,'on');
    gui.handles.detailline = plot(gui.handles.detailaxes,[0 0 1 1 0],[0 1 1 0 0],'k-');
    hold(gui.handles.detailaxes,'off');
    xlim(gui.handles.detailaxes,[0 1]);
    axis(gui.handles.detailaxes,'off');
    
    set(gui.handles.detailtext,'string','');
    set(gui.handles.detailline,'visible','off');
    set(gui.handles.detailpatch,'visible','off');   
    set(gui.handles.detailaxes,'visible','off');
    
    %Ensure that the credentials file exists
    if isempty(transfer.getcredsfile)
      mywarning(...
        ['The credentials file is not installed. This file contains necessary log-in '...
        'information. You should have received this file in email. '...
        'Please contact the study core-lab if you are missing the file, or support@medviso.com. ' ...
        'You will now be prompted to locate the file.']);

      pause(0.1);

      %Ask for file
      g.askforfile();
    end;
    
    if isempty(transfer.getcredsfile)
      myfailed('Please fix problem with missing credentials and restart application.');
    end;

    end
    
    %----------------------
    function inittoolbar(g) %#ok<INUSD>
    %----------------------
    %Initiate toolbars. Overloads mainGUI method.
    end
    
    %---------------------------
    function updateicons(g,mode) %#ok<INUSD>
    %---------------------------
    %Overloads mainGUI method
    
    %Reset Tool structure
    g.Tools = [];
    end
    
    %------------------
    function askforfile(g) %#ok<INUSD>
    %------------------
      %Ask for creds file and copies it to pwd

      [filename,pathname] = myuigetfile('*.transfercreds','Select the file transfercreds file');
      if isequal(pathname,-1) || isequal(filename,-1)
        myfailed('Aborted.')
        return;
      end;

      if ~exist([pathname filesep filename],'file')
        myfailed('File does not exist. Aborted.');
      end;
      ok = copyfile([pathname filesep filename],[pwd filesep filename]);
      if ~ok
        myfailed(dprintf('Could not copy file. Folder %s write protected?',pwd));
      else
        mymsgbox('File copy successfull.');
      end;
    end
    
    %----------------------
    function updatetitle(g)
    %----------------------
    %Updates titleline of main GUI. Overloads method in SegmentGUI
   
    set(g.fig,'Name',['Segment Transfer v' changelog]);
    end
    
  end
  
  methods(Static)
    
    %-----------------------
    function dispwelcometext
    %-----------------------
    %Overloads method in mainGUI
    end
        
    %--------------------
    function versionhello
    %--------------------
    %Overloads method in mainGUI
    disp(['Segment Transfer v' changelog '.'])
    end
    
      
    %----------------------------------------------
    function h = mywaitbarstart(iter,stri,varargin)
    %----------------------------------------------
    global DATA

    %Just store it for future work
    gui = DATA.GUI.Segment;   
    gui.numiter = iter;
    gui.iter = 0;
    gui.detailstri = stri;
    
    %Update text
    set(gui.handles.detailtext,'string',stri);
    
    %Make it visible
    set([...
      gui.handles.detailline...
      gui.handles.detailpatch ...
      gui.handles.detailtext],'visible','on');
    h = []; 
    end;
    
    %--------------------------------
    function mywaitbarclose(varargin)
    %--------------------------------
    global DATA
    
    gui = DATA.GUI.Segment;   
    
    %Make it invisible
    set([...
      gui.handles.detailline...
      gui.handles.detailpatch ...
      gui.handles.detailtext],'visible','off');    
    
    end;
    
    %-------------------------------------
    function h = mywaitbarupdate(varargin)
    %-------------------------------------
    global DATA
    
    gui = DATA.GUI.Segment;
    
    gui.iter = gui.iter+1;
    f = gui.iter/gui.numiter;
    
    set(gui.handles.detailpatch,'xdata',[0 0 f f]);
    drawnow;
    
    h = [];
    end;
    
    
  end
  
end