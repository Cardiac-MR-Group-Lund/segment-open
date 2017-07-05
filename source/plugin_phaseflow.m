function varargout = plugin_phaseflow(fcn,varargin)
%Function to calculate flow from phase images only

%Nils Lundahl

global SET NO

if nargin==0
  myfailed('Expects at least one input argument.');
  return;
end;

switch fcn
  case 'getname'
    varargout = cell(1,1);
    
    %Register callbacks
    uimenu(varargin{1},'Label','Create magnitude stack for flow calculations','Callback','plugin_phaseflow(''go_Callback'')');
    set(varargin{1},'Callback','');

    %Return title
    varargout{1} = 'Flow from phase image plugin';
    
  case 'getdependencies'
    %No dependencies
    varargout = cell(1,4);
    
  case 'go_Callback'
    %Create magnitude stack and modify phase stack for calculations
    
    if ~yesno('This might yield incorrect results for images not from Philips scanners. Are you sure you want to go on?')
      return
    end
        
    %Modify phase image to contain proper flow data
    imout = calcfunctions('calctruedata',SET(NO).IM,NO);
    venc = SET(NO).VENC;
    SET(NO).IntensityScaling = 2*venc;
    SET(NO).IntensityOffset = -venc;
    im = imout/(2*venc)+0.5;
    SET(NO).IM = im;
    
    %Create magnitude image stack
    immag = 2*abs(im-0.5);
    newno = length(SET)+1;
    SET(newno) = SET(NO);
    SET(newno).IM = immag;
    SET(newno).IntensityScaling = venc;
    SET(newno).IntensityOffset = 0;
    
    %Create linking structure
    linkies = [NO newno];
    SET(NO).Linked = linkies;
    SET(newno).Linked = linkies;
    SET(NO).Parent = newno;
    SET(newno).Children = NO;
    
    %Create flow struct
    flowstruct = struct;
    flowstruct.MagnitudeNo = newno;
    flowstruct.PhaseNo = NO;
    flowstruct.PhaseX = [];
    flowstruct.PhaseY = [];
    flowstruct.PhaseCorr = [];
    flowstruct.Angio = [];
    flowstruct.VelMag = [];
    SET(NO).Flow = flowstruct;
    SET(newno).Flow = flowstruct;
    
    drawfunctions('drawthumbnails');
end