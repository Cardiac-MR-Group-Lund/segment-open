function varargout = shuntvalve(varargin)
%Functionality for Qp/Qs and Shunt and Valve analysis
%Nils Lundahl

if nargin == 0
  varargin = {'qpqs'};
end

macro_helper(varargin{:});
[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%-----------------------------------------------------------------
function [svaorta,svpulmo,fwaorta,bwaorta,fwpulmo,bwpulmo] = getsv
%-----------------------------------------------------------------
%Calculate and return stroke, forward and backward volumes both for the 
%aorta and the pulmonary artery.
global SET
%Find magnitude image on flow
[~,~,flowno] = findfunctions('findno');
flowmagnitudeno = [];
for loop = flowno
  if loop == SET(loop).Flow.MagnitudeNo
    flowmagnitudeno = [flowmagnitudeno loop]; %#ok<AGROW>
  end
end

%Find proper vessels for analysis
svaorta = [];
svpulmo = [];
fwaorta = [];
bwaorta = [];
fwpulmo = [];
bwpulmo = [];
for floop =  flowmagnitudeno
  if SET(floop).RoiN > 0
    reportflow('init',floop);
    rois = reportflow('getroiname');
    aaf = find(strcmp('Aortic ascending flow',rois), 1);
    paf = find(strcmp('Pulmonary artery',rois), 1);
    if ~isempty(aaf) || ~isempty(paf)
      tots = reportflow('gettotal');
      fwds = reportflow('getfwdflow');
      bwds = reportflow('getbwdflow');
      if ~isempty(aaf)
        svaorta = tots(aaf);
        fwaorta = fwds(aaf);
        bwaorta = -bwds(aaf);
      end
      if ~isempty(paf)
        svpulmo = tots(paf);
        fwpulmo = fwds(paf);
        bwpulmo = -bwds(paf);
      end
    end
    reportflow('close_Callback');
  end
end

%---------------------
function result = qpqs
%---------------------
%Calculate Qp/Qs and display result in a message box

[svaorta,svpulmo] = getsv;
result = [];
if ~isempty(svaorta) && ~isempty(svpulmo)
  result = svpulmo/svaorta;
end

if nargout == 0
  if ~isempty(result)
    stri = sprintf([...
      'SV aorta: %0.1f ml\n'...
      'SV pulmonalis: %0.1f ml\n'...
      'Qp/Qs ratio: %0.0f %%'],...
      svaorta,svpulmo,100*result);
    msgbox(stri, 'Qp/Qs');
  else
    myfailed('Could not find flows for both aorta and pulmonalis');
  end
end

%----------------------------------------------------
function [mitrfrac,tricfrac,mitrvol,tricvol] = regurg %#ok<DEFNU>
%----------------------------------------------------
%Calculate regurgitant fractions for mitralis and tricusp
global SET
cineshortaxisno = findfunctions('findcineshortaxisno');
if isempty(cineshortaxisno) || SET(cineshortaxisno).SV == 0
  myfailed('SV from cine short-axis image stack is missing.');
  return;
end
[~,~,fwaorta,bwaorta,fwpulmo,bwpulmo] = getsv; %net volumes from flow analysis
stri = '';

if ~isempty(bwaorta)
  stri = [stri sprintf([...
    'Regurgitant volume aorta: %0.0f ml\n'...
    'Regurgitant fraction aorta: %0.0f %%\n'],...
    bwaorta,100*bwaorta/fwaorta)];
end

if ~isempty(bwpulmo)
  stri = [stri sprintf([...
    'Regurgitant volume pulmonary artery: %0.0f ml\n'...
    'Regurgitant fraction pulmonary artery: %0.0f %%\n'],...
    bwpulmo,100*bwpulmo/fwpulmo)];
end

lvsv = SET(cineshortaxisno).SV; %planometric stroke volume
if lvsv > 0 && ~isempty(fwaorta)
  mitrvol = lvsv - fwaorta;
  mitrfrac = mitrvol / lvsv;
  stri = [stri sprintf([...
    'Regurgitant volume mitralis: %0.0f ml\n'...
    'Regurgitant fraction mitralis: %0.0f %%\n'],...
    mitrvol,100*mitrfrac)];
else
  mitrvol = [];
  mitrfrac = [];
end

rvsv = SET(cineshortaxisno).RVSV;
if rvsv > 0 && ~isempty(fwpulmo)
  tricvol = rvsv - fwpulmo;
  tricfrac = tricvol / rvsv;
  stri = [stri sprintf([...
    'Regurgitant volume tricusp: %0.0f ml\n'...
    'Regurgitant fraction tricusp: %0.0f %%'],...
    tricvol,100*tricfrac)];
else
  tricvol = [];
  tricfrac = [];
end

if nargout == 0
  if ~isempty(stri)
    msgbox(stri, 'Shunt and Valve analysis');
  else
    myfailed('Could not find sufficient stroke volumes')
  end
end