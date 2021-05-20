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
eddycheck = false;
isinvisible = true;
for floop =  flowmagnitudeno
  if SET(floop).RoiN > 0
    reportflow('init',floop,eddycheck,isinvisible);
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
function result = qpqs %#ok<DEFNU>
%---------------------
%Calculate Qp/Qs and display result in a message box

[svaorta,svpulmo] = getsv;
result = [];
if ~isempty(svaorta) && ~isempty(svpulmo)
  result = svpulmo/svaorta;
end

if nargout == 0
  if ~isempty(result)
    stri = dprintf(['SV aorta: %0.1f ml\n'...
      'SV pulmonalis: %0.1f ml\n'...
      'Qp/Qs ratio: %0.0f %%'],...
      svaorta,svpulmo,100*result);
    mymsgbox(stri, 'Qp/Qs');
  else
    myfailed('Could not find flows for both aorta and pulmonalis');
  end
end

%----------------------------------------------------
function [mitrfrac,tricfrac,mitrvol,tricvol] = regurg %#ok<DEFNU>
%----------------------------------------------------
%Calculate regurgitant fractions for mitralis and tricusp
global SET
mitrvol = [];
mitrfrac = [];
tricvol = [];
tricfrac = [];
cineshortaxisno = findfunctions('findcineshortaxisno');
if isempty(cineshortaxisno) || SET(cineshortaxisno).SV == 0
  myfailed('SV from cine short-axis image stack is missing.');
  return;
end
[~,~,fwaorta,bwaorta,fwpulmo,bwpulmo] = getsv; %net volumes from flow analysis
stri = '';

if ~isempty(bwaorta)
  stri = [stri dprintf(['Regurgitant volume aorta: %0.0f ml\n'...
    'Regurgitant fraction aorta: %0.0f %%\n'],...
    bwaorta,100*bwaorta/fwaorta)];
end

if ~isempty(bwpulmo)
  stri = [stri dprintf(['Regurgitant volume pulmonary artery: %0.0f ml\n'...
    'Regurgitant fraction pulmonary artery: %0.0f %%\n'],...
    bwpulmo,100*bwpulmo/fwpulmo)];
end

lvsv = SET(cineshortaxisno).SV; %planometric stroke volume
if lvsv > 0 && ~isempty(fwaorta)
  mitrvol = lvsv - fwaorta;
  mitrfrac = mitrvol / lvsv;
  stri = [stri dprintf(['Regurgitant volume mitralis: %0.0f ml\n'...
    'Regurgitant fraction mitralis: %0.0f %%\n'],...
    mitrvol,100*mitrfrac)];
end

rvsv = SET(cineshortaxisno).RVSV;
if rvsv > 0 && ~isempty(fwpulmo)
  tricvol = rvsv - fwpulmo;
  tricfrac = tricvol / rvsv;
  stri = [stri dprintf(['Regurgitant volume tricusp: %0.0f ml\n'...
    'Regurgitant fraction tricusp: %0.0f %%'],...
    tricvol,100*tricfrac)];
end

if nargout == 0
  if ~isempty(stri)
    mymsgbox(stri, 'Shunt and Valve analysis');
  else
    myfailed('Could not find sufficient stroke volumes')
  end
end

%-----------------------------------------------------------------
function [lpafrac,rpafrac,svrpa,svlpa] = getlungflow %#ok<DEFNU> called from report2clipboard
%-----------------------------------------------------------------
%Calculate and return stroke volume and right/left pulmonary blood flow
%fractions
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
svrpa = [];
svlpa = [];
svpa = [];
eddycheck = false;
isinvisible = true;
for floop =  flowmagnitudeno
  if SET(floop).RoiN > 0
    reportflow('init',floop,eddycheck,isinvisible);
    rois = reportflow('getroiname');
    rpaind = find(strcmp('RPA',rois), 1);
    lpaind = find(strcmp('LPA',rois), 1);
    paind = find(strcmp('Pulmonary Artery',rois), 1);
    if ~isempty(rpaind) || ~isempty(lpaind)
      nettoflow = reportflow('gettotal'); % netto flow of all ROIs
      if ~isempty(rpaind)
        svrpa = nettoflow(rpaind);
      end
      if ~isempty(lpaind)
        svlpa = nettoflow(lpaind);
      end
      if ~isempty(paind)
        svpa = nettoflow(paind);
      end
    end
    reportflow('close_Callback');
  end
end
lpafrac = [];
rpafrac = [];
% calculate fraction for RPA and LPA depending whether both
% exist or not
if ~isempty(svrpa) && ~isempty(svlpa)
  svlung = svlpa + svrpa; 
  lpafrac = svlpa/svlung;
  rpafrac = svrpa/svlung;
else
  if ~isempty(svlpa) && ~isempty(svpa)
    lpafrac = svlpa/svpa;
  end
  if ~isempty(svrpa) && ~isempty(svpa)
    rpafrac = svrpa/svpa; 
  end
end


