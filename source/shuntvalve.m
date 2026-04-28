function varargout = shuntvalve(varargin)
%Functionality for Qp/Qs and Shunt and Valve analysis

%Nils Lundahl

%#ok<*GVMIS> 

if nargin == 0
  varargin = {'qpqs'};
end

[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard

%----------------------------------------------------------------------------
function [svaorta,svpulmo,fwaorta,bwaorta,fwpulmo,bwpulmo,noaorta,nopulmo] = getsv
%----------------------------------------------------------------------------
%Calculate and return stroke, forward and backward volumes both for the
%aorta and the pulmonary artery.

global SET
% init output
svaorta = [];
svpulmo = [];
fwaorta = [];
bwaorta = [];
fwpulmo = [];
bwpulmo = [];

[foundstacks,noaorta,nopulmo,roinumaorta,roinumpulmo] = findfunctions('findqpqsno');

if foundstacks(1)
  svaorta = SET(noaorta).Flow.Result(roinumaorta).nettotvol;
  fwaorta = SET(noaorta).Flow.Result(roinumaorta).netforwardvol;
  bwaorta = SET(noaorta).Flow.Result(roinumaorta).netbackwardvol;
  % #3238 report backward as positive value
  bwaorta = -bwaorta;
end

if foundstacks(2)
  svpulmo = SET(nopulmo).Flow.Result(roinumpulmo).nettotvol;
  fwpulmo = SET(nopulmo).Flow.Result(roinumpulmo).netforwardvol;
  bwpulmo = SET(nopulmo).Flow.Result(roinumpulmo).netbackwardvol;
  % #3238 report backward as positive value
  bwpulmo = -bwpulmo;
end

%-----------------------------------------------------------------
function [coaorta,copulm,noaorta,nopulmo] = getco
%-----------------------------------------------------------------
%Calculate and return CO for ROIs aorta and pulm.
coaorta = [];
copulm = [];

[foundstacks,noaorta,nopulmo,roinumaorta,roinumpulmo] = findfunctions('findqpqsno');
if all(foundstacks)
  coaorta = calcfunctions('calcflowco',noaorta,roinumaorta);
  copulm = calcfunctions('calcflowco',nopulmo,roinumpulmo);
end

%---------------------
function [qpqssv,qpqsco] = qpqs(noaorta,nopulmo,roinumaorta,roinumpulmo)
%---------------------
%Calculate Qp/Qs and display result in a message box

global SET
myworkon;
qpqssv = [];
qpqsco = [];

if nargin ~= 4
  [foundstacks,noaorta,nopulmo,roinumaorta,roinumpulmo] = findfunctions('findqpqsno');
else
  % input exist (as in call from 'reporter.reportgenerator('getqpqs')'
  foundstacks = true;
end

if all(foundstacks)
  % QpQs based on stroke volume SV
  svaorta = SET(noaorta).Flow.Result(roinumaorta).nettotvol;
  svpulmo = SET(nopulmo).Flow.Result(roinumpulmo).nettotvol;
  qpqssv = calcfunctions('calcqpqs',svpulmo,svaorta);

  % QpQs value based on cardiac outpur CO
  coaorta = calcfunctions('calcflowco',noaorta,roinumaorta);
  copulm = calcfunctions('calcflowco',nopulmo,roinumpulmo);
  qpqsco = calcfunctions('calcqpqs',copulm,coaorta);
end

myworkoff;
if nargout == 0
  if ~isempty(qpqssv)
    strstack = sprintf('\n%s #%d / #%d\n\n',dprintf('Stack'),nopulmo,noaorta);
    strisv = dprintf(['Net volume aorta: %0.1f ml\n'...
      'Net volume pulmonalis: %0.1f ml\n'...
      'Qp/Qs (%s) ratio: %0.2f'],svaorta,svpulmo,dprintf('SV'),qpqssv);

    stricoao = sprintf('\n\n%s: %0.1f l/min\n',dprintf('Aortic cardiac output'),coaorta);
    stricopul = sprintf('%s: %0.1f l/min\n',dprintf('Pulmonary cardiac output'),copulm);
    costr = dprintf('CO');
    striqpqsco = sprintf('%s: %0.2f',dprintf('Qp/Qs (%s) ratio',costr),qpqsco);
    stri = [strstack,strisv,stricoao,stricopul,striqpqsco];

    mymsgbox(stri, 'Qp/Qs');
    pause(0.1)
  else
    myfailed('Could not find net volumes for both aorta and pulmonalis');
  end
end
myworkoff;

%----------------------------------------------------
function [mitrfrac,tricfrac,mitrvol,tricvol] = regurg 
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
  stri = [stri sprintf('%s: %0.0f ml\n%s: %0.0f %%\n',...
          dprintf('Regurgitant volume aorta'),bwaorta,...
          dprintf('Regurgitant fraction aorta'),100*bwaorta/fwaorta)];
end

if ~isempty(bwpulmo)
   stri = [stri sprintf('%s: %0.0f ml\n%s: %0.0f %%\n',...
          dprintf('Regurgitant volume pulmonary artery'),bwpulmo,...
          dprintf('Regurgitant fraction pulmonary artery'),100*bwpulmo/fwpulmo)];
end

lvsv = SET(cineshortaxisno).SV; %planometric stroke volume
rvsv = SET(cineshortaxisno).RVSV;
[mitrvol,mitrfrac,tricvol,tricfrac] = calcfunctions('calcshuntvalve',lvsv,rvsv,fwaorta,fwpulmo);

if ~isempty(mitrvol) && ~isempty(mitrfrac)
  stri = [stri sprintf('%s: %0.0f ml\n%s: %0.0f %%\n',...
    dprintf('Regurgitant volume mitralis'),mitrvol,...
    dprintf('Regurgitant fraction mitralis'),mitrfrac)];
end

if ~isempty(tricvol) && ~isempty(tricfrac)
  stri = [stri sprintf('%s: %0.0f ml\n%s: %0.0f %%',...
    dprintf('Regurgitant volume tricusp'),tricvol,...
    dprintf('Regurgitant fraction tricusp'),tricfrac)];
end

if nargout == 0
  if ~isempty(stri)
    msgstr = sprintf('\n%s\n\n%s',dprintf('Shunt and Valve analysis'),stri);
    mymsgbox(msgstr, '');
  else
    myfailed('Could not find sufficient stroke volumes')
  end
end
myworkoff;

%-----------------------------------------------------------------
function [lpafrac,rpafrac,svrpa,svlpa] = getlungflow
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
    rpaind = find(strcmpi('RPA',rois), 1);
    lpaind = find(strcmpi('LPA',rois), 1);
    paind = find(strcmpi('Pulmonary Artery',rois), 1);
    if ~isempty(rpaind) || ~isempty(lpaind) || ~isempty(paind)
      nettoflow = reportflow('gettotal'); % netto flow of all ROIs
      if ~isempty(rpaind)
        svrpa = abs(nettoflow(rpaind));
      end
      if ~isempty(lpaind)
        svlpa = abs(nettoflow(lpaind));
      end
      if ~isempty(paind)
        svpa = abs(nettoflow(paind));
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
  rpafrac = svrpa/svlung;
  lpafrac = 1-rpafrac;  
else
  if ~isempty(svlpa) && ~isempty(svpa)
    lpafrac = svlpa/svpa;
    rpafrac = 1-lpafrac;
  end
  if ~isempty(svrpa) && ~isempty(svpa)
    rpafrac = svrpa/svpa; 
    lpafrac = 1-rpafrac;
  end
end


