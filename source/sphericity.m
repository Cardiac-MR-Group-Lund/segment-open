function sphericity
%------------------
%Tool for calculating sphericity of the left ventricle

global SET

cineno = findfunctions('findno');
if isempty(cineno)
  myfailed('No cine stack found')
  return
end
saxno = findfunctions('findcineshortaxisno');
if ~existfunctions('existendo',saxno)
  myfailed('No LV endocardium available in short-axis cine stack.')
  return
end
cineno = setdiff(cineno,saxno);
if isempty(cineno)
  myfailed('Could not find additional cine stack except the short-axis cine stack');
  return;
end

edlviewplanes = {};
edlmeas = [];
eslviewplanes = {};
eslmeas = [];
for noloop = cineno
  imvp = SET(noloop).ImageViewPlane;
  measures = [SET(noloop).Measure];
  for i = 1:numel(measures)
    if strcmp(measures(i).Name,'EDL')
      if ~ismember(imvp,edlviewplanes)
        edlmeas = [edlmeas measures(i)];
        edlviewplanes = [edlviewplanes imvp];
      else
        mywarning(['Found multiple measurements of EDL in the same ' ...
          'image stack. Taking first (arbitrary decision).'])
      end
    elseif strcmp(measures(i).Name,'ESL')
      if ~ismember(imvp,eslviewplanes)
        eslmeas = [eslmeas measures(i)];
        eslviewplanes = [eslviewplanes imvp];
      else
        mywarning(['Found multiple measurements of ESL in the same ' ...
          'image stack. Taking first (arbitrary decision).'])
      end
    end
  end
end

if isempty(edlmeas) || isempty(eslmeas)
  myfailed('Could not find EDL or ESL measures in any cine stack')
  return
end

if SET(saxno).EDT == SET(saxno).EST
  mywarning('EDT occurs at the same time as EST in short-axis cine stack');
end

edmax = maxsaxdiameter(saxno,SET(saxno).EDT,'LV');
esmax = maxsaxdiameter(saxno,SET(saxno).EST,'LV');

c = cell(2,5+numel(edlviewplanes) + numel(eslviewplanes));
c{1,1} = 'Patient name';
c{2,1} = SET(1).PatientInfo.Name;
c{1,2} = 'SAX diameter ED [mm]';
c{2,2} = sprintf('%0.1f',edmax);
c{1,3} = 'SAX diameter ES [mm]';
c{2,3} = sprintf('%0.1f',esmax);
c{1,4} = 'Sphericity ED';
if ~isempty(edlmeas)
  c{2,4} = sprintf('%0.2f',edmax/max([edlmeas.Length]));
else
  c{2,4} = '-';
end
c{1,5} = 'Sphericity ES';
if ~isempty(eslmeas)
  c{2,5} = sprintf('%0.2f',esmax/max([eslmeas.Length]));
else
  c{2,5} = '-';
end
for i = 1:numel(edlviewplanes)
  c{1,5+i} = sprintf('EDL (%s) [mm]',edlviewplanes{i});
  c{2,5+i} = sprintf('%0.1f',edlmeas(i).Length);
end
for j = 1:numel(eslviewplanes)
  c{1,5+i+j} = sprintf('ESL (%s) [mm]',eslviewplanes{j});
  c{2,5+i+j} = sprintf('%0.1f',eslmeas(j).Length);
end

stri = '';
for loop = 2:size(c,2)
  stri = [stri sprintf('%s: %s\n',c{1,loop},c{2,loop})];
end
mymsgbox(stri,'LV Sphericity')
segment('cell2clipboard',c);
%------------------------------------------
function endomax = maxsaxdiameter(no,tf,type)
%------------------------------------------
%Get maximum endocardial diameter in short-axis cine stack
% no - no of stack
% tf - time frame
% type - 'RV', 'LV'
global SET
xres = SET(no).ResolutionX;
yres = SET(no).ResolutionY;

switch type
  case 'LV'
    edxall = squeeze(SET(no).EndoX(:,tf,:));
    edyall = squeeze(SET(no).EndoY(:,tf,:));
  case 'RV'
    edxall = squeeze(SET(no).RVEndoX(:,tf,:));
    edyall = squeeze(SET(no).RVEndoY(:,tf,:));
end
endomax = 0;

for z = 1:SET(no).ZSize
  edx = edxall(:,z);
  if ~isnan(edx(1))
    edy = edyall(:,z);
    edxmat = repmat(edx,1,length(edx));
    edymat = repmat(edy,1,length(edy));
    eddist = sqrt((xres*(edxmat'-edxmat)).^2+ ...
      (yres*(edymat'-edymat)).^2);
    endomax = max(endomax,max(eddist(:)));
  end
end

