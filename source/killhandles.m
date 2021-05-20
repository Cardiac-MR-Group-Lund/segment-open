function mainhandles = killhandles(mainfig,subfigfiles)
%------------------------------------------------------
%Generates a structure of handles, copied from mainfig but where those 
%not present in subfigfile are set to [].

if ~iscell(subfigfiles)
  subfigfiles = {subfigfiles};
end
mainhandles = guihandles(mainfig);

for i = 1:numel(subfigfiles)
  subfigfile = subfigfiles{i};
  subfig = openfig(subfigfile,'new','invisible');
  subhandles = guihandles(subfig);
  
  fnames = fieldnames(subhandles);
  addhandles = addedhandles(subfigfile);
  fnames = [fnames; addhandles];
  
  for loop = 1:length(fnames)
    field = fnames{loop};
    if ~isfield(mainhandles,field)
      mainhandles.(field) = [];
    end
  end
  
  delete(subfig);
end

%---------------------------------------------
function newhandles = addedhandles(subfigfile)
%---------------------------------------------
newhandles=[];

switch subfigfile
	case 'segment.fig'
		newhandles = {
			'perfusionpushbutton'
			't2starpushbutton'
			'lvmenu'
			'rvmenu'
			'marmenu'
			'analysismenu'
      'filesavetopacsmenu'
      'filesavesegdicom'
		};
	case 'segpref.fig'
		
	case 'openfile.fig'
end