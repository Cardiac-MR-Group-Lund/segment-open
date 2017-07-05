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
  if strcmp(subfigfile,'segment.fig');
    icons1 = maintoolbaricons;
    icons2 = segment('initviewtoolbar');
    fnames = [fnames; icons1; icons2];
  end
  
  addhandles = addedhandles(subfigfile);
  fnames = [fnames; addhandles];
  
  for loop = 1:length(fnames);
    field = fnames{loop};
    if ~isfield(mainhandles,field)
      mainhandles.(field) = [];
    end
  end
  
  delete(subfig);
end

%--------------------------------
function icons = maintoolbaricons
%--------------------------------
icons = {'maintoolbar'
        'fileopenicon'
        'filesaveicon'
        'databaseicon'
        'databaseaddicon'
        'pacsicon'
        'pacsaddicon'
        'importfromcdicon'
        'view1panelicon'
        'view2panelicon'
        'view2x1panelicon'
        'view3panelicon'
        'view1x3panelicon'
        'view4panelicon'
        'view6panelicon'
        'view9panelicon'
        'view12panelicon'
        'view16panelicon'
        'panelallicon'
        'saveviewicon'
        'hidepinsicon'
        'hideothercontouricon'
        'hideinterpicon'
        'hidelvicon'
        'hidervicon'
        'hidescaricon'
        'hidemaricon'
        'hideroiicon'
        'hidemeasuresicon'
        'hidepointsicon'
        'hideplusicon'
        'hideintersectionsicon'
        'hidetexticon'
        'hidepapicon'
        'hidesectorgridicon'
        'hideoverlayicon'
        'abouticon'
        'viewpixelyicon'
        'orthoviewicon'
        'mipicon'
        'colorbaricon'
        'autocontrastallicon'};

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