function isopen=isopengui(figname)
%checks if an gui is already open
%input argument figurename, for example 'levelset.fig'

allfigname=get(allchild(0),'filename');

isopen=0;
if iscell(allfigname)
	num=length(allfigname);
	for j=1:num
		[path,name,ext]=fileparts(allfigname{j});
		thisfigname=[name ext];
		if isequal(thisfigname,figname)
			isopen=1;
			return;
		end
	end
else
  [path,name,ext]=fileparts(allfigname);
  thisfigname=[name ext];
	if isequal(thisfigname,figname)
		isopen=1;
		return
	end
end