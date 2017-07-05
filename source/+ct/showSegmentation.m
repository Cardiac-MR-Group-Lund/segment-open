function showSegmentation(call,arg1,arg2)

global setCell helpVar handles 

if nargin==0
	%open GUI
	figName=['+ct' filesep 'showSegmentation.fig'];
	fig=openfig(figName,'reuse');
	handles=guihandles(fig);
	handles.fig=fig;
	
	set(fig,'renderer','zbuffer')
	
	%Set keypress for entire GUI
  segment('recursekeypressfcn',handles.fig,@showSegmentationkeypressed);
	
	helpVar.XSize=setCell{helpVar.showLevel}.size(1);
	helpVar.YSize=setCell{helpVar.showLevel}.size(2);
	helpVar.ZSize=setCell{helpVar.showLevel}.size(3);
	
	helpVar.rslice = round(helpVar.ZSize/2);
	helpVar.gslice = round(helpVar.YSize/2);
	helpVar.bslice = round(helpVar.XSize/2);

	warning('off','MATLAB:contour:ConstantData'); %Ignore warning messages from contour.

	%Install icons
  try
    load('icons.mat','icon');
  catch
    myfailed('Critical error: Could not read icons.');
    return;
  end;
	
	mytoolbar = uitoolbar(handles.fig);
	
	%Selection icon
  props.ClickedCallback = 'ct.showSegmentation(''selecttool'')';
  props.ToolTip = 'Pointer tool.';
  props.CData = icon.selectslices;
  props.Tag = 'selecticon';
  props.separator = 'off';
  handles.selecticon = uitoggletool(mytoolbar,props);
	
	%zoom in icon
  props.ClickedCallback = 'ct.showSegmentation(''zoomin'')';
  props.ToolTip = 'Zoom in';
  props.CData = icon.zoomin;
  props.Tag = 'zoominicon';
  props.separator = 'off';
  handles.zoominicon = uitoggletool(mytoolbar,props);

  %zoom out icon
  props.ClickedCallback = 'ct.showSegmentation(''zoomout'')';
  props.ToolTip = 'Zoom out';
  props.CData = icon.zoomout;
  props.Tag = 'zoomouticon';
  props.separator = 'off';
  handles.zoomouticon = uitoggletool(mytoolbar,props);

  %refresh
  props.ClickedCallback = 'ct.showSegmentation(''refresh'')';
  props.ToolTip = 'Refresh image';
  props.CData = icon.refresh;
  props.Tag = 'refreshicon';
  props.separator = 'on';
  handles.refreshicon = uitoggletool(mytoolbar,props);

	
	%set default values
	showSegmentationDefault
	%set buttons default
	ct.showSegmentation('updatebuttons','default');	
	%set sliders
	ct.showSegmentation('updatesliders');
	%update image
	ct.showSegmentation('update');
	
else
	switch call
		case 'updatebuttons'
			switch arg1
				case 'CTWholeHeartSegmentation'
					set(handles.wholeheartpushbutton,'visible', 'on', 'enable', 'on');
					set(handles.coronaryostiapushbutton,'visible','off','enable', 'off');
					set(handles.bonesegmentationpushbutton,'visible','off','enable','off');
					set(handles.runoffpushbutton,'visible','off','enable','off');
					set(handles.carotidrenalpushbutton,'visible','off','enable','off');
					set(handles.bonethreshslider,'visible','off','enable','off');
					set(handles.bonethreshtext,'visible','off');
				case 'coronaryOstiaSearch'
					switch arg2
						case 'init'
							set(handles.wholeheartpushbutton,'visible', 'on', 'enable', 'off');
							set(handles.coronaryostiapushbutton,'visible','on','enable', 'off');
						case 'segmentation done'
							set(handles.wholeheartpushbutton,'visible', 'on', 'enable', 'on');
							set(handles.coronaryostiapushbutton,'visible','on','enable', 'on');
					end
					set(handles.bonesegmentationpushbutton,'visible','off','enable','off');
					set(handles.runoffpushbutton,'visible','off','enable','off');
					set(handles.carotidrenalpushbutton,'visible','off','enable','off');
					set(handles.bonethreshslider,'visible','off','enable','off');
					set(handles.bonethreshtext,'visible','off');
				case 'CTBoneSegmentation'
					set(handles.wholeheartpushbutton,'visible', 'off', 'enable', 'off');
					set(handles.coronaryostiapushbutton,'visible','off','enable', 'off');
					set(handles.bonesegmentationpushbutton,'visible','on','enable','on');
					set(handles.runoffpushbutton,'visible','off','enable','off');
					set(handles.carotidrenalpushbutton,'visible','off','enable','off');
					set(handles.bonethreshslider,'visible','on','enable','on');
					set(handles.bonethreshtext,'visible','on');
					set(handles.bonethreshslider,'min',-0.5,'max',1.5,'value',0.5,...
						'sliderstep',[0.05 0.1]);
				case 'CTARunOff'
					set(handles.wholeheartpushbutton,'visible', 'off', 'enable', 'off');
					set(handles.coronaryostiapushbutton,'visible','off','enable', 'off');
					set(handles.bonesegmentationpushbutton,'visible','off','enable','off');
					set(handles.runoffpushbutton,'visible','on','enable','on');
					set(handles.carotidrenalpushbutton,'visible','off','enable','off');
					set(handles.bonethreshslider,'visible','off','enable','off');
					set(handles.bonethreshtext,'visible','off');
				case 'CTACarotidRenal'
					set(handles.wholeheartpushbutton,'visible', 'off', 'enable', 'off');
					set(handles.coronaryostiapushbutton,'visible','off','enable', 'off');
					set(handles.bonesegmentationpushbutton,'visible','off','enable','off');
					set(handles.runoffpushbutton,'visible','off','enable','off');
					set(handles.carotidrenalpushbutton,'visible','on','enable','on');
					set(handles.bonethreshslider,'visible','off','enable','off');
					set(handles.bonethreshtext,'visible','off');
				case 'default'
					ct.showSegmentation('updatebuttons','CTWholeHeartSegmentation');
				otherwise
					set(handles.wholeheartpushbutton,'visible', 'on', 'enable', 'off');
					set(handles.coronaryostiapushbutton,'visible','on','enable', 'off');
					set(handles.bonesegmentationpushbutton,'visible','on','enable','off');
					set(handles.runoffpushbutton,'visible','on','enable','off');
					set(handles.carotidrenalpushbutton,'visible','on','enable','off');
					set(handles.bonethreshslider,'visible','off','enable','off');
					set(handles.bonethreshtext,'visible','off');
			end
		case 'updatesliders'
      set(handles.rslider,'min',1,'max',helpVar.ZSize,'value',helpVar.ZSize-helpVar.rslice+1,...
				'sliderstep',[1/helpVar.ZSize 3/helpVar.ZSize]);
			set(handles.gslider,'min',1,'max',helpVar.YSize,'value',helpVar.gslice,...
				'sliderstep',[1/helpVar.YSize 3/helpVar.YSize]);
			set(handles.bslider,'min',1,'max',helpVar.XSize,'value',helpVar.bslice,...
				'sliderstep',[1/helpVar.XSize 3/helpVar.XSize]);
		case 'slice'
			switch arg1
				case 'r'
					helpVar.rslice = max(min(helpVar.rslice-arg2,helpVar.ZSize),1);
				case 'g'
					helpVar.gslice = max(min(helpVar.gslice+arg2,helpVar.YSize),1);
				case 'b'
					helpVar.bslice = max(min(helpVar.bslice+arg2,helpVar.XSize),1);		
			end;
			ct.showSegmentation('update');
		case 'rslider'
			helpVar.rslice = 1+helpVar.ZSize-min(helpVar.ZSize,max(round(mygetslider(handles.rslider)),1));
			set(handles.rslider,'Value',helpVar.ZSize-helpVar.rslice+1);
			ct.showSegmentation('rdraw');
			ct.showSegmentation('sdraw');
			ct.showSegmentation('lineupdate');
		case 'gslider'
			helpVar.gslice = min(helpVar.YSize,max(round(mygetslider(handles.gslider)),1));
			set(handles.gslider,'Value',helpVar.gslice);
			ct.showSegmentation('gdraw');
			ct.showSegmentation('lineupdate');
		case 'bslider'
			helpVar.bslice = min(helpVar.XSize,max(round(mygetslider(handles.bslider)),1));
			set(handles.bslider,'Value',helpVar.bslice);
			ct.showSegmentation('bdraw');
			ct.showSegmentation('lineupdate');
		case 'bonethreshslider'
			thresh=mygetslider(handles.bonethreshslider);
			ct.CTBoneSegmentation('adjustthreshhold',thresh);
		case 'overviewzoom'
			figure(12);
			handles.overviewaxes = gca;
			ct.showSegmentation('update');
			ct.showSegmentation('isosurface');
			axis off image;
			cameratoolbar(12);
		case 'isosurface'
			if isempty(setCell{helpVar.showLevel}.segmentation)
				myfailed('No segmentation exists');
				return;
			end
			fv = isosurface(double(...
				setCell{helpVar.showLevel}.segmentation),...
				0.5);

			if get(handles.reducecheckbox,'value')
				fv = reducepatch(fv,0.2);
			end;

			%Flip in z-direction and xy
			fv.vertices(:,3) = (helpVar.ZSize+1)-fv.vertices(:,3);
			temp = fv.vertices(:,2);
			fv.vertices(:,2) = fv.vertices(:,1);
			fv.vertices(:,1) = temp;

			%--- Display
			hold(handles.overviewaxes,'on');
			h = patch(fv,'parent',handles.overviewaxes);
			hold(handles.overviewaxes,'off');
			set(h,'facecolor',[1 0 0],'edgecolor',[1 1 0],'facealpha',0.5);
			if not(get(handles.polygonscheckbox,'value'))
				set(h,'edgealpha',0);
			end;
		case 'update'			
			ct.showSegmentation('rdraw');
			ct.showSegmentation('gdraw');
			ct.showSegmentation('bdraw');

			ct.showSegmentation('lineupdate');
			ct.showSegmentation('updatesliders');
		case 'selection'
			helpVar.ViewSelection = get(handles.selectioncheckbox,'value');
			ct.showSegmentation('update');
		case 'outline'
			helpVar.ViewOutline = get(handles.outlinecheckbox,'value');
			ct.showSegmentation('update');
		case 'mip'
			helpVar.ViewMIP = get(handles.mipcheckbox,'value');
			ct.showSegmentation('update');
		case 'lineupdate'
			set(handles.grline,'ydata',[helpVar.rslice helpVar.rslice]);
			set(handles.brline,'ydata',[helpVar.rslice helpVar.rslice]);

			set(handles.rgline,'xdata',[helpVar.gslice helpVar.gslice]);
			set(handles.bgline,'xdata',[helpVar.gslice helpVar.gslice]);
			
			set(handles.rbline,'ydata',[helpVar.bslice helpVar.bslice]);
			set(handles.gbline,'xdata',[helpVar.bslice helpVar.bslice]);
		case 'rdraw'
			%--- Red image
			handles.rimage = image(squeeze(setCell{helpVar.showLevel}.IM(:,:,helpVar.rslice)),'parent',handles.raxes);

			%Add objects
			hold(handles.raxes,'on');
			handles.rbline = plot(handles.raxes,[0 helpVar.YSize],[helpVar.bslice helpVar.bslice],'b-');
			handles.rgline = plot(handles.raxes,[helpVar.gslice helpVar.gslice],[0 helpVar.XSize],'g-');
			x1 = helpVar.RZoomState(1)+0.5;
			x2 = helpVar.RZoomState(2)-0.5;
			y1 = helpVar.RZoomState(3)+0.5;
			y2 = helpVar.RZoomState(4)-0.5;
			handles.rimagebox = plot(handles.raxes,[x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'r-');
			set(handles.rimagebox,'linewidth',2);

			if ~isempty(setCell{helpVar.showLevel}.segmentation)
				[c,h] = contour(...
					handles.raxes,...
					double(squeeze(setCell{helpVar.showLevel}.segmentation(:,:,helpVar.rslice))),...
					[0.5 0.5]);
				set(h,'ButtonDownFcn','ct.showSegmentation(''click'',''r'')');
				children = get(h,'children');
				set(children,'edgecolor',[1 1 0]);
				if helpVar.ViewOutline
					set(children,'visible','on');
				else
					set(children,'visible','off');
				end;
			end

			
			%show landmarks
			if ~isempty(setCell{helpVar.showLevel}.landmarks.X)
				for loop=1:length(setCell{helpVar.showLevel}.landmarks.X)
					if round(setCell{helpVar.showLevel}.landmarks.Z(loop))==helpVar.rslice
						handles.pointp(loop) = plot(handles.raxes,...
							setCell{helpVar.showLevel}.landmarks.Y(loop),setCell{helpVar.showLevel}.landmarks.X(loop),'w+');
						handles.pointo(loop) = plot(handles.raxes,...
							setCell{helpVar.showLevel}.landmarks.Y(loop),setCell{helpVar.showLevel}.landmarks.X(loop),'wo');
						handles.pointtext(loop) = text(...
							'position',[setCell{helpVar.showLevel}.landmarks.Y(loop)+2 setCell{helpVar.showLevel}.landmarks.X(loop)],...
							'string',setCell{helpVar.showLevel}.landmarks.Label{loop},...
							'parent',handles.raxes,...
							'color',[1 1 1]);
					end
				end
			end

			hold(handles.raxes,'off');

			%Aspect ratio
			axis(handles.raxes,'off');
			set(handles.raxes,'Clim',[0 1],'dataaspectratio',...
				[1/setCell{helpVar.showLevel}.resolution(2) ...
				1/setCell{helpVar.showLevel}.resolution(1) 1]);

			set(handles.rimage,'ButtonDownFcn','ct.showSegmentation(''click'',''r'')');
			set([handles.rbline handles.rgline],'ButtonDownFcn','ct.showSegmentation(''click'',''r'')');
			ct.showSegmentation('rupdate');
		
		case 'gdraw'
			%--- Green image
			temp = squeeze(setCell{helpVar.showLevel}.IM(:,helpVar.gslice,:))';
			handles.gimage = image(temp,'parent',handles.gaxes);

			%Add objects
			hold(handles.gaxes,'on');

			handles.grline = plot(handles.gaxes,[0 helpVar.YSize],[helpVar.rslice helpVar.rslice],'r-');
			handles.gbline = plot(handles.gaxes,[helpVar.bslice helpVar.bslice],[0 helpVar.ZSize],'b-');

			x1 = helpVar.GZoomState(1)+0.5;
			x2 = helpVar.GZoomState(2)-0.5;
			y1 = helpVar.GZoomState(3)+0.5;
			y2 = helpVar.GZoomState(4)-0.5;
			handles.gimagebox = plot(handles.gaxes,[x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'g-');
			set(handles.gimagebox,'linewidth',2);

			set(handles.gaxes,'Clim',[0 1],'dataaspectratio',...
				[1/setCell{helpVar.showLevel}.resolution(2) ...
				1/setCell{helpVar.showLevel}.resolution(3) 1]);
			
			if ~isempty(setCell{helpVar.showLevel}.segmentation)
				%show segmentation
				[c,h] = contour(...
					handles.gaxes,...
					double(squeeze(setCell{helpVar.showLevel}.segmentation(:,helpVar.gslice,:)))',...
					[0.5 0.5]);
				set(h,'ButtonDownFcn','ct.showSegmentation(''click'',''g'')');
				children = get(h,'children');
				set(children,'edgecolor',[1 1 0]);
				if helpVar.ViewOutline
					set(children,'visible','on');
				else
					set(children,'visible','off');
				end;
			end

			
			%show landmarks
			if ~isempty(setCell{helpVar.showLevel}.landmarks.X)
					for loop=1:length(setCell{helpVar.showLevel}.landmarks.X)
						if round(setCell{helpVar.showLevel}.landmarks.Y(loop))==helpVar.gslice
							handles.pointp(loop) = plot(handles.gaxes,...
								setCell{helpVar.showLevel}.landmarks.X(loop),setCell{helpVar.showLevel}.landmarks.Z(loop),'w+');
							handles.pointo(loop) = plot(handles.gaxes,...
								setCell{helpVar.showLevel}.landmarks.X(loop),setCell{helpVar.showLevel}.landmarks.Z(loop),'wo');
							handles.pointtext(loop) = text(...
								'position',[setCell{helpVar.showLevel}.landmarks.X(loop)+2 setCell{helpVar.showLevel}.landmarks.Z(loop)],...
								'string',setCell{helpVar.showLevel}.landmarks.Label{loop},...
								'parent',handles.gaxes,...
								'color',[1 1 1]);
						end
					end
			end

			hold(handles.gaxes,'off');
			axis(handles.gaxes,'off');
			set(handles.gimage,'ButtonDownFcn','ct.showSegmentation(''click'',''g'')');
			set([handles.grline handles.gbline],'ButtonDownFcn','ct.showSegmentation(''click'',''g'')');
			
			set(handles.gaxes,'Clim',[0 1],'dataaspectratio',...
				[1/setCell{helpVar.showLevel}.resolution(2) ...
				1/setCell{helpVar.showLevel}.resolution(3) 1]);
			ct.showSegmentation('gupdate');
		case 'bdraw'
			
			%--- Blue image
			temp = squeeze(setCell{helpVar.showLevel}.IM(helpVar.bslice,:,:))';
			handles.bimage = image(temp,'parent',handles.baxes);
			
			%Add objects
			hold(handles.baxes,'on');
			handles.brline = plot(handles.baxes,[0 helpVar.XSize],[helpVar.rslice helpVar.rslice],'r-');
			handles.bgline = plot(handles.baxes,[helpVar.gslice helpVar.gslice],[0 helpVar.ZSize],'g-');
			x1 = helpVar.BZoomState(1)+0.5;
			x2 = helpVar.BZoomState(2)-0.5;
			y1 = helpVar.BZoomState(3)+0.5;
			y2 = helpVar.BZoomState(4)-0.5;
			handles.bimagebox = plot(handles.baxes,[x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'b-');
			set(handles.bimagebox,'linewidth',2);
			
			%Add objects
			hold(handles.raxes,'on');
			handles.rbline = plot(handles.raxes,[0 helpVar.YSize],[helpVar.bslice helpVar.bslice],'b-');
			handles.rgline = plot(handles.raxes,[helpVar.gslice helpVar.gslice],[0 helpVar.XSize],'g-');
			x1 = helpVar.RZoomState(1)+0.5;
			x2 = helpVar.RZoomState(2)-0.5;
			y1 = helpVar.RZoomState(3)+0.5;
			y2 = helpVar.RZoomState(4)-0.5;
			handles.rimagebox = plot(handles.raxes,[x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],'r-');
			set(handles.rimagebox,'linewidth',2);

			%show segmentation
			if ~isempty(setCell{helpVar.showLevel}.segmentation)
				[c,h] = contour(...
					handles.baxes,...
					double(squeeze(setCell{helpVar.showLevel}.segmentation(helpVar.bslice,:,:)))',...
					[0.5 0.5]);
				set(h,'ButtonDownFcn','ct.showSegmentation(''click'',''b'')');
				children = get(h,'children');
				set(children,'edgecolor',[1 1 0]);
				if helpVar.ViewOutline
					set(children,'visible','on');
				else
					set(children,'visible','off');
				end;
			end
			
			%show landmarks
			if ~isempty(setCell{helpVar.showLevel}.landmarks.X)
				for loop=1:length(setCell{helpVar.showLevel}.landmarks.X)
					if round(setCell{helpVar.showLevel}.landmarks.X(loop))==helpVar.bslice
						handles.pointp(loop) = plot(handles.baxes,...
							setCell{helpVar.showLevel}.landmarks.Y(loop),setCell{helpVar.showLevel}.landmarks.Z(loop),'w+');
						handles.pointo(loop) = plot(handles.baxes,...
							setCell{helpVar.showLevel}.landmarks.Y(loop),setCell{helpVar.showLevel}.landmarks.Z(loop),'wo');
						handles.pointtext(loop) = text(...
							'position',[setCell{helpVar.showLevel}.landmarks.Y(loop)+2 setCell{helpVar.showLevel}.landmarks.Z(loop)],...
							'string',setCell{helpVar.showLevel}.landmarks.Label{loop},...
							'parent',handles.baxes,...
							'color',[1 1 1]);
					end
				end
			end

			hold(handles.baxes,'off');
			axis(handles.baxes,'off');
			set(handles.bimage,'ButtonDownFcn','ct.showSegmentation(''click'',''b'')');
			set([handles.brline handles.bgline],'ButtonDownFcn','ct.showSegmentation(''click'',''b'')');

			set(handles.baxes,'Clim',[0 1],'dataaspectratio',...
				[1/setCell{helpVar.showLevel}.resolution(1) ...
				1/setCell{helpVar.showLevel}.resolution(3) 1]);
			ct.showSegmentation('bupdate');

		case 'rupdate'
			%--- Update 'red' image
			if helpVar.ViewMIP
				%MIP
				temp = squeeze(max(setCell{helpVar.showLevel}.IM,[],3));
				if ~isempty(setCell{helpVar.showLevel}.segmentation)
					overlaySegmentation  = squeeze(max(setCell{helpVar.showLevel}.segmentation,[],3))>0.5;
				else 
					overlaySegmentation=[];
				end
			else
				%Slice
				temp = squeeze(setCell{helpVar.showLevel}.IM(:,:,helpVar.rslice));
				if ~isempty(setCell{helpVar.showLevel}.segmentation)
					overlaySegmentation  = squeeze(setCell{helpVar.showLevel}.segmentation(:,:,helpVar.rslice))>0.5;
				else
					overlaySegmentation=[];
				end
			end;
			%Update image

			set(handles.rimage,'cdata',...
				remapandoverlaySegmentation(temp,overlaySegmentation));

			%zoom state
			if isempty(helpVar.RZoomState)
				helpVar.RZoomState = [0.5;size(temp,2)-0.5;0.5;size(temp,1)-0.5];
			end;
			set(handles.raxes,...
				'xlim',helpVar.RZoomState(1:2),...
				'ylim',helpVar.RZoomState(3:4));

			x1 = helpVar.RZoomState(1);
			x2 = helpVar.RZoomState(2);
			y1 = helpVar.RZoomState(3);
			y2 = helpVar.RZoomState(4);
			set(handles.rimagebox,'xdata',[x1 x2 x2 x1 x1],'ydata',[y1 y1 y2 y2 y1]);

		case 'gupdate'
			%--- Green image
			if helpVar.ViewMIP
				%MIP
				temp = squeeze(max(setCell{helpVar.showLevel}.IM,[],2))';
				if ~isempty(setCell{helpVar.showLevel}.segmentation)
					overlaySegmentation  = squeeze(max(setCell{helpVar.showLevel}.segmentation,[],2))'>0.5;
				else
					overlaySegmentation=[];
				end
			else
				%Slice
				temp = squeeze(setCell{helpVar.showLevel}.IM(:,helpVar.gslice,:))';	
				if ~isempty(setCell{helpVar.showLevel}.segmentation)
					overlaySegmentation  = squeeze(setCell{helpVar.showLevel}.segmentation(:,helpVar.gslice,:))'>0.5;
				else
					overlaySegmentation=[];
				end
			end;
			%Update image
			set(handles.gimage,'cdata',...
				remapandoverlaySegmentation(temp,overlaySegmentation));

			%zoom state
			if isempty(helpVar.GZoomState)
				helpVar.GZoomState = [0.5;size(temp,2)-0.5;0.5;size(temp,1)-0.5];
			end;
			set(handles.gaxes,...
				'xlim',helpVar.GZoomState(1:2),...
				'ylim',helpVar.GZoomState(3:4));

			x1 = helpVar.GZoomState(1);
			x2 = helpVar.GZoomState(2);
			y1 = helpVar.GZoomState(3);
			y2 = helpVar.GZoomState(4);
			set(handles.gimagebox,'xdata',[x1 x2 x2 x1 x1],'ydata',[y1 y1 y2 y2 y1]);

		case 'bupdate'
			if helpVar.ViewMIP
				%MIP
				temp = squeeze(max(setCell{helpVar.showLevel}.IM,[],1))';
				if ~isempty(setCell{helpVar.showLevel}.segmentation)
					overlaySegmentation = squeeze(max(setCell{helpVar.showLevel}.segmentation,[],1))'>0.5;
				else
					overlaySegmentation=[];
				end
			else
				%Slice
				temp = squeeze(setCell{helpVar.showLevel}.IM(helpVar.bslice,:,:))';
				if ~isempty(setCell{helpVar.showLevel}.segmentation)
					overlaySegmentation  = squeeze(setCell{helpVar.showLevel}.segmentation(helpVar.bslice,:,:))'>0.5;
				else
					overlaySegmentation=[];
				end
			end;
			%Update image
			set(handles.bimage,'cdata',...
				remapandoverlaySegmentation(temp,overlaySegmentation));
			%zoom state
			if isempty(helpVar.BZoomState)
				helpVar.BZoomState = [0.5;size(temp,2)-0.5;0.5;size(temp,1)-0.5];
			end;
			set(handles.baxes,...
				'xlim',helpVar.BZoomState(1:2),...
				'ylim',helpVar.BZoomState(3:4));

			x1 = helpVar.BZoomState(1);
			x2 = helpVar.BZoomState(2);
			y1 = helpVar.BZoomState(3);
			y2 = helpVar.BZoomState(4);
			set(handles.bimagebox,'xdata',[x1 x2 x2 x1 x1],'ydata',[y1 y1 y2 y2 y1]);

		case 'volume'
			if isempty(setCell{helpVar.showLevel}.segmentation)
				myfailed('No segmentation exists');
				return
			end
			num=sum(sum(sum(setCell{helpVar.showLevel}.segmentation>0.5)));
			voxel = cumprod(setCell{helpVar.showLevel}.resolution);
			num = num*voxel(end);
			num = num/1000; %convert to ml == cm^3
			
			msgbox(dprintf('Totalvolume of object is %0.5g [ml]',num));
		case 'zoomin'
			set(handles.zoominicon,'state','off');
			helpVar.Zoom = helpVar.Zoom*1.2;
			helpVar.RZoomState = zoomhelper(helpVar.RZoomState,1.2,...
				[helpVar.YSize helpVar.XSize]);
			helpVar.GZoomState = zoomhelper(helpVar.GZoomState,1.2,...
				[helpVar.XSize helpVar.ZSize]);
			helpVar.BZoomState = zoomhelper(helpVar.BZoomState,1.2,...
				[helpVar.YSize helpVar.ZSize]);
			
			ct.showSegmentation('rupdate');
			ct.showSegmentation('gupdate');
			ct.showSegmentation('bupdate');
		case 'zoomout'
			set(handles.zoomouticon,'state','off');
			helpVar.Zoom = helpVar.Zoom/1.2;
			helpVar.Zoom = helpVar.Zoom*1.2;
			helpVar.RZoomState = zoomhelper(helpVar.RZoomState,1/1.2,[helpVar.YSize helpVar.XSize]);
			helpVar.GZoomState = zoomhelper(helpVar.GZoomState,1/1.2,[helpVar.XSize helpVar.ZSize]);
			helpVar.BZoomState = zoomhelper(helpVar.BZoomState,1/1.2,[helpVar.YSize helpVar.ZSize]);

			ct.showSegmentation('rupdate');
			ct.showSegmentation('gupdate');
			ct.showSegmentation('bupdate');
		case 'selecttool'
			helpVar.CurrentTool = 'select';
			ct.showSegmentation('update');
		case 'click'
			%Called when mouse pressed
      [x,y] = mygetcurrentpoint(gca);

			type = get(gcf,'SelectionType');
			switch helpVar.CurrentTool
				case 'select'
					%--- Selection tool
					switch type
						case 'normal'
							%Mark pen/window
							helpVar.pencolor = arg1;
							selectmotion; %Call once to update

							%Start selection tool if kept down
							set(gcf,'WindowButtonMotionFcn','ct.showSegmentation(''selectmotion'')');
							set(gcf,'WindowButtonUpFcn','ct.showSegmentation(''selectbuttonup'')');
						case 'extend'
							%Start paning tool
							helpVar.pencolor = arg1;
							set(gcf,'WindowButtonMotionFcn','ct.showSegmentation(''panmotion'')');
							set(gcf,'WindowButtonUpFcn','ct.showSegmentation(''panbuttonup'')');
					end;
					
			end; %switch clause
		case 'refresh'
			set(handles.refreshicon,'state','off');
			helpVar.RZoomState = [0.5 helpVar.XSize-0.5 0.5 helpVar.YSize-0.5];
			helpVar.GZoomState = [0.5 helpVar.YSize-0.5 0.5 helpVar.ZSize-0.5];
			helpVar.BZoomState = [0.5 helpVar.XSize-0.5 0.5 helpVar.ZSize-0.5];
			helpVar.LevelSet.Zoom = 1;
			ct.showSegmentation('update');
	end
end

%---------------------------------------
function showSegmentationkeypressed(fignum,evnt)
%---------------------------------------
global setCell helpVar

key = getkey(evnt);

switch key
  case 'z'
    %Zoom in (z)
    ct.showSegmentation('zoomin');
  case 'x'
    %Zoom out (x)
		ct.showSegmentation('zoomout');
	case 'uparrow'
		switch helpVar.pencolor
			case 'r'
				ct.showSegmentation('slice','r',1)
			case 'g'
				ct.showSegmentation('slice','g',1)
			case 'b'
				ct.showSegmentation('slice','b',1)
		end
	case 'downarrow'
		switch helpVar.pencolor
			case 'r'
				ct.showSegmentation('slice','r',-1)
			case 'g'
				ct.showSegmentation('slice','g',-1)
			case 'b'
				ct.showSegmentation('slice','b',-1)
		end
end;

%------------------------------------------------------------------
function im = remapandoverlaySegmentation(temp,overlaySegmentation)
%------------------------------------------------------------------
global setCell helpVar

%Remap to correct colorscale and add overlay if user wants that.

%Remap data
if isequal(helpVar.ViewIm,'magnitude')
  tempr = calcfunctions('remapuint8',temp);
  tempg = tempr;
  tempb = tempr;
else
  temp = temp/4+0.5;
  tempr = segment('remap',temp,helpVar.Colormap(:,1));
  tempg = segment('remap',temp,helpVar.Colormap(:,2));
  tempb = segment('remap',temp,helpVar.Colormap(:,3));
end;
%Add overlay
if helpVar.ViewSelection
  tempr(overlaySegmentation)=1;
  %tempg(overlaybw)=0;
  %tempb(overlaybw)=0; %red
end;
im = cat(3,tempr,tempg,tempb);

%-------------------------------
function showSegmentationDefault
%-------------------------------
global helpVar setCell

helpVar.pencolor = ''; %contains r,g,b depending which window drawing in
helpVar.RZoomState = [0.5 helpVar.XSize-0.5 0.5 helpVar.YSize-0.5];
helpVar.GZoomState = [0.5 helpVar.YSize-0.5 0.5 helpVar.ZSize-0.5];
helpVar.BZoomState = [0.5 helpVar.XSize-0.5 0.5 helpVar.ZSize-0.5];
helpVar.Zoom = 1;
helpVar.CurrentTool = 'select';
helpVar.ViewIm = 'magnitude'; %can also be speed
helpVar.ViewMIP = false; %view mip image
helpVar.ViewInteraction = true; %View manual interaction
helpVar.ViewSelection = false; %View selection
helpVar.ViewOutline = true; %View outline

%-----------------------------------------------------
function [zoomstate] = zoomhelper(zoomstate,f,sz)
%-----------------------------------------------------
%imagesize
xsize=sz(1);
ysize=sz(2);
%Get old position
temp = zoomstate;
oldxspan=temp(2)-temp(1);
oldyspan=temp(4)-temp(3);

xlim = [...
  0.5*(temp(1)+temp(2))-0.5*(temp(2)-temp(1))/f ...
  0.5*(temp(1)+temp(2))+0.5*(temp(2)-temp(1))/f];
%f = f*0.5;
ylim = [...
  0.5*(temp(3)+temp(4))-0.5*(temp(4)-temp(3))/f ...
  0.5*(temp(3)+temp(4))+0.5*(temp(4)-temp(3))/f];

xspan=(xlim(2)-xlim(1));
yspan=(ylim(2)-ylim(1));
if f>1
  if xsize>ysize
    if xspan>oldyspan
      ylim=[0.5 ysize-0.5];
    elseif oldyspan>xspan
      ylim=[0.5*(temp(3)+temp(4)) - 0.5*(xlim(2)-xlim(1))...
        0.5*(temp(3)+temp(4)) + 0.5*(xlim(2)-xlim(1))];
    end
  else
    if yspan>oldxspan
      xlim=[0.5 xsize-0.5];
    elseif oldxspan>yspan
      xlim=[0.5*(temp(1)+temp(2)) - 0.5*(ylim(2)-ylim(1))...
        0.5*(temp(1)+temp(2)) + 0.5*(ylim(2)-ylim(1))];
    end
  end
else
  if xsize>ysize
    if xspan<=xsize
      if oldxspan>yspan || yspan>=ysize
        ylim=[0.5 ysize-0.5];
      elseif yspan>oldxspan
        ylim=[0.5*(temp(3)+temp(4)) - 0.5*(xlim(2)-xlim(1))...
          0.5*(temp(3)+temp(4)) + 0.5*(xlim(2)-xlim(1))];
      end
    end
  else
    if yspan<=ysize
      if oldyspan>xspan ||xspan>=xsize
        xlim=[0.5 xsize-0.5];
      elseif xspan>oldyspan
        xlim=[0.5*(temp(1)+temp(2)) - 0.5*(ylim(2)-ylim(1))...
          0.5*(temp(1)+temp(2)) + 0.5*(ylim(2)-ylim(1))];
      end
    end
  end
end
         
zoomstate = [xlim(:);ylim(:)];
%------------------------------
function selectbuttonup
%------------------------------

set(gcf,'WindowButtonMotionFcn','');
set(gcf,'WindowButtonUpFcn','');
ct.showSegmentation('update');

%----------------------------
function selectmotion
%----------------------------
global helpVar
%rewritten to be able to use in timeresolve volumes (JS)

[x,y] = mygetcurrentpoint(gca);


switch helpVar.pencolor
	case 'r'
		helpVar.gslice = round(x);
		helpVar.bslice = round(y);
	case 'g'
		helpVar.bslice = round(x);
		helpVar.rslice = round(y);
	case 'b'
		helpVar.gslice = round(x);
		helpVar.rslice = round(y);
end;

%Move the line
ct.showSegmentation('update')


%---------------------------
function panbuttonup
%---------------------------
global helpVar

test='i panbuttonup'

set(gcf,'WindowButtonMotionFcn','');
set(gcf,'WindowButtonUpFcn','');

ax = gca;
xlim = get(ax,'xlim');
ylim = get(ax,'ylim');

switch helpVar.pencolor
  case 'r'
    helpVar.RZoomState = [xlim ylim];
  case 'g'
    helpVar.GZoomState = [xlim ylim];
  case 'b'
    helpVar.BZoomState = [xlim ylim];
end;
ct.showSegmentation('update');
panmotion('reset'); %Reset for next time

%--------------------------------
function panmotion(reset)
%--------------------------------
persistent startpos

test='i panmotion'

if nargin==1
  startpos = [];
  return;
end;

%Get mouse position
[p(1),p(2)] = mygetcurrentpoint(gca);

if isempty(startpos)
  startpos = p;
end;

xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
set(gca,...
  'xlim',xlim+startpos(1)-p(1),...
  'ylim',ylim+startpos(2)-p(2));


