%Exports a mesh to Ensight Gold format
function matlab2ensightgolddeformingmesh(pathname,filename,filesuffix,x,y,z,timesteps)
	global NO

% 	[rotatedx, rotatedy] = rotateendorings(x,y);
% 	x(:,:,:) = rotatedx;
% 	y(:,:,:) = rotatedy;

    %Fix timeframes without segmentation issues: 
	%----------------------------------------------------------------------
	isnanendo=isnan(sum(x(:,:,1),1));
	tfvector=(1:size(x,2));
	
	isnotnanendo_indice = tfvector(~isnanendo);
	isnanendo_indice = tfvector(isnanendo)';
	
	if isempty(isnotnanendo_indice)
		mymsgbox('Error: No segmentation detected. Export aborted.');
		return;
	end
	
	if ~isempty(isnanendo_indice)
		%Copy closest segmentation to timeframes missing segmentation:
		[~, closestendo] = min(abs(bsxfun(@minus,isnanendo_indice,isnotnanendo_indice)),[],2);
		x(:,isnanendo_indice,:) = x(:,isnotnanendo_indice(closestendo),:);
		y(:,isnanendo_indice,:) = y(:,isnotnanendo_indice(closestendo),:);
	end
	
	%if singleslice: duplicate it. 
	if ((numel(z)==1) && (size(x,3)==1) && (size(y,3)==1))
		z=[z,1];
		x = x(:,:,ones(2,1));
		y = y(:,:,ones(2,1));
	end
	%----------------------------------------------------------------------
	
	numsteps = length(timesteps); %Number of time steps

	% ---- Create Case file

	%Create and open the file
	disp('Creating case file.');
	casefilefid = fopen([pathname filesep filename filesuffix '.case'],'w');
	if isequal(casefilefid,-1)
		error(sprintf('Could not open %s for output.',[pathname filesep filename '.case']));
	end;

	%Write FORMAT
	fprintf(casefilefid,'FORMAT\n\n');
	fprintf(casefilefid,'type:\tensight gold\n\n');

	%Write GEOMETRY section
	if numsteps>1
		wildcard = '****';
	else
		wildcard = '';
	end;
	fprintf(casefilefid,'GEOMETRY\n\n');
	fprintf(casefilefid,'model:\t\t1\t%s%s.geo%s\n\n', filename, filesuffix, wildcard);

	%Write TIME section
	fprintf(casefilefid,'TIME\n');
	fprintf(casefilefid,'time set:\t\t1\tModel\n'); %default 1
	fprintf(casefilefid,'number of steps:\t%d\n',numsteps);
	fprintf(casefilefid,'filename start number:\t1\n'); %default 1
	fprintf(casefilefid,'filename increment:\t1\n'); %default 1
	fprintf(casefilefid,'time values:\t%12.5e ',timesteps(1));

	if numsteps>1
		for loop=2:numsteps
			fprintf(casefilefid,'%12.5e ',timesteps(loop)); %default 1
			if mod(loop,4)==0
				fprintf(casefilefid,'\n'); 
			end;
		end;
		fprintf(casefilefid,'\n');
	end;

	%Close file
	fclose(casefilefid);

	%%% ---- Create geometry file
	for timeframe = 1:numsteps
		newname = sprintf('%s%s%s%s.geo%04d',pathname,filesep,filename,filesuffix,timeframe);
		geofilefid = fopen(newname,'w'); %  geofilefid = fopen(newname,'w');
		if isequal(geofilefid,-1)
			error(sprintf('Could not open %s for output.',newname));
		end;

        % Print Header
        fprintf(geofilefid, '%s\n', '3D and 2D parts');
        fprintf(geofilefid, '%s\n', sprintf('Exported by %s', mfilename));
        fprintf(geofilefid, '%s\n', 'node id given');
        fprintf(geofilefid, '%s\n', 'element id given');

        %Start with the part
        fprintf(geofilefid, '%s\n', 'part');

        %Write part number
        fprintf(geofilefid, '         1\n');

        %Write description
        fprintf(geofilefid,'description for part1\n');

        fprintf(geofilefid,'coordinates\n');

        %Write number of nodes
        nodes = 0;
        for nodecoordinates = 1:size(x, 3)
            if not(isnan(x(1, timeframe, nodecoordinates)))
                nodes = nodes+size(x, 1);
            end;
        end;
        fprintf(geofilefid,'%10d\n', nodes);

        % write node ids and set up lookup structure
        nodeindexlookup = java.util.Hashtable;
        nodeid = 1;
        for currentcolumn = 1:size(x, 3)
            nodeindexlookup.put(currentcolumn, java.util.Hashtable);
            for currentrow = 1:size(x, 1)
                if not(isnan(x(currentrow, timeframe, currentcolumn)))
                    nodeindexlookup.get(currentcolumn).put(currentrow, nodeid);
                    fprintf(geofilefid,'%10d\n', nodeid);
                    nodeid=nodeid+1;
                end;
            end;
        end;

        % write x node coordinates
        for currentcolumn = 1:size(x, 3)
            for currentrow = 1:size(x, 1)
                if not(isnan(x(currentrow, timeframe, currentcolumn)))
                    coords = calcfunctions('xyz2rlapfh', NO, x(currentrow, timeframe, currentcolumn), y(currentrow, timeframe, currentcolumn), z(currentcolumn));
                    expstring = getcrossplatformexponential(coords(1)/1000);
                    fprintf(geofilefid,'%s\n', expstring);
                end;
            end;
        end;

        % write y node coordinates
        for currentcolumn = 1:size(x, 3)
            for currentrow = 1:size(x, 1)
                if not(isnan(x(currentrow, timeframe, currentcolumn)))
                    coords = calcfunctions('xyz2rlapfh', NO,x(currentrow, timeframe, currentcolumn), y(currentrow, timeframe, currentcolumn), z(currentcolumn));
                    expstring = getcrossplatformexponential(coords(2)/1000);
                    fprintf(geofilefid,'%s\n', expstring);
                end;
            end;
        end;

        % write z node coordinates
        for currentcolumn = 1:size(x, 3)
            for currentrow = 1:size(x, 1)
                if not(isnan(x(currentrow, timeframe, currentcolumn)))
                    coords = calcfunctions('xyz2rlapfh', NO,x(currentrow, timeframe, currentcolumn), y(currentrow, timeframe, currentcolumn), z(currentcolumn));
                    expstring = getcrossplatformexponential(coords(3)/1000);
                    fprintf(geofilefid,'%s\n', expstring);
                end;
            end;
        end;

        % write number of quads
        % |NaN|*|*|NaN|
        % |NaN|*|*|NaN|       <------- example with 5 quads
        % |NaN|*|*|NaN|        every node that has a node on the 
        % |NaN|*|*|NaN|        next column and one on the next row 
        % |NaN|*|*|NaN|        is a quad
        % |NaN|*|*|NaN|
        quads = 0;
        for currentcolumn = 1:size(x, 3)-1 % go through all columns untill the second to last one
            for currentrow = 1:size(x, 1)-1
                nextrow = not(isnan(x(currentrow+1, timeframe, currentcolumn)));
                nextcolumn = not(isnan(x(currentrow, timeframe, currentcolumn+1)));
                if nextrow && nextcolumn
                    quads=quads+1;
                end;
            end;
            % count the quad between the beginning and end as well
            if nextrow && nextcolumn 
                quads=quads+1;
            end;
        end;
        fprintf(geofilefid,'quad4\n');
        fprintf(geofilefid,'%10d\n', quads);
        for quadid=1:quads
            fprintf(geofilefid,'%10d\n', quadid);
        end;

        % write the quad topology 
        for currentcolumn = 1:size(x, 3)-1 % go through all columns untill the second to last one
            for currentrow = 1:size(x, 1)-1
                nextrow = not(isnan(x(currentrow+1, timeframe, currentcolumn)));
                nextcolumn = not(isnan(x(currentrow, timeframe, currentcolumn+1)));
                if nextrow && nextcolumn
                  fprintf(geofilefid,'%10d', nodeindexlookup.get(currentcolumn).get(currentrow));
                  fprintf(geofilefid,'%10d', nodeindexlookup.get(currentcolumn+1).get(currentrow));
                  fprintf(geofilefid,'%10d', nodeindexlookup.get(currentcolumn+1).get(currentrow+1));
                  fprintf(geofilefid,'%10d\n', nodeindexlookup.get(currentcolumn).get(currentrow+1));
                end;
            end;
            % add the quad between the last and first row as well
            if nextrow && nextcolumn
              fprintf(geofilefid,'%10d', nodeindexlookup.get(currentcolumn).get(size(x, 1)));
              fprintf(geofilefid,'%10d', nodeindexlookup.get(currentcolumn+1).get(size(x, 1)));
              fprintf(geofilefid,'%10d', nodeindexlookup.get(currentcolumn+1).get(1));
              fprintf(geofilefid,'%10d\n', nodeindexlookup.get(currentcolumn).get(1));
            end;
        end;

        % Close geo file
        fclose(geofilefid);
    end;

function returnedExponentialString = getcrossplatformexponential(number)
    if(number >= 0)
        returnedExponentialString = sprintf('%12.5e', number);
    else
        returnedExponentialString = sprintf('%12.4e', number);
    end;
    %if ispc
    %    returnedExponentialString = strrep(returnedExponentialString, 'e+0','e+');
    %    returnedExponentialString = strrep(returnedExponentialString, 'e-0','e-');
    %end;
    
function [rotatedx, rotatedy] = rotatering(referencex, referencey, xring, yring) 
	rotatedx = zeros(size(xring));
    rotatedy = zeros(size(yring));
	smallestdistance = sqrt((referencex-xring(1))^2+(referencey-yring(1))^2);
	smallestposition = 1;
	for distanceloopcounter = 2:size(xring)
		newdistance = sqrt((referencex-xring(distanceloopcounter))^2+(referencey-yring(distanceloopcounter))^2);
		if(newdistance < smallestdistance)
			smallestdistance = newdistance;
			smallestposition = distanceloopcounter;
		end;
	end;
	
	writtenelements = 1;
	for i = smallestposition:size(xring)
		rotatedx(writtenelements) = xring(i);
		rotatedy(writtenelements) = yring(i);
		writtenelements = writtenelements+1;
	end;

    if smallestposition ~= 1
        for i = 1:(smallestposition-1)
            rotatedx(writtenelements) = xring(i);
            rotatedy(writtenelements) = yring(i);
            writtenelements = writtenelements+1;
        end;
    end;

function [outdatax, outdatay] = rotateendorings(indatax, indatay)
	outdatax = zeros(size(indatax));
    outdatax(:,:,:) = NaN;
    outdatay = zeros(size(indatay));
    outdatay(:,:,:) = NaN;
	
	% Go through the entire input dataset and rotate the endo rings so
	% they are all oriented the same
	for currenttime = 1:size(indatax, 2)
		% Find the first column that is not full of NaNs.
		% The first element of this column will be our reference position
		for referencecolumn = 1:size(indatax, 3)
			if not(isnan(indatax(1, currenttime, referencecolumn)))
				referencex = indatax(1, currenttime, referencecolumn);
				referencey = indatay(1, currenttime, referencecolumn);
				break;
			end;
		end;
		
		% Rotate the next column ring after first column and then the next after that
		for currentcolumn = (referencecolumn+1):size(indatax, 3)
			if not(isnan(indatax(1, currenttime, currentcolumn)))
				[rotatedringx, rotatedringy] = rotatering(referencex, referencey, indatax(:,currenttime, currentcolumn), indatay(:,currenttime, currentcolumn));
				outdatax(:,currenttime, currentcolumn) = rotatedringx;
				outdatay(:,currenttime, currentcolumn) = rotatedringy;
				referencex = outdatax(1, currenttime, currentcolumn);
				referencey = outdatay(1, currenttime, currentcolumn);
			end;
		end;
	end;
	
	function [outdatax, outdatay] = rotateendorings2(indatax, indatay)
	outdatax = nan(size(indatax));
    outdatay = nan(size(indatay));
	
	% Go through the entire input dataset and rotate the endo rings so
	% they are all oriented the same
	for currenttime = 1:size(indatax, 2)
		% Find the first column that is not full of NaNs.
		% The first element of this column will be our reference position
		for referencecolumn = 1:size(indatax, 3)
			if not(isnan(indatax(1, currenttime, referencecolumn)))
				referencex = indatax(1, currenttime, referencecolumn);
				referencey = indatay(1, currenttime, referencecolumn);
				break;
			end;
		end;
		
		% Rotate the next column ring after first column and then the next after that
		for currentcolumn = (referencecolumn+1):size(indatax, 3)
			if not(isnan(indatax(1, currenttime, currentcolumn)))
				[rotatedringx, rotatedringy] = rotatering(referencex, referencey, indatax(:,currenttime, currentcolumn), indatay(:,currenttime, currentcolumn));
				outdatax(:,currenttime, currentcolumn) = rotatedringx;
				outdatay(:,currenttime, currentcolumn) = rotatedringy;
				referencex = outdatax(1, currenttime, currentcolumn);
				referencey = outdatay(1, currenttime, currentcolumn);
			end;
		end;
	end;
	
