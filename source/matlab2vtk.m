function matlab2vtk(pathname, filename, datastruct)
    function [delta, originalDelta, originPre, originPost, rotation] = getInformation()
        originalDelta = datastruct(1).Delta;
        delta = abs(originalDelta);
        originPre = double(datastruct(1).Origin);
        originPost = double(datastruct(1).Origin);
        dirdim = double(datastruct(1).DirDim);
        vx = dirdim(2,:);
        vy = dirdim(1,:);
        vz = dirdim(3,:);
        [rotation(1), rotation(2), rotation(3)] = getangles((vx/norm(vx))', (vy/norm(vy))', (vz/norm(vz))');
        rotation = rotation*(180/pi);
        if datastruct(1).Delta(1) < 0
            originPost(1) = originPost(1)-delta(1)*(datastruct(1).Size(1)-1);
            rotation(1) = rotation(1);
            disp(1);
        end
        if datastruct(1).Delta(2) < 0
            originPost(2) = originPost(2)-delta(2)*(datastruct(1).Size(2)-1);
            rotation(2) = rotation(2);
            disp(2);
        end
        if datastruct(1).Delta(3) < 0
            originPost(3) = originPost(3)-delta(3)*(datastruct(1).Size(3)-1);
            rotation(3) = rotation(3);
            disp(3);
        end
        rotTmp = rotation(1);
        rotation(1) = rotation(2);
        rotation(2) = rotTmp;
        
        disp(originPre);
        disp(originPost);
        disp(rotation);
        disp(delta);
        disp('------');
    end
    function [L, P, H, delta1, delta2, delta3, vx, vy, vz] = get_lph(partno)
        dirdim = double(datastruct(partno).DirDim);
        dd1 = dirdim(:,1);
        dd2 = dirdim(:,2);
        dd3 = dirdim(:,3);

        voxspace = double(datastruct(partno).Delta); % voxel spacing
        origin = double(datastruct(partno).Origin);

        %  Transfer coordinates to patient reference system
        varsize = datastruct(partno).Size;
        if length(varsize)>2 %3D
            [i j k] = ndgrid(0:(varsize(1)-1), 0:(varsize(2)-1), 0:(varsize(3)-1));
        else %2D
            [i j k] = ndgrid(0:(varsize(1)-1), 0:(varsize(2)-1), -.5:.5);
            voxspace(3) = 0.001;
            varsize(3) = 1;
        end

        L = i*voxspace(1)*dd1(1) + j*voxspace(2)*dd2(1) + k*voxspace(3)*dd3(1) + origin(1);
        P = i*voxspace(1)*dd1(2) + j*voxspace(2)*dd2(2) + k*voxspace(3)*dd3(2) + origin(2);
        H = i*voxspace(1)*dd1(3) + j*voxspace(2)*dd2(3) + k*voxspace(3)*dd3(3) + origin(3);
        
        
        vx = [L(2,1,1)-L(1,1,1), P(2,1,1)-P(1,1,1), H(2,1,1)-H(1,1,1)];
        vy = [L(1,2,1)-L(1,1,1), P(1,2,1)-P(1,1,1), H(1,2,1)-H(1,1,1)];
        vz = [L(1,1,2)-L(1,1,1), P(1,1,2)-P(1,1,1), H(1,1,2)-H(1,1,1)];
        vx = vx/norm(vx);
        vy = vy/norm(vy);
        vz = vz/norm(vz);
        
        if voxspace(1) < 0.0
            vx = vx*-1;
        end
        if voxspace(2) < 0.0
            vy = vy*-1;
        end
        if voxspace(3) < 0.0
            vz = vz*-1;
        end
        clear i j k
        L = L(:);
        P = P(:);
        H = H(:);
        delta1 = datastruct(partno).Delta(1);
        delta2 = datastruct(partno).Delta(2);
        delta3 = datastruct(partno).Delta(3);
    end

    function [x1, x2, y1, y2, z1, z2, delta1, delta2, delta3, L, P, H, vx, vy, vz] = get_whole_extents_of_part(partno)
        [L, P, H, delta1, delta2, delta3, vx, vy, vz] = get_lph(partno);
        x1 = 0;
        x2 = datastruct(partno).Size(1)-1;
        y1 = 0;
        y2 = datastruct(partno).Size(2)-1;
        z1 = 0;
        z2 = datastruct(partno).Size(3)-1;
    end
    
    function export_vector_array(partno, geofilefid, current_time)
        fprintf(geofilefid, '\t\t\t\t<DataArray type="Float32" Name="VectorArray" format="ascii">');
        fprintf(geofilefid, '% f ', datastruct(partno).Data(:, :, :, current_time, 1));
        fprintf(geofilefid, '% f ', datastruct(partno).Data(:, :, :, current_time, 2));
        fprintf(geofilefid, '% f ', datastruct(partno).Data(:, :, :, current_time, 3));
        fprintf(geofilefid, '</DataArray>\n');
    end
    
    function export_2d_magnitude_array(partno, geofilefid, current_time)
        fprintf(geofilefid, '\t\t\t\t<DataArray type="Float32" Name="MagnitudeArray" format="ascii">');
        fprintf(geofilefid, '% f ', datastruct(partno).Data(:, :, :, current_time, 1));
        fprintf(geofilefid, '</DataArray>\n');
    end

    function export_points(L, P, H)
        fprintf(geofilefid, '\t\t\t<Points>\n');
            fprintf(geofilefid, '\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">');
            %fprintf(geofilefid, '% f ', reshape([L;P;H],1,[]));
            fprintf(geofilefid, '% f ', L);
            fprintf(geofilefid, '% f ', P);
            fprintf(geofilefid, '% f ', H);
            fprintf(geofilefid, '</DataArray>\n');
        fprintf(geofilefid, '\t\t\t</Points>\n');
    end

    function export_magnitude_array(partno, geofilefid, current_time, originalDelta)
        fprintf(geofilefid, '\t\t\t\t<DataArray type="Float32" Name="MagnitudeArray" format="ascii">');
        x = datastruct(partno).Data(:, :, :, current_time, 2);
        if originalDelta(2) < 0
            x = flip(z, 2);
        end
        y = datastruct(partno).Data(:, :, :, current_time, 1);
        if originalDelta(1) < 0
            y = flip(y, 1);
        end
        z = datastruct(partno).Data(:, :, :, current_time, 3);
        if originalDelta(3) < 0
            z = flip(z, 3);
        end
        
        length = sqrt(x.*x+y.*y+z.*z);
        length = permute(length,[2 1 3]);
        fprintf(geofilefid, '% f ', length);
        fprintf(geofilefid, '</DataArray>\n');
    end

    function export_maximum_magnitude_array(partno, geofilefid, current_time)
        fprintf(geofilefid, '\t\t\t\t<DataArray type="Float32" Name="MaximumMagnitudeArray" format="ascii">');
        outValues = zeros(1, size(datastruct(partno).Data(), 3)*size(datastruct(partno).Data(), 2)*size(datastruct(partno).Data(), 1)*3);
        outValuesIndex = 1;
        for k=fliplr(1:size(datastruct(partno).Data(), 3))
            for j=1:size(datastruct(partno).Data(), 1) %Switched 1 & 2, Einar 2016-08-02
                for i=1:size(datastruct(partno).Data(), 2)
                    looptime = 1:datastruct.NumberOfTimeSteps;
                    x = max(squeeze(datastruct(partno).Data(j, i, k, looptime, 1)));
                    y = max(squeeze(datastruct(partno).Data(j, i, k, looptime, 2)));
                    z = max(squeeze(datastruct(partno).Data(j, i, k, looptime, 3)));
                    maximumvalue = sqrt(x*x+y*y+z*z);
                    outValues(outValuesIndex) = maximumvalue;
                    outValuesIndex = outValuesIndex+1;
                    fprintf(geofilefid, '%f ', maximumvalue);
                end
            end
        end
        x = datastruct(partno).Data(:, :, :, current_time, 1);
        y = datastruct(partno).Data(:, :, :, current_time, 2);
        z = datastruct(partno).Data(:, :, :, current_time, 3);
        length = sqrt(x.*x+y.*y+z.*z);
        fprintf(geofilefid, '% f ', outValues);
        clear outValues;
        fprintf(geofilefid, '</DataArray>\n');                
    end
    
    function export_parts(partno, geofilefid, current_time, originalDelta)
        fprintf(geofilefid, '\t\t\t<PointData>\n');
            %if length(datastruct(partno).Data(1,1,1,1,:)) == 1
            %    export_2d_magnitude_array(partno, geofilefid, current_time);
            %end;
            if length(datastruct(partno).Data(1,1,1,1,:)) == 3
                % Make vector data to scalar, used for volume rendering
                %export_vector_array(partno, geofilefid, current_time);
                %export_magnitude_array(partno, geofilefid, current_time, originalDelta);
                export_maximum_magnitude_array(partno, geofilefid, current_time);
            end
        fprintf(geofilefid, '\t\t\t</PointData>\n');
        %export_points(L, P, H);
    end

    numparts = length(datastruct);
    %Create and open the file
    disp('Creating vtk file.');
    casefilefid = fopen([pathname filesep filename '.pvd'],'w');
    if isequal(casefilefid,-1)
      error(sprintf('Could not open %s for output.',[pathname filesep filename ...
        '.case']));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [delta, originalDelta, originPositionPre, originPositionpost, rotation] = getInformation();
    %Write FORMAT
    fprintf(casefilefid,'<?xml version="1.0"?>\n');
    fprintf(casefilefid,'<!--%f %f %f-->\n', originPositionPre(1), originPositionPre(2), originPositionPre(3));
    fprintf(casefilefid,'<!--%f %f %f-->\n', rotation(1), rotation(2), rotation(3));
    fprintf(casefilefid,'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n');

    fprintf(casefilefid,'\t<Collection>\n');
    %for loop=1:datastruct.NumberOfTimeSteps
        %fprintf(casefilefid,'\t\t<DataSet timestep="%12.5e" group="" part="0" file="%s%d.vti"/>\n', datastruct.TimeSteps(loop), filename, loop); %default 1
        fprintf(casefilefid,'\t\t<DataSet timestep="%12.5e" group="" part="0" file="%s%d.vti"/>\n', datastruct.TimeSteps(1), filename, 1); %default 1
    %end;
    fprintf(casefilefid,'\t</Collection>\n');

    fprintf(casefilefid,'</VTKFile>');

    %Close file
    fclose(casefilefid);

    %for current_time=1:datastruct.NumberOfTimeSteps
        current_time=1;
        newname = sprintf('%s%s%s%d.vti',pathname,filesep,filename,current_time);
        geofilefid = fopen(newname,'w'); %  geofilefid = fopen(newname,'wb');
        if isequal(geofilefid,-1)
          error(sprintf('Could not open %s for output.',newname));
        end

        fprintf(geofilefid,'<?xml version="1.0"?>\n');
        fprintf(geofilefid,'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n');
        [y1, y2, x1, x2, z1, z2] = get_whole_extents_of_part(1);
        fprintf(geofilefid, '\t<ImageData WholeExtent="%d %d %d %d %d %d" Origin="%f %f %f" Spacing="%f %f %f">\n', x1, x2, y1, y2, z1, z2, originPositionpost(1), originPositionpost(2), originPositionpost(3), delta(1), delta(2), delta(3));
        for partno = 1:numparts
            fprintf(geofilefid, '\t\t<Piece Extent="%d %d %d %d %d %d">\n', x1, x2, y1, y2, z1, z2);
            export_parts(partno, geofilefid, current_time, originalDelta);
            fprintf(geofilefid, '\t\t</Piece>\n');
        end
        fprintf(geofilefid, '\t</ImageData>\n');
        fprintf(geofilefid,'</VTKFile>\n');
        fclose(geofilefid);
    %end;
    disp('Finished exporting to vtk');
end
