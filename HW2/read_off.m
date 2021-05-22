function [vertices, faces] = read_off(off_filepath)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parses .off file into 2 matrices - vertices, faces.                %
    % Assuming that the .off file:                                       %
    %   begins with 'OFF'                                                %
    %   doesn't contain comments                                         %
    %   doesn't contain colors                                           %
    %                                                                    %
    % off_filepath = filepath to .off file                               %
    % vertices = [Vx3] matrix containing vertices coordinates            %
    % faces = [Fx3] matrix containing the vertices of each face          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('>> Reading .off file - %s\n', off_filepath)
    fid = fopen(off_filepath,'r');
    str = fgets(fid);   % -1 if eof
    if ~strcmp(str(1:3), 'OFF')
        error('File assumed to begin with OFF, hence cannot be read.');    
    end
    
    str = fgets(fid);
    
    [tok,str] = strtok(str); nv = str2num(tok);
    [tok,str] = strtok(str); nf = str2num(tok);
    
    [vertices,cnt] = fscanf(fid,'%f %f %f', 3*nv);
    if cnt~=3*nv
        warning('Problem in reading vertices.');
    end
    vertices = reshape(vertices, 3, cnt/3);
    vertices = vertices';
    
    [faces,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nf);
    if cnt~=4*nf
        warning('Problem in reading faces.');
    end
    faces = reshape(faces, 4, cnt/4);
    faces = faces(2:4,:)+1;
    faces = faces';
    
    fclose(fid);

end