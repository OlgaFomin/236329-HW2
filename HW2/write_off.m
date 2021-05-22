function write_off(vertices, faces, out_filepath)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Writes .off file out of 2 matrices - vertices, faces.              %
    %                                                                    %
    % out_filepath = filepath output to .off file                        %
    % vertices = [Vx3] matrix containing vertices coordinates            %
    % faces = [Fx3] matrix containing the vertices of each face          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('>> Writing .off file - %s\n', out_filepath)
    nv = size(vertices,1);
    nf = size(faces,1);
    
    fid = fopen(out_filepath,'w');
    fprintf(fid,'%s\n','OFF'); 
    fprintf(fid,'%u %u 0\n',[nv nf]);
    fprintf(fid,'%f %f %f\n',vertices');
    fprintf(fid,'3 %u %u %u\n',faces');
    fclose(fid);
    
end