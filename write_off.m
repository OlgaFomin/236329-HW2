function write_off(vertices, triangles, out)

    num_vertices = size(vertices,1);
    num_triangles = size(triangles,1);
    
    fid = fopen(out,'w');
    fprintf(fid,'%s\n','OFF'); 
    fprintf(fid,'%u %u 0\n',[num_vertices num_triangles]);
    fprintf(fid,'%u %u %u\n',vertices);
    fprintf(fid,'3 %u %u %u\n',triangles);
    fclose(fid);
    
end
