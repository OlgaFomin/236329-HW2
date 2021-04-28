function visualize_fun(vertices, triangles, fun, type)

    if strcmp(type,"face")
        figure()
        patch('Faces',triangles,'Vertices',vertices,'FaceVertexCData',fun,'FaceColor','flat');
        colorbar
    elseif strcmp(type,"vertex")
        figure()
        patch('Faces',triangles,'Vertices',vertices,'FaceVertexCData',fun,'EdgeColor','interp','FaceColor','none');
        colorbar
    else
        disp 'unvalid type. type should be face/vertex'
    end
    
end