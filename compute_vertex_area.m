function vertex_area = compute_vertex_area(vertex_array, traingles_array)

    vertex_area = [];
    for i = 1:length(vertex_array)
        [face_indices, dim] = find(traingles_array==i);
        adjacent_triangles = traingles_array(face_indices, :);
        area = compute_triangle_area(vertex_array, adjacent_triangles);
        vertex_area = [vertex_area; sum(area)/3];       
    
    end

end