function faces_near_edges_mat = faces_near_edges(adjacency_matrix, triangles)

    
    [v1_vec,v2_vec] = find(adjacency_matrix==1);
    
    faces_near_edges_mat = sparse(length(v1_vec));
    
    for i=1:length(v1_vec)
       v1 = v1_vec(i);
       v2 = v2_vec(i);
       
       [f1,~] = find(triangles==v1);
       [f2,~] = find(triangles==v2);
       
       faces_near_edges_mat(v1,v2) = length(intersect(f1,f2));
        
    end
    
end