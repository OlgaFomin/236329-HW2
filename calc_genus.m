function [genus, boundary_edges] = calc_genus(adjacency_matrix, triangles, vertices)

    
    [v1_vec,v2_vec] = find(adjacency_matrix==1);
    
    faces_near_edges_mat = sparse(length(v1_vec));
    
    for i=1:length(v1_vec)
       v1 = v1_vec(i);
       v2 = v2_vec(i);
       
       if v1 == v2
           continue
       end
       
       [f1,~] = find(triangles==v1);
       [f2,~] = find(triangles==v2);
       
       faces_near_edges_mat(v1,v2) = length(intersect(f1,f2));
        
    end
    
    boundary_edges = length(find(faces_near_edges_mat==1));
    
    V = size(vertices,1);
    F = size(triangles,1);
    E = (length(v1_vec) - size(adjacency_matrix,1))/2;
    
    chi = V - E + F;
    genus = 0.5*(2 - boundary_edges - chi);
    
end