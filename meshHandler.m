classdef meshHandler < handle
    
    properties
        
        adjacency_matrix
        shared_matrix
        
    end
    
    methods
        
        function obj = meshHandler(vertices_array, triangles_array)
            
            nv = size(vertices_array,1);
            nf = size(triangles_array,1);
            
            % Create shared matrix:
            shared_matrix.vertices_array = vertices_array;
            shared_matrix.triangles_array = triangles_array;
            obj.shared_matrix = shared_matrix;
            
            % Create adjacency matrix:
            adjacency_matrix = sparse(nv,nv);
            for i = 1:nf
                adjacency_matrix(triangles_array(i,1), triangles_array(i,2)) = 1;
                adjacency_matrix(triangles_array(i,2), triangles_array(i,3)) = 1;
                adjacency_matrix(triangles_array(i,3), triangles_array(i,1)) = 1;
            end
            obj.adjacency_matrix = adjacency_matrix + speye(size(adjacency_matrix));
            
        end
        
    end
    
end
