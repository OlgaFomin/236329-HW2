classdef mesh < handle
    
    properties 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mesh Class                                                    %
        %                                                               %
        % adjacency_matrix - [VxV] sparse matrix. Contains '1' for      %
        % adjacent vertices.                                            %
        % shared_matrix.vertices = [Vx3] matrix containing vertices     % 
        % coordinates                                                   %
        % shared_matrix.faces = [Fx3] matrix containing the vertices    %
        % of each face                                                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vertices 
        faces
        v_adj
        vf_adj
    end
    
    methods  
        function obj = mesh(vertices, faces)
            
            nv = size(vertices,1);
            nf = size(faces,1);
            
            % Create shared matrix:
            obj.faces     = faces;
            obj.vertices  = vertices;
            
            % Create vertex-face adjacency matrix [VxF]:
            vf_adj      = sparse(faces, repmat((1:nf)', 1, 3), ones(nf,3), nv, nf);
            obj.vf_adj  = vf_adj;
            
            % Create adjacency matrix [VxV]:
            v_adj       = vf_adj*vf_adj';
            v_adj       = double((v_adj + v_adj') > 0);
            obj.v_adj   = v_adj;
            
        end 
        function plot_mesh(obj)
            figure()
            patch('Faces',obj.faces,'Vertices',obj.vertices,'EdgeColor','black','FaceColor','none') 
        end    
        function faces_area = get_faces_area(obj) 
            v1 = obj.vertices(obj.faces(:,1),:);
            v2 = obj.vertices(obj.faces(:,2),:);
            v3 = obj.vertices(obj.faces(:,3),:);
            faces_area = calc_triangles_area(v1, v2, v3);
        end     
        function vertices_area = get_vertices_area(obj) 
            faces_area    = obj.get_faces_area();
            vf_area_adj   = obj.vf_adj.*faces_area';
            vertices_area = sum(vf_area_adj,2)/3;
        end
        function visualize(obj, fun, type) 
            if strcmp(type,"faces")
                figure()
                patch('Faces',obj.faces,'Vertices',obj.vertices,'FaceVertexCData',full(fun),'FaceColor','flat');
                colorbar
                if min(fun)~=max(fun)
                    caxis([min(fun) max(fun)])
                end
            elseif strcmp(type,"vertices")
                figure()
                patch('Faces',obj.faces,'Vertices',obj.vertices,'FaceVertexCData',full(fun),'EdgeColor','interp','FaceColor','none');
                colorbar
                if min(fun)~=max(fun)
                    caxis([min(fun) max(fun)])
                end
            else
                fprintf('>> unvalid type. type should be faces/vertices\n')
            end
    
        end
        function Iftov = iterp_face_to_vertex(obj)
            
            Af = obj.get_faces_area();
            Av = obj.get_vertices_area();
            
            Iftov = 1/3*obj.vf_adj.*Af'./Av;
            
        end
        function Ivtof = iterp_vertex_to_face(obj)
            Af = obj.get_faces_area();
            Av = obj.get_vertices_area();
            Iftov = obj.iterp_face_to_vertex();
            Ivtof = 1./Af.*Iftov'.*Av';
            
        end
        function [boundary_edges, boundary_vertices] = get_boundary_edges(obj)

            faces_near_edges = obj.vf_adj*obj.vf_adj';        
            boundary_edges = length(find(faces_near_edges==1));
            
            [v, ~] = find(faces_near_edges==1);
            boundary_vertices = unique(v);
            
            
        end
        function genus = get_genus(obj)
            
            boundary_edges = obj.get_boundary_edges();
            fprintf('>> Found %u boundary edges\n', boundary_edges)
            V = size(obj.vertices,1);
            F = size(obj.faces,1);
            E = (length(find(obj.v_adj == 1)) - V)/2;
    
            chi = V - E + F;
            fprintf('>> Chi is %u \n', chi)
            genus = 0.5*(2 - boundary_edges - chi);
        end
        function avg_edge_length = get_avg_edge_length(obj)
           
            [v1,v2] = find(tril(obj.v_adj,-1) == 1);
            edge_vectors = obj.vertices(v1,:) - obj.vertices(v2,:);
            lengths = vecnorm(edge_vectors,2,2);
            avg_edge_length = mean(lengths);
            
        end
        function valence = calc_valence(obj)
            valence = full(sum(obj.v_adj,2) - 1);
        end
        function N = calc_normal(obj)
            
            % Vertices of faces
            v1 = obj.vertices(obj.faces(:,1),:);
            v2 = obj.vertices(obj.faces(:,2),:);
            v3 = obj.vertices(obj.faces(:,3),:);
    
            % Calc normal
            vec1 = v1 - v2;
            vec2 = v1 - v3;

            N = cross(vec1, vec2, 2);
            N = N./vecnorm(N,2,2);
            Nx = N(:,1);
            Ny = N(:,2);
            Nz = N(:,3);

            % Calc centers of faces
            X = (v1(:,1) + v2(:,1) + v3(:,1))/3;
            Y = (v1(:,2) + v2(:,2) + v3(:,2))/3;
            Z = (v1(:,3) + v2(:,3) + v3(:,3))/3;
            
            % Plot Normal
            figure()
            quiver3(X, Y, Z, Nx, Ny, Nz, 'color', '#0072BD')
            hold on
            patch('Faces',obj.faces,'Vertices',obj.vertices,'EdgeColor','black','FaceColor','#4DBEEE') 
            axis equal
 
        end   
        function G = calc_gauss_curv(obj)
            
            v1 = obj.vertices(obj.faces(:,1),:);
            v2 = obj.vertices(obj.faces(:,2),:);
            v3 = obj.vertices(obj.faces(:,3),:);

            % Calc thetas 
            % For v1
            u12 = v2 - v1; w13 = v3 - v1;
            cos_theta1 = max(min(dot(u12,w13,2)./(vecnorm(u12,2,2).*vecnorm(w13,2,2)),1),-1);
            theta1 = real(acos(cos_theta1));

            % For v2
            u23 = v3 - v2; w21 = v1 - v2;
            cos_theta2 = max(min(dot(u23,w21,2)./(vecnorm(u23,2,2).*vecnorm(w21,2,2)),1),-1);
            theta2 = real(acos(cos_theta2));

            % For v3
            u32 = v2 - v3; w31 = v1 - v3;
            cos_theta3 = max(min(dot(u32,w31,2)./(vecnorm(u32,2,2).*vecnorm(w31,2,2)),1),-1);
            theta3 = real(acos(cos_theta3));

            thetas = [theta1 theta2 theta3];

            % construct vertex-face adjacency with thetas instead of '1's
            [i,j] = find(obj.vf_adj);
            facesj = obj.faces(j,:);
            [r,c] = find(i == facesj);
            % i(r) = facesj(r,c);

            ind_of_thetas = sub2ind(size(facesj),r,c);
            thetasj = thetas(j,:);
            sorted_thetas = thetasj(ind_of_thetas);

            ind_to_replace = sub2ind(size(obj.vf_adj),i(r),j(r));
            vtheta_adj = obj.vf_adj;
            vtheta_adj(ind_to_replace) = sorted_thetas;

            Av = obj.get_vertices_area();
            G = (2*pi - sum(vtheta_adj,2))./Av;

        end
    end
end
