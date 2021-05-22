function Err = tri_error_unit_sphere(obj)

    vertices = obj.vertices;
    faces = obj.faces;
    
    v1 = vertices(faces(:,1),:);
    v2 = vertices(faces(:,2),:);
    v3 = vertices(faces(:,3),:);
    
    Xx = (v1(:,1) + v2(:,1) + v3(:,1))/3;
    Xy = (v1(:,2) + v2(:,2) + v3(:,2))/3;
    Xz = (v1(:,3) + v2(:,3) + v3(:,3))/3;
    
    X = [Xx Xy Xz];
    S = X./vecnorm(X,2,2);
    
    Err = 1/size(faces,1)*sum(vecnorm(X - S,2,2));
    
end