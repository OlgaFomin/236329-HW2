function area = compute_triangle_area(vertex_array, traingles_array)

v1 = vertex_array(traingles_array(:,1));
v2 = vertex_array(traingles_array(:,2));
v3 = vertex_array(traingles_array(:,3));

u = v1 - v2;
v = v1 - v3;

CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
ThetaInDegrees = real(acosd(CosTheta));

area = 0.5 * vecnorm(u,2,2) * vecnorm(v2,2) * sind(ThetaInDegrees); 

end