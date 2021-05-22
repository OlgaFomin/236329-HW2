function triangles_area = calc_triangles_area(v1, v2, v3)

    u = v1 - v2;
    w = v1 - v3;

    cos_theta = max(min(dot(u,w,2)./(vecnorm(u,2,2).*vecnorm(w,2,2)),1),-1);
    theta_degrees = real(acosd(cos_theta));

    triangles_area = 0.5 * vecnorm(u,2,2) .* vecnorm(w,2,2) .* sind(theta_degrees); 

end