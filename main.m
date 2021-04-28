off_filepath = 'hw2_data/bunny2.off';

%1
[vertices, triangles] = read_off(off_filepath);

%2
obj = meshHandler(vertices, triangles);

%3
plot_obj(vertices, triangles);

%4
traingles_area = compute_triangle_area(vertices, triangles);
vertex_area = compute_vertex_area(vertices, triangles);

%5
visualize_fun(vertices, triangles, traingles_area, "face")
visualize_fun(vertices, triangles, vertex_area, "vertex")




