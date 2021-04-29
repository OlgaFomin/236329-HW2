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

%7
[genus, boundary_edges] = calc_genus(obj.adjacency_matrix, triangles, vertices);

% Analysis - 1
off_filepath = 'hw2_data/sphere_s0.off';
[vertices, triangles] = read_off(off_filepath);
plot_obj(vertices, triangles);

off_filepath = 'hw2_data/sphere_s2.off';
[vertices, triangles] = read_off(off_filepath);
plot_obj(vertices, triangles);

off_filepath = 'hw2_data/sphere_s4.off';
[vertices, triangles] = read_off(off_filepath);
plot_obj(vertices, triangles);


