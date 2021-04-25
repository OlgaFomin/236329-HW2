off_filepath = 'hw2_data/bunny2.off';
[vertices, triangles] = read_off(off_filepath);
plot_obj(vertices, triangles);

obj = meshHandler(vertices, triangles);
traingles_area = compute_triangle_area(vertices, triangles);
vertex_area = compute_vertex_area(vertices, triangles);