%% Digital Geometry Processing - 236329
% HW 2
% Olga Fomin 316948694
% Asaf Levy 201631349
clear all; close all; clc


%% Objectives

% 1
off_filepath = 'hw2_data/sphere_s0.off';
[vertices, faces] = read_off(off_filepath);

out_filepath = 'test_write_off.off';
write_off(vertices, faces, out_filepath);

% 2
obj = mesh(vertices, faces);

% 3
fprintf('>> Plotting mesh\n')
obj.plot_mesh();

% 4
fprintf('>> Calculating faces area\n')
faces_area = obj.get_faces_area();
fprintf('>> Calculating vertices area\n')
vertices_area = obj.get_vertices_area();

% 5
fprintf('>> Plotting function on faces\n')
obj.visualize(faces_area, 'faces');
fprintf('>> Plotting function on vertices\n')
obj.visualize(vertices_area, 'vertices');

% 6
fprintf('>> Calculating Interpulation from faces to vertices\n')
Iftov = obj.iterp_face_to_vertex();
fprintf('>> Calculating number of boundary edges\n')
Ivtof = obj.iterp_vertex_to_face();

% 7
fprintf('>> Calculating number of boundary edges\n')
[boundary_edges, boundary_vertices] = obj.get_boundary_edges();
fprintf('>> Calculating genus\n')
genus = obj.get_genus();

%% Analysis

% 1
[vertices_s0, faces_s0] = read_off('hw2_data/sphere_s0.off'); shpere0 = mesh(vertices_s0, faces_s0); shpere0.plot_mesh();
Err0 = tri_error_unit_sphere(shpere0);
avg_length0 = shpere0.get_avg_edge_length();

[vertices_s1, faces_s1] = read_off('hw2_data/sphere_s1.off'); shpere1 = mesh(vertices_s1, faces_s1); 
Err1 = tri_error_unit_sphere(shpere1);
avg_length1 = shpere1.get_avg_edge_length();

[vertices_s2, faces_s2] = read_off('hw2_data/sphere_s2.off'); shpere2 = mesh(vertices_s2, faces_s2); shpere2.plot_mesh();
Err2 = tri_error_unit_sphere(shpere2);
avg_length2 = shpere2.get_avg_edge_length();

[vertices_s3, faces_s3] = read_off('hw2_data/sphere_s3.off'); shpere3 = mesh(vertices_s3, faces_s3);
Err3 = tri_error_unit_sphere(shpere3);
avg_length3 = shpere3.get_avg_edge_length();

[vertices_s4, faces_s4] = read_off('hw2_data/sphere_s4.off'); shpere4 = mesh(vertices_s4, faces_s4);
Err4 = tri_error_unit_sphere(shpere4);
avg_length4 = shpere4.get_avg_edge_length();

[vertices_s5, faces_s5] = read_off('hw2_data/sphere_s5.off'); shpere5 = mesh(vertices_s5, faces_s5);
Err5 = tri_error_unit_sphere(shpere5);
avg_length5 = shpere5.get_avg_edge_length();


Err = [Err0 Err1 Err2 Err3 Err4 Err5];
avg_length = [avg_length0 avg_length1 avg_length2 avg_length3 avg_length4 avg_length5];
figure()
plot(avg_length, Err, 'r')
xlabel('AVG edge length');
ylabel('Err');
title('Triangulation error VS AVG edge length');

% 2
off_filepath = 'hw2_data/disk.off';[vertices, faces] = read_off(off_filepath);obj = mesh(vertices, faces);
valence = obj.calc_valence();obj.visualize(valence, 'vertices');
fprintf('>> AVG valence for %s is %f\n', off_filepath, mean(valence))

% 3
off_filepath = 'hw2_data/sphere_s0.off';[vertices, faces] = read_off(off_filepath); obj = mesh(vertices, faces);
Iftov = obj.iterp_face_to_vertex();
Ivtof = obj.iterp_vertex_to_face();
kernel_ftov = null(full(Iftov));
kernel_vtof = null(full(Ivtof));

% 4 - Calc and plot normals
off_filepath = 'hw2_data/sphere_s0.off';[vertices, faces] = read_off(off_filepath);obj = mesh(vertices, faces);
N = obj.calc_normal();

off_filepath = 'hw2_data/torus_fat_r2.off';[vertices, faces] = read_off(off_filepath);obj = mesh(vertices, faces);
N = obj.calc_normal();

off_filepath = 'hw2_data/disk.off';[vertices, faces] = read_off(off_filepath);obj = mesh(vertices, faces);
N = obj.calc_normal();

% 5 -  Calc and plot gaussian curvature
off_filepath = 'hw2_data/disk.off';[vertices, faces] = read_off(off_filepath);obj = mesh(vertices, faces);obj.plot_mesh();
G = obj.calc_gauss_curv();
obj.visualize(G, 'vertices');






