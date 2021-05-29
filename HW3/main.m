%% Digital Geometry Processing - 236329
% HW 3
% Olga Fomin 316948694
% Asaf Levy 201631349
clear all; close all; clc


%% Objectives

% 1
off_filepath = 'hw2_data/sphere_s1.off';
[vertices, faces] = read_off(off_filepath);
mesh_sphere = MeshHandle(vertices, faces);


[Nf,coord] = mesh_sphere.face_normal();
mesh_sphere.visualize_vec(coord,Nf,0,'none');

faces_area = mesh_sphere.get_faces_area();
mesh_sphere.visualize_vec(coord,Nf,faces_area,'faces');

vertices_area = mesh_sphere.get_vertices_area();
mesh_sphere.visualize_vec(coord,Nf,vertices_area,'vertices');

% 3
[Nv,coord] = mesh_sphere.vertex_normal();
mesh_sphere.visualize_vec(coord,Nv,0,'none');


% 4
G = mesh_sphere.calc_gauss_curv();
mesh_sphere.visualize_fun(G, 'vertices');

H = mesh_sphere.calc_mean_curv();
mesh_sphere.visualize_fun(H, 'vertices');

% 2a
% Gradiant
vertices = mesh_sphere.vertices;
faces = mesh_sphere.faces;

v1 = vertices(faces(:,1),:);
v2 = vertices(faces(:,2),:);
v3 = vertices(faces(:,3),:);

xc = (v1(:,1) + v2(:,1) + v3(:,1))/3;
yc = (v1(:,2) + v2(:,2) + v3(:,2))/3;
zc = (v1(:,3) + v2(:,3) + v3(:,3))/3;

grad = mesh_sphere.calc_grad();
f = vertices(:,3);
gradf = grad*f;
F = size(mesh_sphere.faces,1);
gradf = [gradf(1:F) gradf(F+1:2*F) gradf(2*F+1:3*F)];

coord = [xc, yc, zc];
mesh_sphere.visualize_vec(coord,gradf./vecnorm(gradf,2,2),0,'none');

% Divergence
div = mesh_sphere.calc_div();
divf = div*[gradf(:,1);gradf(:,2);gradf(:,3)];
mesh_sphere.visualize_fun(divf, 'vertices');

% Laplacian
L = mesh_sphere.calc_laplas();
mesh_sphere.visualize_fun(diag(L), 'vertices');

% 2b
% Comparing Laplace operator with cot-laplce implementation:
cotL = mesh_sphere.calc_laplas_cot();
L = mesh_sphere.calc_laplas();
assert(abs(max(max(cotL - L)) ) < 1e-8) 


%% Analysis

[vertices, faces] = read_off('hw2_data/disk.off'); mesh_disk = MeshHandle(vertices, faces); 
analyze_W(mesh_disk, 'Disk', 'vertices height', mesh_disk.vertices(:,3), 'gaussian curvature', mesh_disk.calc_gauss_curv(), true)

[vertices, faces] = read_off('hw2_data/sphere_s0.off'); mesh_sphere = MeshHandle(vertices, faces); 
analyze_W(mesh_sphere, 'Sphere', '"1"s vector', ones(size(mesh_sphere.vertices,1) , 1), 'valence', mesh_sphere.calc_valence(), false)

[vertices, faces] = read_off('hw2_data/vase.off'); mesh_vase = MeshHandle(vertices, faces); 
analyze_W(mesh_vase, 'Vase', 'valence', mesh_vase.calc_valence(), 'vertices area', mesh_vase.get_vertices_area(), false)









