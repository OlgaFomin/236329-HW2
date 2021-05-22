%% Digital Geometry Processing - 236329
% HW 3
% Olga Fomin 316948694
% Asaf Levy 201631349
clear all; close all; clc


%% Objectives

% 1
off_filepath = 'hw2_data/torus_fat_r2.off';
[vertices, faces] = read_off(off_filepath);
obj = mesh(vertices, faces);


[Nf,coord] = obj.face_normal();
obj.visualize_vec(coord,Nf,0,'none');

faces_area = obj.get_faces_area();
obj.visualize_vec(coord,Nf,faces_area,'faces');

vertices_area = obj.get_vertices_area();
obj.visualize_vec(coord,Nf,vertices_area,'vertices');

% 3
[Nv,coord] = obj.vertex_normal();
obj.visualize_vec(coord,Nv,0,'none');


% 4
G = obj.calc_gauss_curv();
obj.visualize_fun(G, 'vertices');

H = obj.calc_mean_curv();
obj.visualize_fun(H, 'vertices');

% 2a
% Gradiant
gradf = calc_gradf(obj, H);
vertices = obj.vertices;
faces = obj.faces;

v1 = vertices(faces(:,1),:);
v2 = vertices(faces(:,2),:);
v3 = vertices(faces(:,3),:);

xc = (v1(:,1) + v2(:,1) + v3(:,1))/3;
yc = (v1(:,2) + v2(:,2) + v3(:,2))/3;
zc = (v1(:,3) + v2(:,3) + v3(:,3))/3;

coord = [xc, yc, zc];
obj.visualize_vec(coord,gradf./vecnorm(gradf,2,2),G,'vertices');

% Divergence
Fx = gradf(:,1);
Fy = gradf(:,2);
Fz = gradf(:,3);
X = xc;
Y = yc;
Z = zc;

div = divergence(X,Y,Z,Fx,Fy,Fz);


