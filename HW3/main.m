%% Digital Geometry Processing - 236329
% HW 3
% Olga Fomin 316948694
% Asaf Levy 201631349
clear all; close all; clc


%% Objectives

% 1
off_filepath = 'hw2_data/sphere_s0.off';
[vertices, faces] = read_off(off_filepath);
obj = mesh(vertices, faces);


[Nf,coord] = obj.face_normal();
obj.visualize_vec(Nf,coord,0,'none');

faces_area = obj.get_faces_area();
obj.visualize_vec(Nf,coord,faces_area,'faces');

vertices_area = obj.get_vertices_area();
obj.visualize_vec(Nf,coord,vertices_area,'vertices');



% 2
% NEED TO IMPLEMENT!!!!!!!


% 3
[Nv,coord] = obj.vertex_normal();
obj.visualize_vec(Nv,coord,0,'none');


% 4
G = obj.calc_gauss_curv();
obj.visualize(G, 'vertices');








