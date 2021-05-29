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

% 1 
% 1.1 Null
off_filepath = 'hw2_data/disk.off'; [vertices, faces] = read_off(off_filepath); mesh_disk = MeshHandle(vertices, faces); W_disk = mesh_disk.calc_W();
fdisk_const = mesh_disk.vertices(:,3); % Vertices height
fdisk_nonconst = mesh_disk.calc_gauss_curv(); % Gaussian curvature
const_null_disk = vecnorm(W_disk*fdisk_const); nonconst_null_disk = vecnorm(W_disk*fdisk_nonconst);
disp(['Disk - ||W*f|| for constant f (vertices height) = ', num2str(const_null_disk)])
disp(['Disk - ||W*f|| for non constant f (gaussian curvature) = ', num2str(nonconst_null_disk)])
mesh_disk.visualize_fun(fdisk_const, 'vertices');hold on;title(['f = height, ||W*f|| = ', num2str(const_null_disk)]);hold off
mesh_disk.visualize_fun(fdisk_nonconst, 'vertices');hold on;title(['f = gaussian curv, ||W*f|| = ', num2str(nonconst_null_disk)]);hold off



off_filepath = 'hw2_data/sphere_s0.off'; [vertices, faces] = read_off(off_filepath); mesh_sphere = MeshHandle(vertices, faces); W_sphere = mesh_sphere.calc_W();
fsphere_const = ones(size(mesh_sphere.vertices,1) , 1); % '1's vector
fsphere_nonconst = mesh_sphere.calc_valence(); % Valence
const_null_sphere = vecnorm(W_sphere*fsphere_const); nonconst_null_sphere = vecnorm(W_sphere*fsphere_nonconst);
disp(['Sphere - ||W*f|| for constant f ("1"s vector) = ', num2str(const_null_sphere)])
disp(['Sphere - ||W*f|| for non constant f (valence) = ', num2str(nonconst_null_sphere)])
mesh_sphere.visualize_fun(fsphere_const, 'vertices');hold on;title(['f = 1s vector, ||W*f|| = ', num2str(const_null_sphere)]);hold off
mesh_sphere.visualize_fun(fsphere_nonconst, 'vertices');hold on;title(['f = valence, ||W*f|| = ', num2str(nonconst_null_sphere)]);hold off



off_filepath = 'hw2_data/vase.off'; [vertices, faces] = read_off(off_filepath); mesh_vase = MeshHandle(vertices, faces); W_vase = mesh_vase.calc_W();
fvase_const = mesh_vase.calc_valence(); % Valence
fvase_nonconst = mesh_vase.get_vertices_area(); % Vertices area
const_null_vase = vecnorm(W_vase*fvase_const); nonconst_null_vase = vecnorm(W_vase*fvase_nonconst);
disp(['Vase - ||W*f|| for constant f (valence) = ', num2str(const_null_vase)])
disp(['Vase - ||W*f|| for non constant f (vertices area) = ', num2str(nonconst_null_vase)])
mesh_vase.visualize_fun(fvase_const, 'vertices');hold on;title(['f = valence, ||W*f|| = ', num2str(const_null_vase)]);hold off
mesh_vase.visualize_fun(fvase_nonconst, 'vertices');hold on;title(['f = v-area, ||W*f|| = ', num2str(nonconst_null_vase)]);hold off



% 1.2 - Symmetry
disp(['Disk Symmetry - || W - WT ||F =  ', num2str(norm(W_disk - transpose(W_disk),'fro'))])
disp(['Sphere Symmetry - || W - WT ||F =  ', num2str(norm(W_sphere - transpose(W_sphere),'fro'))])
disp(['Vase Symmetry - || W - WT ||F =  ', num2str(norm(W_vase - transpose(W_vase),'fro'))])



% 1.3 - Localization
disp(['Disk - nnz = ', num2str(nnz(W_disk)), ' #edges = ', num2str(size(mesh_disk.get_edges(),1)) , ' #vertices = ', num2str(size(mesh_disk.vertices,1))])
disp(['Sphere - nnz = ', num2str(nnz(W_sphere)), ' #edges = ', num2str(size(mesh_sphere.get_edges(),1)) , ' #vertices = ', num2str(size(mesh_sphere.vertices,1))])
disp(['Vase - nnz = ', num2str(nnz(W_vase)), ' #edges = ', num2str(size(mesh_vase.get_edges(),1)) , ' #vertices = ', num2str(size(mesh_vase.vertices,1))])


% 1.4 - Positivity
pos_disk = find(W_disk < -10^-10);
if ~isempty(pos_disk)
    disp('W_disk is not positive')
else
    disp('W_disk is positive')
end


pos_sphere = find(W_sphere < -10^-10);
if ~isempty(pos_sphere)
    disp('W_sphere is not positive')
else
    disp('W_sphere is positive')
end


pos_vase = find(W_vase < -10^-10);
if ~isempty(pos_vase)
    disp('W_vase is not positive')
else
    disp('W_vase is positive')
end


% 1.5 - Positive semi-definite
seig_disk = eigs(W_disk, 1, 'sm'); leig_disk = eigs(W_disk, 1, 'lm');
if seig_disk >= 0 && leig_disk >= 0
    disp (['Disk - Positive semi-definite. Smallest eigenvalue = ', num2str(seig_disk), ', Largest eigenvalue = ', num2str(leig_disk)])
else
    disp('Disk - Non positive semi-definite.');
end
    
seig_sphere = eigs(W_sphere, 1, 'sm'); leig_sphere = eigs(W_sphere, 1, 'lm');
if seig_sphere >= 0 && leig_sphere >= 0
    disp (['Sphere - Positive semi-definite. Smallest eigenvalue = ', num2str(seig_sphere), ', Largest eigenvalue = ', num2str(leig_sphere)])
else
    disp('Sphere - Non positive semi-definite.');
end
    
seig_vase = eigs(W_vase, 1, 'sm'); leig_vase = eigs(W_vase, 1, 'lm');
if seig_vase >= 0 && leig_vase >= 0
    disp (['Vase - Positive semi-definite. Smallest eigenvalue = ', num2str(seig_vase), ', Largest eigenvalue = ', num2str(leig_vase)])
else
    disp('Vase - Non positive semi-definite.');
end
    

% 2
% Visualize 9 eigenvectors
k = 9;
[B_disk, eig_disk] = eigs(W_disk, k, 'sm');
figure()
for i = 1:k
    subplot(3, 3, i);
    patch('Faces',mesh_disk.faces,'Vertices',mesh_disk.vertices,'FaceVertexCData',full(B_disk(:,i)),'EdgeColor','interp','FaceColor','none');
    colorbar
    title(['Eigenvalue = ', num2str(eig_disk(i, i))]);
end


k = 9;
[B_sphere, eig_sphere] = eigs(W_sphere, k, 'sm');
figure()
for i = 1:k
    subplot(3, 3, i);
    patch('Faces',mesh_sphere.faces,'Vertices',mesh_sphere.vertices,'FaceVertexCData',full(B_sphere(:,i)),'EdgeColor','interp','FaceColor','none');
    colorbar
    title(['Eigenvalue = ', num2str(eig_sphere(i, i))]);
end

k = 9;
[B_vase, eig_vase] = eigs(W_vase, k, 'sm');
figure()
for i = 1:k
    subplot(3, 3, i);
    patch('Faces',mesh_vase.faces,'Vertices',mesh_vase.vertices,'FaceVertexCData',full(B_vase(:,i)),'EdgeColor','interp','FaceColor','none');
    colorbar
    title(['Eigenvalue = ', num2str(eig_vase(i, i))]);
end



% Visualize approximations
fhat = zeros(size(mesh_disk.vertices,1),1); fhat(size(mesh_disk.vertices,1)) = 1; f = fhat;
mesh_disk.visualize_fun(fhat, 'vertices'); hold on; title('fhat'); hold off

ks = 10:20:300;
[B, ~] = eigs(W_disk, 300, 'sm');
norms = zeros(length(ks), 1);
for i = 1:length(ks)
    k = ks(i);
    Bi = B(:, 1:k);
    gi = Bi * (transpose(Bi) * f);
    norms(i) = vecnorm(gi - f);
    if k == 70 || k == 250
        mesh_disk.visualize_fun(gi, 'vertices'); hold on; title(['k = ', num2str(k)]); hold off
    end
   
end

% Visualize error of approximation as a function of k
figure();
plot(ks, norms);
title('||gi - f|| as a function of k');
xlabel('ki');
ylabel('||gi - f||');









