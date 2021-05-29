function analyze_W(mesh, obj, const_fun_name, f_const, non_const_fun_name, f_nonconst, approx)
    
    % Calc W
    W = mesh.calc_W();
    
    % 1.1 - Null
    const_null = vecnorm(W*f_const); 
    nonconst_null = vecnorm(W*f_nonconst);
    disp([obj, ' - ||W*f|| for constant f (', const_fun_name, ') = ', num2str(const_null)])
    disp([obj, ' - ||W*f|| for non constant f (', non_const_fun_name, ') = ',  num2str(nonconst_null)])
    mesh.visualize_fun(f_const, 'vertices');hold on;title(['f = ', const_fun_name, ', ||W*f|| = ', num2str(const_null)]);hold off
    mesh.visualize_fun(f_nonconst, 'vertices');hold on;title(['f = ', non_const_fun_name, ', ||W*f|| = ', num2str(nonconst_null)]);hold off

    % 1.2 - Symmetry
    disp([obj, ' Symmetry - || W - WT ||F =  ', num2str(norm(W - transpose(W),'fro'))])
    
    % 1.3 - Localization
    disp([obj, ' - nnz = ', num2str(nnz(W)), ' #edges = ', num2str(size(mesh.get_edges(),1)) , ' #vertices = ', num2str(size(mesh.vertices,1))])

    % 1.4 - Positivity
    pos = find(W < -10^-10);
    if ~isempty(pos)
        disp([obj, ' - W is not positive'])
    else
        disp([obj, ' - W is positive'])
    end
    
    % 1.5 - Positive semi-definite
    seig = eigs(W, 1, 'sm'); 
    leig = eigs(W, 1, 'lm');
    if seig >= 0 && leig >= 0
        disp ([obj, ' - Positive semi-definite. Smallest eigenvalue = ', num2str(seig), ', Largest eigenvalue = ', num2str(leig)])
    else
        disp([obj, ' - Non positive semi-definite.']);
    end
    
    
    % 2- Visualize 9 eigenvectors
    k = 9;
    [Bi, eig] = eigs(W, k, 'sm');
    figure()
    for i = 1:k
        subplot(3, 3, i);
        patch('Faces',mesh.faces,'Vertices',mesh.vertices,'FaceVertexCData',full(Bi(:,i)),'EdgeColor','interp','FaceColor','none');
        colorbar
        title(['Eigenvalue = ', num2str(eig(i, i))]);
    end

    if approx 
        % 2 - Function approximations
        fhat = zeros(size(mesh.vertices,1),1); 
        fhat(size(mesh.vertices,1)) = 1; 
        f = fhat;
        mesh.visualize_fun(f, 'vertices'); hold on; title('fhat'); hold off

        ks = 10:20:300;
        [B, ~] = eigs(W, 300, 'sm');
        norms = zeros(length(ks), 1);
        for i = 1:length(ks)
            k = ks(i);
            Bi = B(:, 1:k);
            gi = Bi * (transpose(Bi) * f);
            norms(i) = vecnorm(gi - f);
            if k == 70 || k == 250
                mesh.visualize_fun(gi, 'vertices'); hold on; title(['k = ', num2str(k)]); hold off
            end
        end


        % Error of approximation as a function of k
        figure();
        plot(ks, norms);
        title('||gi - f|| as a function of k');
        xlabel('k');
        ylabel('||gi - f||');
    end

end