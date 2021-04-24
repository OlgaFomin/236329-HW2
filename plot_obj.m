function plot_obj(vertices, triangles)

    showVertices = 0;
    disp '>> Plot beginning';
    figure(); grid on; grid minor; axis equal;
    hold on;
    if(showVertices)
        for i=1:size(vertices,2)
            currV = vertices(i);
            plot3(currV(1),currV(2),currV(3),'ob');

        end
    end
    for j=1:size(triangles,1)
        
        currF = triangles(j,:);
        
        xCoo = [];
        yCoo = [];
        zCoo = [];
        for k=1:size(currF,2)
            xCoo = [xCoo vertices(currF(k)+1,1)];
            yCoo = [yCoo vertices(currF(k)+1,2)];
            zCoo = [zCoo vertices(currF(k)+1,3)];
        end
        xCoo = [xCoo xCoo(1)];
        yCoo = [yCoo yCoo(1)];
        zCoo = [zCoo zCoo(1)];

        plot3(xCoo, zCoo, yCoo); %<----

    end
    view(-140,12);
    hold off;
    
end