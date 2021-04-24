function plot_obj(obj)

    vertexList = obj.vertices;
    faceList = obj.faces;
    
    showVertices = 0;
    disp '>> Plot beginning';
    figure(); grid on; grid minor; axis equal;
    hold on;
    if(showVertices)
        for i=1:size(vertexList,2)
            %currV = vertexList{i};
            currV = vertexList(i);
            plot3(currV(1),currV(2),currV(3),'ob');

        end
    end
    %for j=1:size(faceList,2)
    for j=1:size(faceList,1)
        
        %currF = faceList{j};
        currF = faceList(j,:);
        
        xCoo = [];
        yCoo = [];
        zCoo = [];
        for k=1:size(currF,2)
%             xCoo = [xCoo vertexList{currF(k)+1}(1)];
%             yCoo = [yCoo vertexList{currF(k)+1}(2)];
%             zCoo = [zCoo vertexList{currF(k)+1}(3)];
            xCoo = [xCoo vertexList(currF(k)+1,1)];
            yCoo = [yCoo vertexList(currF(k)+1,2)];
            zCoo = [zCoo vertexList(currF(k)+1,3)];
        end
        xCoo = [xCoo xCoo(1)];
        yCoo = [yCoo yCoo(1)];
        zCoo = [zCoo zCoo(1)];

        plot3(xCoo, zCoo, yCoo); %<----

    end
    view(-140,12);
    hold off;
    
end