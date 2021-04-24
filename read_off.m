function [obj] = read_off(off_filepath)

    fileID = fopen(off_filepath);
      
    disp '>> Storing the structure...';
    dimRow = true;
    rowCount = 1;
    %vertexList{1} = [];
    %faceList{1} = [];
    vertexList = [];
    faceList = [];
    
    while (~feof(fileID))

        currRow = textscan(fileID,'%s',1,'Delimiter','\n');
        splittedRow = strsplit(char(currRow{1}),' ');

        if (~strcmp(splittedRow(1),'#') && ~strcmp(splittedRow(1),'OFF'))
            splittedRow = str2double(splittedRow);

            if(dimRow)
                dimRow = false;
                vertexCount = splittedRow(1);
                faceCount = splittedRow(2);
                edgeCount = splittedRow(3);
            else
                if(rowCount <= vertexCount)
                    %vertexList{rowCount} = [splittedRow(1) splittedRow(2) splittedRow(3)];
                    vertexList = [vertexList; splittedRow(1) splittedRow(2) splittedRow(3)];
                end
                
                if(vertexCount < rowCount && (rowCount-vertexCount) <= faceCount)  
                    if(splittedRow(1) == 3)
                        %faceList{rowCount-vertexCount} = [splittedRow(2) splittedRow(3) splittedRow(4)];
                        faceList = [faceList; splittedRow(2) splittedRow(3) splittedRow(4)];
                    end
                    
                    if(splittedRow(1) == 4)
                        %faceList{rowCount-vertexCount} = [splittedRow(2) splittedRow(3) splittedRow(4) splittedRow(5)];
                        faceList = [faceList; splittedRow(2) splittedRow(3) splittedRow(4) splittedRow(5)];
                    end
                end

                rowCount = rowCount +1;

                % progress
                if(mod(rowCount,1000)==0)
                    disp('.');
                end
            end
        end
    end
    
    obj.vertices = vertexList;
    obj.faces = faceList;
    obj.vertexCount = vertexCount;
    obj.faceCount = faceCount;
    obj.edgeCount = edgeCount;
    
    disp '>> Strucure stored';
end