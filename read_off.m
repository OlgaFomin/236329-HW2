function [vertices, triangles] = read_off(off_filepath)

    fileID = fopen(off_filepath);
      
    % disp '>> Storing the structure...';
    fprintf('>> Storing the structure...')
    dimRow = true;
    rowCount = 1;
    vertices = [];
    triangles = [];
    
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
                    vertices = [vertices; splittedRow(1) splittedRow(2) splittedRow(3)];
                end
                
                if(vertexCount < rowCount && (rowCount-vertexCount) <= faceCount)  
                    triangles = [triangles; splittedRow(2) + 1, splittedRow(3) + 1, splittedRow(4) + 1];
                end

                rowCount = rowCount +1;

                % progress
                if(mod(rowCount,1000)==0)
                    fprintf('.');
                end
            end
        end
    end
    
    fprintf('>> Strucure stored');
end