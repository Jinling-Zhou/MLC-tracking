% Read out DICOM plan file into individual data files. 
% Jinling Zhou, 3/25/2023.

% specify the path to the DICOM RP file
path = 'RP.QA303005.test.dcm';

% read the DICOM file
dicomInfo = dicominfo(path);

% call the recursive function to traverse the DICOM metadata structure
outputPath = 'RP.Dicom';
mkdir(outputPath);
traverseMetadata(dicomInfo, outputPath);

% recursive function to traverse the DICOM metadata structure
function traverseMetadata(metadata, outputPath)
    % get the field names of the metadata structure
    fields = fieldnames(metadata);
    
    % loop through each field
    for i = 1:numel(fields)
        % get the field name and value
        fieldName = fields{i};
        fieldValue = metadata.(fieldName);
        
        % check if the field value is a numerical or string type
        if isnumeric(fieldValue) || ischar(fieldValue)
            % if the field value is a numerical or string type, write the full information to an ASCII file
            outputFilePath = fullfile(outputPath, [fieldName '.txt']);
            fileID = fopen(outputFilePath, 'w');
            fprintf(fileID, '%s\n', fieldName);
         %   fprintf(fileID, '%s\n', num2str(fieldValue));
            for j = 1:length(fieldValue)
                  fprintf(fileID, '%s\n', num2str(fieldValue(j)));
            end

            fclose(fileID);
        else
            % if the field value is not a numerical or string type, recursively traverse the field
            newOutputPath = fullfile(outputPath, fieldName);
            mkdir(newOutputPath);
            traverseMetadata(fieldValue, newOutputPath);
        end
    end
end


