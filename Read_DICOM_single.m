% Read out a DICOM plan file to a single ASCII file (.txt).
% Jinling Zhou, 3/27/2023.
clc;
clear all;
close all;

% specify the path to the DICOM RP file
path ='RP.QA303005.test.dcm';

% read the DICOM file
dicomInfo = dicominfo(path);

% call the recursive function to traverse the DICOM metadata structure
% outputPath = fullfile('Single', 'RP.txt');
fileID = fopen('RP.txt', 'a');
traverseMetadata(dicomInfo, fileID);
fclose(fileID);

% recursive function to traverse the DICOM metadata structure
function traverseMetadata(metadata, fID)
    % get the field names of the metadata structure
    fields = fieldnames(metadata);
    % fID = fopen('outputPath\RP.Dicom.txt', 'a');

    % loop through each field
    for i = 1:numel(fields)
        % get the field name and value
        fieldName = fields{i};
        fieldValue = metadata.(fieldName);
        
        % check if the field value is a numerical or string type
        if isnumeric(fieldValue) || ischar(fieldValue)
            % if the field value is a numerical or string type, write the full information to an ASCII file 
            fprintf(fID,'\n');
            fprintf(fID, '%s', fieldName);
         %   fprintf(fileID, '%s\n', num2str(fieldValue));
            for j = 1:length(fieldValue)
                  fprintf(fID, ' %s', num2str(fieldValue(j)));
                  fprintf('\n'); % this line might be useless, as it doesn't write in the ASCII file (fID).
            end

        else
            % if the field value is not a numerical or string type, recursively traverse the field
            % newOutputPath = fullfile(outputPath, fieldName);
            % mkdir(newOutputPath);
            traverseMetadata(fieldValue, fID);
        end
    end
end




