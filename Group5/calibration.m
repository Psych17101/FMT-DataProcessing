clc; clear variables; close all; 

%%  Config
urange = 0:2:20; 
fileName = 'Calibration_%03d'; 

for ii = 1:numel(urange)
    fileIn = sprintf(fileName,urange(ii))
    delimiter = ' ';
    startRow = 23;
    formatSpec = '%s';
    try
        fileID = fopen(fileIn,'r');
    catch
        fileID = fopen(fileIn{1},'r');
    end
    tmp = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    tmp2 = strrep(tmp{1},',','.');
    tmp2 = str2mat(tmp2);
    vel = tmp2(:,10:17);
    u(ii) = mean(str2num(vel));
    
end