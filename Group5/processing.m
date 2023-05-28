clc; clear; close all;

%%

fileIn = 'CorrelationTest';

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
t = tmp2(:,1:8);
t = str2num(t);
u= tmp2(:,10:17);
u = str2num(u);
u_new = u - mean(u);

C = xcorr(u, u);
C_new = xcorr(u_new, u_new);
y = ifft(C_new);
figure()
plot(y)

figure(2)
t = linspace(0,10,length(C_new));
plot(t,C_new)


C = xcorr(u, u);
C_new = C/max(C);
figure
plot(C_new)