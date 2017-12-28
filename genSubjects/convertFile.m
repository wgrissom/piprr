% Convert the *.mat file into text files.
load('normal_head_dist');

% %% Duke
% subStr = 'Duke';
% Str1 = 'facs';
% fileID_X = fopen('DukeX.txt','w');
% fileID_Y = fopen('DukeY.txt','w');
% fileID_Z = fopen('DukeZ.txt','w');
% for n = 1:length(DukeAPfacs)
%     fprintf(fileID_X,'%f\r\n',round(DukeLRfacs(n),4));
%     fprintf(fileID_Y,'%f\r\n',round(DukeAPfacs(n),4));
%     fprintf(fileID_Z,'%f\r\n',round(DukeHFfacs(n),4));
% end
% fclose(fileID_X);
% fclose(fileID_Y);
% fclose(fileID_Z);

%% Ella
subStr = 'Ella';
Str1 = 'facs';
fileID_X = fopen('EllaX.txt','w');
fileID_Y = fopen('EllaY.txt','w');
fileID_Z = fopen('EllaZ.txt','w');
for n = 1:length(DukeAPfacs)
    fprintf(fileID_X,'%f\r\n',round(EllaLRfacs(n),4));
    fprintf(fileID_Y,'%f\r\n',round(EllaAPfacs(n),4));
    fprintf(fileID_Z,'%f\r\n',round(EllaHFfacs(n),4));
end
fclose(fileID_X);
fclose(fileID_Y);
fclose(fileID_Z);
