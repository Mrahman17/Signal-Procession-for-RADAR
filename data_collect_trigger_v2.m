clear; clc; close all; 

subject = '13 jan ridvan ademola sean akthar emre';
root_dir='/media/rspl-admin/Seagate Backup Plus Drive/Gait Data/Data_Martelli_10sept/';
tot_dir=dir(root_dir);
%for dd=1:18
data = [root_dir,'*.bin'];
 %RDout = ['D:\Gait Data\Summer 2021 Data\'];
mDout = '/media/rspl-admin/Seagate Backup Plus Drive/Gait Data/Martelli_mD_out_legs/';
% DOAout = ['/mnt/HDD01/rspl-admin/DATASETS/Fall Sequential/Outputs/' subject '/rangeDOA/'];

if ~exist(mDout, 'dir')
       mkdir(mDout)
end
% if ~exist(DOAout, 'dir')
%        mkdir(DOAout)
% end

files = dir(data);

seqPerRecord = 1; 

filenames2 = {files.name};

for z = 1:length(filenames2)
        temp{1,z} = filenames2{z}(1:end-10);
end
uniqs = unique(temp);
for j = 1:length(uniqs)
        match = strfind(filenames2,uniqs{j}); % find matches
        idx = find(~cellfun(@isempty,match)); % find non-empty indices
        RDC = [];
        % concat RDCs with same names
        for r = 1:length(idx)
                fname = fullfile(files(idx(r)).folder,files(idx(r)).name);
                temp2 = RDC_extract(fname);
                RDC = [RDC temp2];
        end
        % divide into sub RDCs
        numChirps = floor(size(RDC,2)/seqPerRecord);
        for r =1:seqPerRecord
                tic
                msg = ['Processing: Subject ''' subject ''', File: ' int2str(j) ' of ' int2str(length(uniqs)) ', Part ' ...
                        num2str(r) '/' num2str(seqPerRecord)];   % loading message
                disp(msg);
                subRDC = RDC(:,(r-1)*numChirps+1:r*numChirps,:);
                mD_Out = [mDout uniqs{j} '_' num2str(r) '.png'];
%                  RD_Out = [RDout uniqs{j} '_' num2str(r) '.avi'];
%                 DOA_Out = [DOAout uniqs{j} '_' num2str(r) '.avi'];
%                 [cfar_bins] = RDC_to_rangeDopp(subRDC, RD_Out);
                RDC_to_microDopp(subRDC, mD_Out)
%                 RDC_to_rangeDOA_AWR1642(subRDC, DOA_Out)
                toc
        end

end
%end
