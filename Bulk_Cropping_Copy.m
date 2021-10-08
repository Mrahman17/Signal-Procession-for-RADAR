clc; clear; warning off
saveColor = '/media/rspl-admin/Seagate Backup Plus Drive/20_words_Imitation Fall 2020/Spectrograms/Cropped_imit/';
%saveGray = 'C:\Users\mrahman17\Desktop\ASL Training Data Output\77 GHz\Front\microDoppler\Gray\Cropped\';
datapath = '/media/rspl-admin/Seagate Backup Plus Drive/20_words_Imitation Fall 2020/Spectrogram_oct/Emre 24 Oct/';

pattern = strcat(datapath,'*.fig');
files = dir(pattern);

% All Signs
Folder_name = {'You';'Hello';'Walk';'Drink';'Friend';'Knife';'Well';'Car';'Engineer';'Mountain';'Lawyer';'Hospital';'Health';'Earthquake';'Breath';'Help';'Push';'Go';'Come';'Write'};



I_MAX = numel(files); % # of files in "files" 
m = 100; % set 1 value below th desired, cropping dimensions - width - 168
n = 370; % set 1 value below th desired, cropping dimensions - height - 441
      
for ii =1:I_MAX 
    fname = strcat(datapath, files(2*ii-1).name);
    openfig(fname,'invisible');
    set(gcf, 'units','normalized','outerposition',[0 0 1.2 .65]);
    F = getframe(gca);
    [img,~] = frame2im(F);
    colormap(gray)
    Fgray = getframe(gca);
    [imgray,~] = frame2im(Fgray);
    
      OUT_DIR=strcat(saveColor,Folder_name{ii},'/');
         if 2~=exist(OUT_DIR,'dir')
                mkdir(OUT_DIR);
        end
    for jj=1:15
        msg = strcat('Crop (', int2str(ii),'/',int2str(I_MAX), ') |  ',...
            int2str(jj),'. image (',Folder_name(ii),')');
        disp(msg)
      
        name_color = strcat(OUT_DIR, Folder_name(ii),'_E_Oct_',num2str(jj), '.png');
        name_gray = strcat('/media/rspl-admin/Seagate Backup Plus Drive/20_ASL_Native/new_imit_128/',Folder_name(ii),'/', Folder_name(ii),'_E_Oct_',num2str(jj), '.png');
        imshow(img,'Border','tight')
%         set(gcf, 'units','normalized','outerposition',[0 .4 .5 .7]);
        h = imrect(gca,[m*jj 30 m n]);
        position = wait(h); %[leftto right, top to bottom, x length, y length]
        delete(h); 
        cropColor = imcrop(img, position);
        %cropGray = imcrop(imgray, position);
        imwrite(cropColor,name_color{1,1}) %change jj with the outer loop index while processing multiple images
        cropColor2=imresize(cropColor,[128 128]);
        imwrite(cropColor2,name_gray{1,1}) %change jj with the outer loop index while processing multiple images
        
    end
    close all
end