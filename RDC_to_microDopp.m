function [] = RDC_to_microDopp( RDC, fOut)
        
        numADCBits = 16; % number of ADC bits per sample
        SweepTime = 40e-3; % Time for 1 frame
        NTS = size(RDC,1); %256 Number of time samples per sweep
        numADCSamples = NTS;
        numTX = 2; % '1' for 1 TX, '2' for BPM
        NoC = 128;%128; % Number of chirp loops
        NPpF = numTX*NoC; % Number of pulses per frame
        numRX = 4;
        
        numLanes = 2; % do not change. number of lanes is always 4 even if only 1 lane is used. unused lanes
        % NoF = fileSize/2/NPpF/numRX/NTS; % Number of frames
        numChirps = size(RDC,2);
        NoF = round(numChirps/NPpF); % Number of frames, 4 channels, I&Q channels (2)
        dT = SweepTime/NPpF; %
        prf = 1/dT; %
        RDC1=RDC(:,:,1);
        rp = fft(RDC1);
         figure; imagesc(20*log10(abs(rp)))
       %% MTI v2
            [b,a]=butter(1, 0.01, 'high'); %  4th order is 24dB/octave slope, 6dB/octave per order of n
        %                                      [B,A] = butter(N,Wn, 'high') where N filter order, b (numerator), a (denominator), ...
        %                                      highpass, Wn is cutoff freq (half the sample rate)
            [m,n]=size(rp);
            rngpro=zeros(m,n);
            for k=1:size(rp,1)
                rngpro(k,:)=filter(b,a,rp(k,:,1));
            end

        rBin=60:170;
 

        nfft = 2^12;window = 200;noverlap = 120;shift = window - noverlap;
        sx = myspecgramnew(sum(rngpro(rBin,:)),window,nfft,shift); % mti filter and IQ correction
      

        sx2 = abs(flipud(fftshift(sx,1)));
        %% Spectrogram
        timeAxis = [1:NPpF*NoF]*SweepTime/NPpF*numTX/2 ; % Time
        freqAxis = linspace(-prf/2,prf/2,nfft); % Frequency Axis
        
        figure('visible','on');
        C=colormap(jet(256));

      imagesc(timeAxis,((3e8)*[-prf/2 prf/2])/(2*77e9),20*log10(sx2./max(sx2(:))));
      ylim([-3.19 2.89])
        %set(gcf,'units','normalized','outerposition',[0,0,1,1]);
            axis xy
            set(gca,'FontSize',10)
            title(['RBin: ',num2str(rBin)]);
        title(fOut(end-21:end-6))
             xlabel('Time (sec)');
             ylabel('Velocity (M/S)');
             ylim([-3.19 2.89])
        caxis([-34 0]) % 40
        set(gca, 'YDir','normal')
        %     colorbar;
        
        
        
        Limit=((3e8)*(1500)/(2*77e9));
        axis([0 timeAxis(end) -Limit Limit])
      
        set(gca,'xtick',[],'ytick',[])
        %saveas(gcf,'myfig','png')
        
        
        %% trim Sx2
        
               figure('visible','on');
        colormap(jet(256));
        trimY_L=1000;
        trimY_U=3000;
        tot=trimY_U-trimY_L;
        sx4=flipud(sx2(trimY_L:trimY_U,:));
      imagesc(timeAxis,[-3.19 2.89],20*log10(sx4./max(sx4(:))));
     
        %set(gcf,'units','normalized','outerposition',[0,0,1,1]);
            axis xy
            set(gca,'FontSize',10)
            title(['RBin: ',num2str(rBin)]);
        title(fOut(end-21:end-6))
             xlabel('Time (sec)');
             ylabel('Velocity (M/S)');
        caxis([-34 0]) % 40
        set(gca, 'YDir','normal')
    
        img2=20*log10(sx4./max(sx4(:)));
        img2(img2<-34)=-100;
        figure;
        colormap(jet(256))
        imagesc(img2);
     
     %%   
        
     %G = 20*log10(sx4./max(sx4(:)));
        sx2 = abs(flipud(fftshift(sx,1)));
         sx4=flipud(sx2(trimY_L:trimY_U,:));
         img2=20*log10(sx4./max(sx4(:)));
        img2(img2<-34)=-100;
     G=img2;

% Now make an RGB image that matches display from IMAGESC:
C = colormap;  % Get the figure's colormap.
L = size(C,1);
% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(G(:)),max(G(:)),L),1:L,G));

MS=Gs;
MS(MS<=1)=0;
figure;
imagesc(MS)  
% Does this image match the other one?
title('IMAGE (MxNx3)')   

%%


[upper_env, central_env, lower_env] = env_find(MS);
figure;  plot(lower_env); hold on; plot(upper_env);hold on; plot(central_env)
a=upper_env;
b= central_env;
c=lower_env;
new_env=zeros(1,length(a));
for idx=1:length(central_env)
        if b(idx)<1000 
                if  a(idx)< 1000
                        
                        new_env(idx)=a(idx);
                else
                        new_env(idx)=b(idx);
                end
        end
        if b(idx)>1000
                if c(idx)> 1000
                        new_env(idx)=c(idx);
                else
                        new_env(idx)=b(idx);
                end
        end
       
end

  %figure;  plot(lower_env); hold on; plot(upper_env);hold on; plot(central_env);  hold on; plot(new_env,'r', 'LineWidth',2);               

 figure; imagesc(MS); hold on; plot(upper_env,'LineWidth',2);hold on; plot(lower_env,'LineWidth',2);hold on; plot(central_env,'LineWidth',2); hold on; plot(new_env,'r', 'LineWidth',1.3)

 y_vel=linspace(-6.2338,6.2338,4096);
 y_act=new_env+1000;
for kk=1:length(new_env) 
        y_new(kk)=y_vel(round(y_act(kk)));
end
leg_vel=-1*y_new;

% Inerpolate
x=1:length(new_env);
v=leg_vel;
xq=linspace(1,18000,18000);
ML_walk=interp1(x,v,xq);
 

figure; plot(leg_vel,'r'); hold on; plot(ML_walk,'--g');

 
        
%         frame = frame2im(getframe(gca));
%         imwrite(frame,[fOut(1:end-4) '.png']);
%         close all
        
end