%reading in the "caked data"

clear all
close all

root = '/Volumes/FHofmann Oxf/I22 Tendons/processing/output';
folders = ['i22-207088_processed';
           'i22-207091_processed';
           'i22-207092_processed';
           'i22-207093_processed';
           'i22-207094_processed'];
fileroot = 'data';

filenums = [0: 20]';

peak_pos_top = [55 118     245 309 372     499 563 626     752]'; % top peak positions roughly
peak_pos_bot = [   118 181 245 308 372 435     562 626 689]';     %bottom peak positions

q_cent_top = zeros(size(filenums,1), size(folders,1), size(peak_pos_top,1));
intensity_top = zeros(size(filenums,1), size(folders,1), size(peak_pos_top,1));
fwhm_top = zeros(size(filenums,1), size(folders,1), size(peak_pos_top,1));
q_cent_bot = zeros(size(filenums,1), size(folders,1), size(peak_pos_top,1));
intensity_bot = zeros(size(filenums,1), size(folders,1), size(peak_pos_top,1));
fwhm_bot = zeros(size(filenums,1), size(folders,1), size(peak_pos_top,1));
q_scale = 0.001430552:0.00014652308578388:0.159528961676387;

% background1 = double(imread([root,'/' ,folders(1,:),'/', fileroot, num2str(filenums(1), '%03u'), '.tif']));
% background2 = double(imread([root,'/' ,folders(2,:),'/', fileroot, num2str(filenums(1), '%03u'), '.tif']));
% background = (background1 + background2)./2;

%compute background...

for ii = 1:size(folders,1);
    ii
    for jj = 1;
        foldername = dir([root,'/' ,folders(ii,:), '*']);
        filename = dir([root,'/' ,foldername.name,'/', '*', fileroot, num2str(filenums(jj), '%03u'), '*.tif']);
        filepath = [root,'/' ,foldername.name,'/',filename.name]
        background_st(:,:,ii) = double(imread(filepath));
    end
end

background = mean(background_st,3);


for ii = 1 : size(folders,1);
    ii
    for jj = 1:size(filenums,1);
        foldername = dir([root,'/' ,folders(ii,:), '*']);
        filename = dir([root,'/' ,foldername.name,'/', '*', fileroot, num2str(filenums(jj), '%03u'), '*.tif']);
        filepath = [root,'/' ,foldername.name,'/',filename.name]
        I = double(imread(filepath))-background;
        I(I<0)=0;
        I_line = mean(I(253:293,:),1);
%         plot(I_line);
%         pause;
        %now pick out one peak and fit... 
        for kk = 1: size(peak_pos_top,1);
            I_dat = I_line(peak_pos_top(kk)-20:peak_pos_top(kk)+20);
            q_pos = q_scale(peak_pos_top(kk)-20:peak_pos_top(kk)+20)';
            w= 0.4 + sqrt(I_dat./max(I_dat));
            [Alt,k0]= max(I_dat);
            x0(3) = Alt-min(I_dat);
            x0(1)=q_pos(k0);
            x0(2)=0.0001;
            x0(4)=min(I_dat);
            
            [yfit,x,~,~,~,~]=leasqr(q_pos,I_dat,x0,'gauss',0.0000001,40, w);
            
            %Storing the results
            if (x(1)> min(q_pos)) && (x(1) < max(q_pos))
                q_cent_top(jj,ii,kk)=real(x(1)); %1st index - position, 2nd index slice, 3rd index reflection order
                intensity_top(jj,ii,kk)=real(x(3)); %1st index - position, 2nd index slice, 3rd index reflection order
                fwhm_top(jj,ii,kk) = x(2); %1st index - position, 2nd index slice, 3rd index reflection order
            else
                q_cent_top(jj,ii,kk)=0; %1st index - position, 2nd index slice, 3rd index reflection order
                intensity_top(jj,ii,kk)=0; %1st index - position, 2nd index slice, 3rd index reflection order
                fwhm_top(jj,ii,kk) = 0;%1st index - position, 2nd index slice, 3rd index reflection order
            end
            
%             figure(1)
%             plot(q_pos, I_dat, 'go'); hold on
%             plot(q_pos, yfit, 'b-'); hold off
%             title(['***TOP*** image no:', num2str(jj), ' peak no:', num2str(kk)])
%             pause;
        end
        I_line = mean(I(794:834,:),1);
        for kk = 1: size(peak_pos_bot,1);
            I_dat = I_line(peak_pos_bot(kk)-20:peak_pos_bot(kk)+20);
            q_pos = q_scale(peak_pos_bot(kk)-20:peak_pos_bot(kk)+20)';
            w= 0.4 + sqrt(I_dat./max(I_dat));
            [Alt,k0]= max(I_dat);
            x0(3) = Alt-min(I_dat);
            x0(1)=q_pos(k0);
            x0(2)=0.0001;
            x0(4)=min(I_dat);
            
            [yfit,x,~,~,~,~]=leasqr(q_pos,I_dat,x0,'gauss',0.0000001,40, w);
            
            %Storing the results
            if (x(1)> min(q_pos)) && (x(1) < max(q_pos))
                q_cent_bot(jj,ii,kk)=real(x(1)); %1st index - position, 2nd index slice, 3rd index reflection order
                intensity_bot(jj,ii,kk)=real(x(3)); %1st index - position, 2nd index slice, 3rd index reflection order
                fwhm_bot(jj,ii,kk) = x(2);%1st index - position, 2nd index slice, 3rd index reflection order
            else
                q_cent_bot(jj,ii,kk)=0; %1st index - position, 2nd index slice, 3rd index reflection order
                intensity_bot(jj,ii,kk)=0; %1st index - position, 2nd index slice, 3rd index reflection order
                fwhm_bot(jj,ii,kk) = 0; %1st index - position, 2nd index slice, 3rd index reflection order
            end
            
%             figure(1)
%             plot(q_pos, I_dat, 'go'); hold on
%             plot(q_pos, yfit, 'b-'); hold off
%             title(['***BOT*** image no:', num2str(jj), ' peak no:', num2str(kk)])
%             pause;
        end        
    end
end

cmap = [0:0.01:1; zeros(1,101); 1:-0.01:0;]';

intensity = squeeze(sum(intensity_top,3)+sum(intensity_bot,3));
intensity = intensity./max(max(intensity));

% figure; 
% colormap(cmap)
% imagesc(intensity, 'AlphaData', intensity./max(max(intensity))); axis equal
% 
% figure; plot(q_cent_top(:,1,1),'r'); hold on; plot(q_cent_top(:,2,1),'g');
% figure; plot(intensity(:,1),'r'); hold on; plot(intensity(:,2),'g');
% 
% q0 = sum(q_cent(:,1).*intensity(:,1))/sum(intensity(:,1));
% 
% strain = q0./q_cent - 1;

alphamap = intensity./max(max(intensity));

alphamap = intensity;
alphamap(alphamap <0.1)=0; alphamap(alphamap >0.09)=1;
alpha_clean = bwareaopen(alphamap, 3);
% %imagesc(strain, 'AlphaData', alpha_clean); axis equal
% colormap(cmap)
% caxis([-0.002 0.015])
% colorbar
% title('tendon 62')

intensity_clean_nn = intensity.* alpha_clean;
intensity_clean = intensity_clean_nn./max(max(intensity_clean_nn)); 
for kk = 1: size(peak_pos_top,1);
intensity_clean_3D_top(:,:,kk) = intensity_clean;
end
for kk = 1: size(peak_pos_bot,1);
intensity_clean_3D_bot(:,:,kk) = intensity_clean;
end

subplot(1,2,1); imagesc(intensity);
subplot(1,2,2); imagesc(intensity_clean);

format short e
%calculate itensity weighted mean q...
q_mean_top = squeeze((sum(q_cent_top.*intensity_clean_3D_top,1)./sum(intensity_clean_3D_top,1)));
q_mean_bot = squeeze((sum(q_cent_bot.*intensity_clean_3D_bot,1)./sum(intensity_clean_3D_bot,1)));

%calculate q normalised by order....
order_mat_top = (ones(size(folders,1),1)*[1 2 4 5 6 8 9 10 12]);
q_mean_top_norm = q_mean_top./order_mat_top

order_mat_bot = (ones(size(folders,1),1)*[2 3 4 5 6 7 9 10 11]);
q_mean_bot_norm = q_mean_bot./order_mat_bot

%calculate intensity weighted fwhm...
fwhm_mean_top = squeeze((sum(fwhm_top.*intensity_clean_3D_top,1)./sum(intensity_clean_3D_top,1)))
fwhm_mean_bot = squeeze((sum(fwhm_bot.*intensity_clean_3D_bot,1)./sum(intensity_clean_3D_bot,1)))


%calculate summed intensity
intensity_sum_top = squeeze((sum(intensity_top,1)))
intensity_sum_bot = squeeze((sum(intensity_bot,1)))

%out = [q_mean_top_norm, q_mean_bot_norm, fwhm_mean_top, fwhm_mean_bot, intensity_sum_top, intensity_sum_bot]


