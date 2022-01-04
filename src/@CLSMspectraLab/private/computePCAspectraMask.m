function [obj, varargout] = computePCAspectraMask(obj,plotFlag)
% GETPCASPECTRAMASK computes the domains with similar spectra and stores a
% mask
%
% The first principal component of the image/volume is computed and
% thresholded using Otsu's method. The method can be extented to
% multithesholding, but this has not been implemented yet.
%
% Usage:
% obj = getPCAspectraMask(obj) computes and stores the mask in obj.Imask.
%
% [obj, threshold] = getPCAspectraMask(obj) also gives the threshold of the
% first eigenimage, indexing: [z, number of thresholds].

% hardcoded variables
n_ths = 1; % number of thresholds (only one now, check commented sections for multithreshold behavior)
medfiltersize = [7 7 7];

if nargin < 2
    plotFlag = false;
end

% First delete old mask
obj = obj.deleteMask;

% extract full data volume: X x Y x spectrum x z-index
Vr = obj.getReshapedImage();

m = obj.imageX;
n = obj.imageY;
z = obj.stackSize;
chn = obj.nChannels;

% PCA
wbh=waitbar(0,'Doing PCA... please wait');

I = nan(n,m,z); % preallocation
for ii=1:z
    % get z slice: x y ch
    Vxyc = reshape(squeeze(Vr(:,:,:,ii)),n*m,chn);
    % do a PCA:
    [~, score, ~] = pca(Vxyc,'Centered','on','Economy','on');
    % get first eigenimage
    I(:,:,ii) = reshape(score(:,1),n,m);
    % get thresholds of eigenimage using Otsu's method
    th(ii,:) = multithresh(I(:,:,ii),n_ths);
    % update wait bar
    waitbar(ii/z,wbh);
end
close(wbh);

if plotFlag
    figure();
    plot(th);
    title('threshold(s) for each slice');
    xlabel('Slice #')
    ylabel('Threshold')
    
    [nR,nC] = subplotDimensions(z);
    figure
    for ii = 1:z
        thisI = I(:,:,ii);
        subplot(nR,nC,ii);
        histogram(thisI(:))
        hold on
        for jj = 1:n_ths
            xline(th(ii,jj));
        end
        title('Eigenimage threshold')
        xlabel('Score')
        ylabel('Counts')
    end
end

if obj.isVolume
    % use average of last 10 values (check plot if OK)
    if ii>9
        for ii=1:n_ths
            th_f(1,ii) = mean(th(end-10:end,ii));
        end
    else
        for ii=1:n_ths
            th_f(1,ii) = mean(th(:,ii));
        end
    end
    % and binarize the whole stack:
    Ibin = imbinarize(I,th_f(1,1));

    % remove salt and pepper noise:
    Ibin_mdf = medfilt3(Ibin,medfiltersize);
    % fill holes:
    Ibin_mdf_filled = imfill(Ibin_mdf,18,'holes'); 
else
    Ibin_mdf_filled = nan(n,m,z);
    for ii = 1:z
        % and binarize the whole stack:
        Ibin = imbinarize(I(:,:,ii),th(ii,1));

        % remove salt and pepper noise:
        Ibin_mdf = medfilt2(Ibin,medfiltersize(1:2));
        % fill holes:
        Ibin_mdf_filled(:,:,ii) = imfill(Ibin_mdf,8,'holes');
    end
end

if plotFlag
    if obj.isVolume
        figure
        channelInt = nan(obj.nChannels,2);
        for ii = 1:obj.nChannels
            Ichannel = squeeze(Vr(:,:,ii,:));
            % background
            channelInt(ii,1) = nanmean(Ichannel(Ibin_mdf_filled == 0)); % empty gives NaN
            % foreground
            channelInt(ii,2) = nanmean(Ichannel(Ibin_mdf_filled > 0)); % empty gives NaN
        end
        plot((1:obj.nChannels),channelInt(:,1),'.-');
        hold on
        plot((1:obj.nChannels),channelInt(:,2),'.-');
        xlabel('Channel')
        ylabel('Mean intensity')
        legend({'Background','Foreground'})
    else
        figure
        for jj = 1:z
            channelInt = nan(obj.nChannels,2);
            for ii = 1:obj.nChannels
                Ichannel = squeeze(Vr(:,:,ii,jj));
                % background
                channelInt(ii,1) = nanmean(Ichannel(Ibin_mdf_filled(:,:,jj) == 0)); % empty gives NaN
                % foreground
                channelInt(ii,2) = nanmean(Ichannel(Ibin_mdf_filled(:,:,jj) > 0)); % empty gives NaN
            end
            subplot(nR,nC,jj);
            plot((1:obj.nChannels),channelInt(:,1),'.-');
            hold on
            plot((1:obj.nChannels),channelInt(:,2),'.-');
            xlabel('Channel')
            ylabel('Mean intensity')
            legend({'Background','Foreground'})
        end

    end
end

% store in object
obj.Imask = Ibin_mdf_filled;

if nargout > 1
    varargout{1} = th;
end

% % for plotting
% [nR,nC] = subplotDimensions(z);
% 
% 
% % fix for n_ths > 1
% % if n_ths>1
% %     for ii=2:n_ths
% %         Ibin = Ibin + imbinarize(I,th_f(1,ii)).*ii; % get mask with 0 for I < th_f(1); 1 for th_f(1) < I < th_f(2); etc.. 3, 6, 10
% %         % why *ii ?? to get 0,1,3,6,10 scaling
% %     end
% % end
% 
% 
% % if UserMask
% %     % let user select region of interest:
% %     Ibin_sum = sum(Ibin,3);
% %     figure('Name', 'please select the ROI'); imagesc(Ibin_sum); axis image;
% %     h = imfreehand('Closed','true');
% %     UserMask = repmat(h.createMask,1,1,z);
% %     Ibin = Ibin.*UserMask;
% % end
% 
% % WARNING: this piece of code is weird for n_ths > 1, because the medfilt3
% % function applies to the binarized image with irregular scaling per
% % threshold, i.e. 0,1,3,6,10,... Let's for now stick to one threshold only.
% 
% % remove salt and pepper noise:
% Ibin_mdf = medfilt3(Ibin,[7 7 7]);
% % fill holes:
% Ibin_mdf_filled = imfill(Ibin_mdf,18,'holes'); 
% 
% % for debugging
% if plot_flag
%     if z<10
%         figure('name','TPV');
%         for ii=1:z
%            subplot(nR,nC,ii); 
%            imagesc(Ibin_mdf_filled(:,:,ii)); 
%            axis image;
%         end
%     else
%         volumeViewer(Ibin_mdf_filled);
%         input('Displaying TPV. Ready to continue? Press Enter... (it is a good idea to first close the volume viewer)','s');
%     end
% end
% 
% % bin data:
% I = imresize(I,binfact); % only resizes first two dimensions
% Ibin_mdf_filled = imresize(Ibin_mdf_filled,binfact);
% [n,m,z]=size(I);
% 
% % for debugging:
% if plot_flag
%     if z<10
%         figure('name','binned data');
%         for ii=1:z
%            subplot(nR,nC,ii);
%            imagesc(I(:,:,ii)); 
%            axis image;
%         end
%     else
%         volumeViewer(I);
%         input('Displaying binned data. Ready to continue? Press Enter... (it is a good idea to first close the volume viewer)','s');
%     end
% end
% 
% % mask volume:
% Imasked = I.*Ibin_mdf_filled; % weird behaviour if multiple thresholds are selected. The magnitude of the intensity is increased by the cluster it belongs to, because Ibin =/= [0 1]
% th = zeros(z,2);
% if plot_flag
%     if z<10
%         figure('name','segmented sample (TPV)');
%         for ii=1:z
%            subplot(nR,nC,ii);
%            imagesc(Imasked(:,:,ii)); 
%            axis image;
%            % within masked region do again thresholding:
%            th(ii,:) = multithresh(Imasked(:,:,ii),2);
%         end
%     else
%         volumeViewer(Imasked);
%         input('Displaying segmented sample (TPV). Ready to continue? Press Enter... (it is a good idea to first close the volume viewer)','s');
%     end
% 
%     figure();
%     plot(th);
%     title('threshold for each slice (within TPV)');
% end
% % use average of last 10 values (check plot if OK)
% %th_f = mean(th(end-10:end,1));
% % use maximum of all th's:
% th_f = max(th(:,1));
% % and binarize the whole stack:
% Imasked_bin = imbinarize(Imasked,th_f);
% % figure();
% % for ii=1:z
% %    subplot(nR,nC,ii);
% %    imagesc(Imasked_bin(:,:,ii)); 
% %    axis image;
% % end
% % use average of last 10 values (check plot if OK)
% %th_f = mean(th(end-10:end,2));
% % use maximum of all th's:
% th_f = max(th(:,2));
% % and binarize the whole stack:
% Imasked_bin = imbinarize(Imasked,th_f);
% 
% % remove salt and pepper noise:
% If2 = medfilt3(Imasked_bin,[3 3 3]);
% % fill holes:
% %If2 = imfill(If1,18,'holes'); 
% 
% if plot_flag
%     if z<10
%         figure('name','segmented regions of highest intensity');
%         for ii=1:z
%            subplot(nR,nC,ii); 
%            imagesc(If2(:,:,ii)); 
%            axis image;
%         end
%     else
%         volumeViewer(If2);
%         input('Displaying segmented regions of highest intensity. Ready to continue? Press Enter... (it is a good idea to first close the volume viewer)','s');
%     end
% end
% 
% % put this into one result matrix:
% R = If2 + Ibin_mdf_filled;
% % figure();
% % vol3d('cdata',R,'texture','3D');
% % view(3); 
% % grid on;
% 
% % % % save for Avizo:
% % % save(['LabeledVol_',sfname,'.mat'],'R');
% % % Rbin = R;
% % % Rbin(R>0)=1;
% % % save(['BinaryVol_',sfname,'.mat'],'Rbin');
% 
% 
% % get the spectrum of the 3 regions from binned data:
% 
% Vr_binn = imresize(Vr,binfact); % only resizes first two dimensions
% Vr_binn_Chnsum = squeeze(sum(Vr_binn,3));
% Vr_segmented = Vr_binn_Chnsum.*Ibin_mdf_filled;
% % % save(['SegmentedVol_',sfname],'Vr_segmented');
% 
% for ii=1:chn
%     % for the same channel get pixels of the 3 regions:
%     Vchn = squeeze(Vr_binn(:,:,ii,:)); % this is a volume of size R for one channel
%     Vchn1(:,ii) = Vchn(R==0); % all voxels outside sample
%     Vchn2(:,ii) = Vchn(R==1); % all voxels outside sample
%     Vchn3(:,ii) = Vchn(R==2); % all voxels outside sample
% end
% 
% np1 = size(Vchn1,1);
% np2 = size(Vchn2,1);
% np3 = size(Vchn3,1);
% 
% output_text{1} = sprintf('Nr. of pixel of material 1: %i vol. fraction: %.2f',np1,np1/(np1+np2+np3));
% output_text{2} = sprintf('Nr. of pixel of material 2: %i vol. fraction: %.2f',np2,np2/(np1+np2+np3));
% output_text{3} = sprintf('Nr. of pixel of material 3: %i vol. fraction: %.2f',np3,np3/(np1+np2+np3));
% 
% if plot_flag
%     figure(); 
%     title('Spectra of regions');
%     plot(mean(Vchn1,1)); hold on;
%     plot(mean(Vchn2,1)); hold on;
%     plot(mean(Vchn3,1));
%     xlabel ('channel');
%     ylabel('intensity');
%     legend({'Average spectrum of material 1','Average spectrum of material 2','Average spectrum of material 3'});
% 
%     figure(); 
%     title('Spectra of regions, normalized');
%     plot(mean(Vchn1,1)./(mean(mean(Vchn1,1),2))); hold on;
%     plot(mean(Vchn2,1)./(mean(mean(Vchn2,1),2))); hold on;
%     plot(mean(Vchn3,1)./(mean(mean(Vchn3,1),2)));
%     xlabel ('channel');
%     ylabel('intensity');
%     legend({'Average spectrum of material 1 (divided by average intensity)',...
%         'Average spectrum of material 2 (divided by average intensity)',...
%         'Average spectrum of material 3 (divided by average intensity)'});
% end