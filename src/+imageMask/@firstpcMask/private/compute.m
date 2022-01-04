function obj = compute(obj,imageStack)
%COMPUTE Compute the mask
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
medfiltersize = obj.medianFilterSize;

if nargin < 2
    plotFlag = false;
end

% First delete old mask
obj = obj.deleteMask;

% extract full data volume: X x Y x spectrum x z-index
Vr = imageStack.getReshapedImage('x','y','c','z');

m = obj.imageX;
n = obj.imageY;
z = obj.stackSize;
chn = obj.nChannels;

% PCA
wbh=waitbar(0,'Doing PCA... please wait');

I = nan(n,m,z); % preallocation
for ii = 1:z
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

end
