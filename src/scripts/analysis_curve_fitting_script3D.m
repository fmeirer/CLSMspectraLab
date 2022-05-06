fp = 'C:\Users\3971627\OneDrive - Universiteit Utrecht\CFM image fo Acidity Map';

fn = {...
    'F_Al_uZstained12151allmirrors.nd2'...
    };

% cd(fp)
% allfiles = dir('*/*.nd2');
% allfiles = dir('*.nd2');

lambda = 410:10:720;

ch1 = [27 28];
ch2 = [24 25];
nChannels = 32;
skipPerc = 0;
printFlag = true;

p = printFig('printFlag',printFlag,'plotExt',{'-dpng','-dsvg'}); % 
nExpectedChannels = numel(lambda);

for ii = 1:numel(fn) % numel(allfiles)
    
    [path,name,ext] = fileparts(fn{ii});
    p.savepath = fullfile(fp,path);
    
    
%     [~,name,ext] = fileparts(allfiles(ii).name);
%     p.savepath = allfiles(ii).folder;
    
    I = bfopen(fullfile(fp,fn{ii})); % bfopen(fullfile(allfiles(ii).folder,allfiles(ii).name));
    
    % file has been opened, now  save name without periods (. --> -)
    name = strrep(name,'.','-');
   
    I = cat(3,I{1}{:,1});
    nZ = size(I,3)/nChannels;
    I = reshape(I,[size(I,1) size(I,2) nChannels nZ]); % X Y C Z
    
    Imean = uint16(squeeze(mean(I,3)));
    
    for jj = 1:nZ
        pname = [name '_mean_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2)) '_3D.tiff'];
        if jj == 1
            imwrite(Imean(:,:,jj),fullfile(fp,pname))
        else
            imwrite(Imean(:,:,jj),fullfile(fp,pname),'WriteMode','append')
        end
    end
    
    figure
    plot(squeeze(mean(mean(mean(I,4),2),1)),'-o')
    xlabel('Spectral bin')
    ylabel('Mean intensity')

   
    p.print([name '_spectrum_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2)) '_3D'])


    pname = [name '_Rch1_Gch2_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2)) '_3D.tiff'];
    pname_auto = [name '_Rch1_Gch2_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2)) '_3D_autoscale.tiff'];
    plotRGcolormap3D(I,ch1,ch2,fullfile(fp,pname_auto),fullfile(fp,pname))
    
    pname = [name '_ch1-div-ch2_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2)) '_skipp' strrep(num2str(skipPerc),'.','') '_3D.tiff'];
    plotRatio3D(I,ch1,ch2,skipPerc,p,fullfile(fp,pname))
    
    pname = [name '_log10_ch1-div-ch2_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2)) '_skipp' strrep(num2str(skipPerc),'.','') '_3D.tiff'];
    plotRatio3DLog10(I,ch1,ch2,skipPerc,p,fullfile(fp,pname))
    close all
end

function plotRGcolormap3D(I,ch1,ch2,fname_auto,fname)

I = double(I);

% autoscale
for ii = 1:size(I,4)
    Irgb = zeros(size(I,1),size(I,2),3);
    
    % red channel
    Irgb(:,:,1) = rescale(sum(squeeze(I(:,:,ch1,ii)),3));
    
    % green
    Irgb(:,:,2) = rescale(sum(squeeze(I(:,:,ch2,ii)),3));
    
    if ii == 1
        imwrite(Irgb,fname_auto)
    else
        imwrite(Irgb,fname_auto,'WriteMode','append')
    end
    
end

% non-autoscaled
Irgb = zeros(size(I,1),size(I,2),3,size(I,4));
for ii = 1:size(I,4)
    
    % red channel
    Irgb(:,:,1,ii) = sum(squeeze(I(:,:,ch1,ii)),3);
    
    % green
    Irgb(:,:,2,ii) = sum(squeeze(I(:,:,ch2,ii)),3);
    
end

Irgb(:,:,1,:) = rescale(Irgb(:,:,1,:)); %
Irgb(:,:,2,:) = rescale(Irgb(:,:,2,:));
    
for ii = 1:size(I,4)
    
    if ii == 1
        imwrite(Irgb(:,:,:,ii),fname)
    else
        imwrite(Irgb(:,:,:,ii),fname,'WriteMode','append')
    end
    
end

end

function plotRatio3D(I,ch1,ch2,skipPerc,p,fname)

I = double(I);

[path,name,~] = fileparts(fname);

path = fullfile(path,['ratio_' name]);
if ~exist(path,'dir') % is folder
    mkdir(path) % make folder
end
p.savepath = path;

for ii = 1:size(I,4)
    
    Iratio = sum(I(:,:,ch1,ii),3)./sum(I(:,:,ch2,ii),3);
    
    figure
    imagesc(Iratio,"AlphaData",rescale(sum(cat(3,I(:,:,ch1,ii),I(:,:,ch2,ii)),3)));
    axis image
    
    % skip brightest pixels
    if skipPerc ~= 0
        s = sort(Iratio(:));
        upperlim = s(round(numel(s)*(100-skipPerc)/100));
        if ~isnan(upperlim)
            l = caxis;
            caxis([l(1) upperlim])
        end
    end
    colormap jet
    colorbar
    
    p.print([name '_' num2str(ii)])
    
end

end

function plotRatio3DLog10(I,ch1,ch2,skipPerc,p,fname)

I = double(I);

[path,name,~] = fileparts(fname);

path = fullfile(path,['ratio_' name]);
if ~exist(path,'dir') % is folder
    mkdir(path) % make folder
end
p.savepath = path;

for ii = 1:size(I,4)
    
    Iratio = log10(sum(I(:,:,ch1,ii),3)./sum(I(:,:,ch2,ii),3));
    
    figure
    imagesc(Iratio,"AlphaData",rescale(sum(cat(3,I(:,:,ch1,ii),I(:,:,ch2,ii)),3)));
    axis image
    
    % skip brightest pixels
    if skipPerc ~= 0
        s = sort(Iratio(:));
        upperlim = s(round(numel(s)*(100-skipPerc)/100));
        if ~isnan(upperlim)
            l = caxis;
            caxis([l(1) upperlim])
        end
    end
    colormap jet
    colorbar
    
    p.print([name '_' num2str(ii)])
    
end

end