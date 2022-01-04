fp = 'C:\Users\3971627\OneDrive - Universiteit Utrecht\CFM image fo Acidity Map';

fn = {...
    '12-15-1-3.5kDa-mixture-5-4.nd2',...
    'FCC-Si-12-15-1-Allmirrors2.nd2',...
    };

lambda = 410:10:720;

ch1 = [27 28];
ch2 = [24 25];
skipPerc = 0;
printFlag = true;

p = printFig('printFlag',printFlag,'plotExt',{'-dpng','-dsvg'}); % 
nExpectedChannels = numel(lambda);


for ii = 1:numel(fn)
    
    [path,name,ext] = fileparts(fn{ii});
    p.savepath = fullfile(fp,path);
    I = bfopen(fullfile(fp,fn{ii}));
    
% for ii = 1:numel(allfiles)
%     [~,name,ext] = fileparts(allfiles(ii).name);
%     p.savepath = allfiles(ii).folder;
%     I = bfopen(fullfile(allfiles(ii).folder,allfiles(ii).name));    

    I = cat(3,I{1}{:,1});
    
    % file has been opened, now  save name without periods (. --> -)
    name = strrep(name,'.','-');
    
    nChannels = size(I,3);
    
    figure
    imagesc(mean(I,3))
    axis image
    colormap viridis
    
    if nChannels == nExpectedChannels
        p.print([name '_mean_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2))])
    else
        p.print([name '_mean'])
    end
    
    figure
    plot(squeeze(mean(mean(I,2),1)),'-o')
    xlabel('Spectral bin')
    ylabel('Mean intensity')
    
    if nChannels == nExpectedChannels
        p.print([name '_spectrum_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2))])
    else
        p.print([name '_spectrum'])
    end

    if nChannels == nExpectedChannels
        plotRGcolormap(I,ch1,ch2)
        p.print([name '_Rch1_Gch2_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2))])
        
        plotRatio(I,ch1,ch2,skipPerc)
        p.print([name '_ch1-div-ch2_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2)) '_skipp' strrep(num2str(skipPerc),'.','')])
        
        plotRatioLog10(I,ch1,ch2,skipPerc)
        p.print([name '_log10_ch1-div-ch2_ch1-' num2str(min(ch1)) '-' num2str(max(ch1)) '_ch2-' num2str(min(ch2)) '-' num2str(max(ch2)) '_skipp' strrep(num2str(skipPerc),'.','')])
    end
    
    close all
end

%%
function plotRGcolormap(I,ch1,ch2)

I = double(I);
Irgb = zeros(size(I,1),size(I,2),3);

% red channel
Irgb(:,:,1) = rescale(sum(I(:,:,ch1),3));

% green
Irgb(:,:,2) = rescale(sum(I(:,:,ch2),3));

imshow(Irgb)
axis image

end

function plotRatio(I,ch1,ch2,skipPerc)

I = double(I);

Iratio = sum(I(:,:,ch1),3)./sum(I(:,:,ch2),3);

figure
imagesc(Iratio,"AlphaData",rescale(sum(cat(3,I(:,:,ch1),I(:,:,ch2)),3)));
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

end

function plotRatioLog10(I,ch1,ch2,skipPerc)

I = double(I);

Iratio = log10(sum(I(:,:,ch1),3)./sum(I(:,:,ch2),3));

figure
imagesc(Iratio,"AlphaData",rescale(sum(cat(3,I(:,:,ch1),I(:,:,ch2)),3)));
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

end