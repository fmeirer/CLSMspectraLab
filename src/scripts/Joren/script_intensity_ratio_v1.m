%% Intensity Ratio Script for Joren
% Erik Maris 2022-10-27
% This script plots the normalised intensity ratio of three channels as
% histogram and heatmap. The mean and std of the intensity ratio 
% distribution is saved as .csv file.

%% Variables

% define channels
ch1 = 13:18;
ch2 = 19:24;
ch3 = 25:30;

% size of the bin in pixels. 2 gives a 2x2 bin.
bin_factor = 2;

% first bin wavelength (nm)
start_wl = 450;
% last bin wavelength (nm)
end_wl = 750;
% number of bins
nBins = 30;

% threshold of the mask (use GUI to find appropriate value); use [] for
% automatic thresholding
mask_threshold = [0,0,0];

% save figures? true/false
printFlag = true;

% normalization of histograms. Options: 
% https://www.mathworks.com/help/releases/R2021b/matlab/ref/matlab.graphics.chart.primitive.histogram-properties.html
% and check Normalization options
normalization = 'probability';

% filepath to data (relative to folder script)
fp = '..\CFM data\2023-08-23-Withoilsession';
% fp = 'C:\Users\3971627\OneDrive - Universiteit Utrecht\Joren\CFM data\2022-10-11';
fn = {'2023-08-23-CNT-22-2383-Pu8-COTPP-oil.nd2',...
      '2023-08-23-CNT-22-2434-Pu7-COTPP-oil.nd2',...
      '2023-08-23-CNT-22-2797-Pu7-COTPP-oil.nd2'};

% file path to dark (relative to folder script)
fp_dark = '..\CFM data\2022-10-11 - Background Measurement 2';
% fp_dark = 'C:\Users\3971627\OneDrive - Universiteit Utrecht\Joren\CFM data\2022-10-11';
fn_dark = '2022-10-11-Background-1.nd2';

% filepath to save folder (relative to folder script)
fp_save = 'results from Erics software\Extra CNT';

% file path to CLSM spectra lab code
fp_code = 'MATLAB CLSM Spectra Lab';
% file path to printFig code
fp_pf = 'MATLAB printFig';

%% add code to path
[fp_current,~,~] = fileparts(mfilename("fullpath"));
cd(fp_current)
addpath(genpath(fp_code))
addpath(genpath(fp_pf))

%% set up printFig
p = printFig('savepath',fp_save,'plotExt',{'-dpng','fig'},'printFlag',printFlag);
savenames = cell(numel(fn),1);
for ii = 1:numel(savenames)
    [~,n,~] = fileparts(fn{ii});
    savenames{ii} = n;
end

%% load data and bin
cl = CLSMspectraLab(fullfile(fp,fn));
cl_bin = cl.getBinnedObj(bin_factor);

%% perform background correction
cl_bin.bgCorrection = bgCorrection.reference;
Iref = imageStack.import(fullfile(fp_dark,fn_dark)); % import background image
Iref_bin = Iref.binImageFirst2Dims(bin_factor,{'x','y','z','t','c'});
cl_bin.bgCorrection(1).Iref = deal(Iref_bin);
cl_bin.bgCorrection(1).uniformBackgroundFlag = false;
cl_bin = cl_bin.computeBgCorrection;

%% apply mask
if ~isempty(mask_threshold)
    if cl_bin.nInput ~= numel(mask_threshold)
        error('The number of thresholds does not match the number of input files')
    end
    for ii = 1:cl_bin.nInput
        cl_bin.mask(ii).thresholds = mask_threshold(ii);
    end
end
cl_bin = cl_bin.computeMask;

%% Plot results
% preallocation
ratio_mean = nan(cl.nInput,5);
ratio_std = nan(cl.nInput,5);

rgb_bins = getRGBbins(start_wl,end_wl,nBins);

for ii = 1:cl.nInput
    % plot masked image
    figure; imagesc(cl_bin.getImageProcessed(ii,'mean','x','y')); colorbar
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_meanIntensity'])

    % plot background-corrected spectrum
    figure; subplot(1,2,1); cl_bin.bgCorrectionFlag = false; cl_bin.plotChannels([],ii); cl_bin.bgCorrectionFlag = true; 
    title('Raw spectrum'); xlabel('Bin'); ylabel('Counts')
    subplot(1,2,2); cl_bin.plotChannels([],ii)
    title('Background-corrected spectrum'); xlabel('Bin'); ylabel('Counts')
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_bgCorrectionSpectrum'])
 
    % plot channel ratio histogram
    [Ich1,Ich2,Ich3] = getIbinsProcessed(cl_bin,ii,'sum','c',{ch1,ch2,ch3},'x','y','c');
    Nf = sum(Ich1,3)+sum(Ich2,3)+sum(Ich3,3); % normalization factor

    Rch1 = sum(Ich1,3)./Nf; Rch1(isnan(Rch1)) = []; % div by 0/0=nan (masked out); removing nan creates vector
    figure; histogram(Rch1,'Normalization',normalization); xlim([0 1]); caxis([0 1])
    xlabel('ch1/(ch1+ch2+ch3)'); ylabel(normalization)
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_histogramPeak1_123'])
    ratio_mean(ii,1) = mean(Rch1(:));
    ratio_std(ii,1) = std(Rch1(:));

    Rch2 = sum(Ich2,3)./Nf; Rch2(isnan(Rch2)) = [];
    figure; histogram(Rch2,'Normalization',normalization); xlim([0 1]); caxis([0 1])
    xlabel('ch2/(ch1+ch2+ch3)'); ylabel(normalization)
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_histogramPeak2_123'])
    ratio_mean(ii,2) = mean(Rch2(:));
    ratio_std(ii,2) = std(Rch2(:));

    Rch3 = sum(Ich3,3)./Nf; Rch3(isnan(Rch3)) = [];
    figure; histogram(Rch3,'Normalization',normalization); xlim([0 1]); caxis([0 1])
    xlabel('ch3/(ch1+ch2+ch3)'); ylabel(normalization)
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_histogramPeak3_123'])
    ratio_mean(ii,3) = mean(Rch3(:));
    ratio_std(ii,3) = std(Rch3(:));

    Rch_13 = sum(Ich1,3)./(sum(Ich1,3)+sum(Ich3,3)); Rch_13(isnan(Rch_13)) = [];
    figure; histogram(Rch_13,'Normalization',normalization); xlim([0 1]); caxis([0 1])
    xlabel('ch1/(ch1+ch3)'); ylabel(normalization)
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_histogramPeak1_13'])
    ratio_mean(ii,4) = mean(Rch_13(:));
    ratio_std(ii,4) = std(Rch_13(:));

    Rch_tot = Nf(:); Rch_tot(Rch_tot==0) = []; % is masked, = 0; removing nan creates vector
    figure; histogram(Rch_tot,'Normalization',normalization)
    xlabel('ch1+ch2+ch3'); ylabel(normalization)
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_histogramPeak123'])
    ratio_mean(ii,5) = mean(Rch_tot(:));
    ratio_std(ii,5) = std(Rch_tot(:));

    % plot ratio map
    figure
    cl_bin.plotRatio_2D([],ii,1,1,ch1,[ch1 ch2 ch3],false); caxis([0 1])
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_map1_123'])
    figure
    cl_bin.plotRatio_2D([],ii,1,1,ch1,[ch1 ch2 ch3],true); caxis([0 1])
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_map1_123_norm'])

    figure
    cl_bin.plotRatio_2D([],ii,1,1,ch2,[ch1 ch2 ch3],false); caxis([0 1])
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_map2_123'])
    figure
    cl_bin.plotRatio_2D([],ii,1,1,ch2,[ch1 ch2 ch3],true); caxis([0 1])
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_map2_123_norm'])

    figure
    cl_bin.plotRatio_2D([],ii,1,1,ch3,[ch1 ch2 ch3],false); caxis([0 1])
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_map3_123'])
    figure
    cl_bin.plotRatio_2D([],ii,1,1,ch3,[ch1 ch2 ch3],true); caxis([0 1])
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_map3_123_norm'])

    figure
    cl_bin.plotRatio_2D([],ii,1,1,ch1,[ch1 ch3],false); caxis([0 1])
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_map1_13'])
    figure
    cl_bin.plotRatio_2D([],ii,1,1,ch1,[ch1 ch3],true); caxis([0 1])
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_map1_13_norm'])
    
    % plot true color image
    figure
    cl_bin.plotTrueColor2D([],ii,1,1,rgb_bins)
    p.print([savenames{ii} '_bf' num2str(bin_factor) '_trueColor'])
    
    % clear all plots
    close all
end

%% save csv with mean and std of distribution
if printFlag
    T = cell2table(num2cell(ratio_mean),'RowNames',savenames,'VariableNames',{'1/(1+2+3)','2/(1+2+3)','3/(1+2+3)','1/(1+3)','1+2+3'});
    writetable(T,fullfile(fp_save,'ratio_mean.csv'),'WriteRowNames',true)
    
    T = cell2table(num2cell(ratio_std),'RowNames',savenames,'VariableNames',{'1/(1+2+3)','2/(1+2+3)','3/(1+2+3)','1/(1+3)','1+2+3'});
    writetable(T,fullfile(fp_save,'ratio_std.csv'),'WriteRowNames',true)
end
%% plot distributions overlay

% go to View --> Plot Browser in the Figure panel to select the data you
% wish to display
figure
for ii = 1:cl.nInput
    [Ich1,Ich2,Ich3] = getIbinsProcessed(cl_bin,ii,'sum','c',{ch1,ch2,ch3},'x','y','c');
    Nf = sum(Ich1,3)+sum(Ich2,3)+sum(Ich3,3); % normalization factor
    Rch1 = sum(Ich1,3)./Nf; Rch1(isnan(Rch1)) = [];
    Rch2 = sum(Ich2,3)./Nf; Rch2(isnan(Rch2)) = [];
    Rch3 = sum(Ich3,3)./Nf; Rch3(isnan(Rch3)) = [];
    
    if ii == 1
        subplot(1,3,1); histogram(Rch1,'Normalization',normalization,'DisplayName',savenames{1}); hold on
        [~,edges1] = histcounts(Rch1);
        subplot(1,3,2); histogram(Rch2,'Normalization',normalization,'DisplayName',savenames{1}); hold on
        [~,edges2] = histcounts(Rch2);
        subplot(1,3,3); histogram(Rch3,'Normalization',normalization,'DisplayName',savenames{1}); hold on
        [~,edges3] = histcounts(Rch3);
    else
        [Ich1,Ich2,Ich3] = getIbinsProcessed(cl_bin,ii,'sum','c',{ch1,ch2,ch3},'x','y','c');
        Rch1 = sum(Ich1,3)./Nf; Rch1(isnan(Rch1)) = [];
        Rch2 = sum(Ich2,3)./Nf; Rch2(isnan(Rch2)) = [];
        Rch3 = sum(Ich3,3)./Nf; Rch3(isnan(Rch3)) = [];
        
        subplot(1,3,1); histogram(Rch1,edges1,'Normalization',normalization,'DisplayName',savenames{ii}); xlim([0 1])
        subplot(1,3,2); histogram(Rch2,edges2,'Normalization',normalization,'DisplayName',savenames{ii}); xlim([0 1])
        subplot(1,3,3); histogram(Rch3,edges3,'Normalization',normalization,'DisplayName',savenames{ii}); xlim([0 1])
    end

end
subplot(1,3,1); xlabel('ch1/(ch1+ch2+ch3)'); ylabel(normalization)
subplot(1,3,2); xlabel('ch2/(ch1+ch2+ch3)'); ylabel(normalization)
subplot(1,3,3); xlabel('ch3/(ch1+ch2+ch3)'); ylabel(normalization)
legend(string(1:cl_bin.nInput))
p.print(['bf' num2str(bin_factor) '_histogramPeaks_123_all'])

close all