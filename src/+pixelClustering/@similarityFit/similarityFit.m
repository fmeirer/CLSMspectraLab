classdef similarityFit < pixelClustering.abstract
    %SIMILARITYFIT Computes the clustering by assigning each pixel to 
    %   Detailed explanation goes here
    
    properties(Constant)
        name = 'similarityFit';
    end
    
    properties
        loadFilepath = []; % filepath to the reference spectra (should be csv file with columns with spectra)
        nClusters = []; % number of clusters
        minPixGroupSize = 2; % minimum number of pixels in group to be considered a particle for in the processed cluster map
        diagonalPixConnected = true; % if true, connectivity of 8; otherwise, 6
    end
    
    properties(Hidden)
    end
    
    methods
        function obj = similarityFit(varargin)
            %SIMILARITYFIT Construct an instance of this class
            %
            % Usage:
            %   obj = similarityFit(varargin) with varargin Name-Value pairs,
            %   e.g. if we set the property 'prop' to 5: 
            %   obj = similarityFit('prop',5).
            if nargin > 0
                nIn = numel(varargin);
                if mod(nIn,2) ~= 0
                    error('Variables can only be parsed as Name-Value pairs')
                end
                for ii = 1:nIn/2
                    try
                        obj.(varargin{2*ii-1}) = varargin{2*ii};
                    catch
                        error(['The property ''' varargin{2*ii-1} ''' is not recognised.'])
                    end
                end
            end
        end
        
        function obj = compute(obj,iS)
            %COMPUTE Loads the reference files and
            %    and stores the cluster map as obj.I
            
            if isempty(obj.loadFilepath)
               [file,path] = uigetfile('*.csv',...
                    'Select reference spectra file', ...
                    'MultiSelect', 'off');
               if file == 0
                   error('No file loaded.')
               end
               filename = fullfile(path,file);
            else
               filename = obj.loadFilepath;
            end
            
            % open files
            refSpectra = readmatrix(filename); % provide address of CSV file for reference spectra with order as specified in template CSV file
            refSpectraNorm = normalizes(refSpectra);

            obj.nClusters = size(refSpectraNorm,2);

            bgk=refSpectra(:,1);
            Clay=refSpectra(:,2);
            USY=refSpectra(:,3);
            Alumina=refSpectra(:,4);
            SiO2=refSpectra(:,5);
            
%             Clay=normalizes(Clay);
%             Alumina=normalizes(Alumina);
%             USY=normalizes(USY);
%             SiO2=normalizes(SiO2);
%             bgk=normalizes(bgk);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%             subplot(2,2,1)
%             plot(bgk,'k')
%             hold on
%             plot(Clay,'r')
%             plot(USY,'b')
%             plot(Alumina,'g')
%             plot(SiO2,'y')
%             ylim([-0.1 1.5])
%             legend('bgk','clay','USY','alumina','Silica')
%             title('reference spectra')
%             %%%%%%%%%%%%%%%%%%%%%%%%
%             subplot(2,2,2)
            
            if any(ismember(iS.dimLabel,{'z','t'}))
                error('Dimensions ''z'' and ''t'' are not (yet) supported');
            end
            
            I = iS.getReducedImage('reshaped','c','x','y');
            I = double(I); % convert to double to compare with reference spectra
            
            % data = bfopen;
            % S=data{1,1};
            % S1=size(S{1,1});
            % n=zeros(length(S),S1(1),S1(2));
            % for i=1:length(S)
            %  n(i,:,:)=S{i,1};
            % end
            % 
            % binx1=ceil(S1(1)*binscale);
            % biny1=ceil(S1(2)*binscale);
            % nnew=zeros(length(S),binx1,biny1);
            % for i=1:length(S)
            %     sss=squeeze(n(i,:,:));
            %     sss2=imresize(sss,binscale);
            %     ind2=size(sss2);
            %     sss2 = reshape(sss2,[1,ind2(1),ind2(2)]);
            %     nnew(i,:,:)=sss2;
            % end
            
            binx1 = size(I,2);
            biny1 = size(I,3);
            
            n=I;
            n2=normalizes(I);
            nave=squeeze(sum(n));
            nave=normalizesall(nave);

            % pre alloc
            outp=zeros(binx1,biny1);
            energym=zeros(binx1,biny1,5);
            imagef=zeros(binx1,biny1,3);
            minm=zeros(binx1,biny1);
            
            residuals = nan(size(I,2),size(I,3),obj.nClusters);
            for ii = 1:obj.nClusters
                residuals(:,:,ii) = sqrt(sum((n2-refSpectraNorm(:,ii)).^2));
            end
            [~,cluster] = min(residuals,[],3);
            for i=1:binx1
                for j=1:biny1

                    Al=(sum((n2(:,i,j)-Alumina).^2))^0.5;
                    bgk=(sum((n2(:,i,j)-bgk).^2))^0.5;
                    Cl=(sum((n2(:,i,j)-Clay).^2))^0.5;
                    USy=(sum((n2(:,i,j)-USY).^2))^0.5;
                    Si=(sum((n2(:,i,j)-SiO2).^2))^0.5;        
                    tot=[Al Cl USy Si bgk];        
                    energym(i,j,:)=tot;        
                    min1=min(tot);
                    ind=find(tot==min1);
                    ind=ind(1);
                    minm(i,j)=min1;
                    if (ind==1) && (nave(i,j)<thresholdal)
                        outp(i,j)=20;
                        imagef(i,j,:)=[0 255 0];
                    elseif ind==2
                        outp(i,j)=50;
                        imagef(i,j,:)=[255 0 0];
                    elseif ind==3
                        outp(i,j)=100;
                        imagef(i,j,:)=[0 0 255];
                    elseif ind==4
                        outp(i,j)=80;  
                        imagef(i,j,:)=[255 255 0];
                    else
                        outp(i,j)=1;
                        imagef(i,j,:)=[0 0 0];
                    end
                    if (nave(i,j)<thresholdbgk) && outp(i,j)~=20
                        outp(i,j)=1;
                        imagef(i,j,:)=[0 0 0];
                    end
                 end
            end
            ind=isinf(energym);
            energym(ind)=5;
            imagesc(outp);
            colorbar
            title('colormap of different part 0=bgk 20=alumina 50=clay 100=zeolite 80=silica')
            subplot(2,2,3)
            imshow(imagef);
            title('resolved image, black=bgk green=alumina/bgk red=clay blue=zeolite yellow=silica')
            subplot(2,2,4)
            ss=mean(n(:,:,:),[1]);
            ss=squeeze(ss);
            imagesc(ss);
            colorbar
            title('image based on intensity profile')
            obj = obj.addImage(clustermap,'x','y'); % 3D is not yet supported
            obj.nClusters = max(obj.I,[],'all'); % starts counting at 0

            function y1=normalizes(y)
                ymin=min(y,[],1);
                ymax=max(y,[],1);
                y1=(y-ymin)./(ymax-ymin);
            end
            
            function y1=normalizesall(y)
                ymin=min(y,[],3);
                ymax=max(y,[],3);
                y1=(y-ymin)./(ymax-ymin);
            end

%             function y1=normalizesall(y)
%                 ymin=min(y,[],'all');
%                 ymax=max(y,[],'all');
%                 y1=(y-ymin)./(ymax-ymin);
%             end
        end

        function obj = setLoadFilenames(obj,filepath)
            % SETLOADFILENAME allows to set the filepath of the reference
            % spectra.
            % If this function is not run, a file selection dialog box will
            % prompt for a file path when running compute.

            [obj.loadFilepath] = deal(filepath);
        end

    end

end