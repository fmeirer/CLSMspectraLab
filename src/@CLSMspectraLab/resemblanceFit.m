function out = resemblanceFit(obj,id,thresholdbgk,thresholdal,pathReference)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alumina = green, Zeolite=blue, bgk=black, clay=red, silica=yellow
% Run the code "imagesc(minm)" to see residual plot 
%       ***     main Variables   ***
% imagef == reconstracted image (different parts)
% minm = residul matrix
% outp = output, reconstracted image without tresholding
% n=original data matrix
% n2=normalized data matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1 = csvread(pathReference,1,0); % provide address of CSV file for reference spectra with order as specified in template CSV file

bgk=data1(:,1);
Clay=data1(:,2);
USY=data1(:,3);
Alumina=data1(:,4);
SiO2=data1(:,5);

Clay=normalizes(Clay);
Alumina=normalizes(Alumina);
USY=normalizes(USY);
SiO2=normalizes(SiO2);
bgk=normalizes(bgk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
subplot(2,2,1)
plot(bgk,'k')
hold on
plot(Clay,'r')
plot(USY,'b')
plot(Alumina,'g')
plot(SiO2,'y')
ylim([-0.1 1.5])
legend('bgk','clay','USY','alumina','Silica')
title('reference spectra')
%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)

if any(ismember(obj.input(id).dimLabel,{'z','t'}))
    error('Dimensions ''z'' and ''t'' are not (yet) supported');
end

I = obj.getImageProcessed(id,'reshaped','c','x','y');
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
outp=zeros(binx1,biny1);
energym=zeros(binx1,biny1,5);
imagef=zeros(binx1,biny1,3);
minm=zeros(binx1,biny1);
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


end

function y1=normalizes(y)
ymin=min(y);
ymax=max(y);
y1=(y-ymin)./(ymax-ymin);
end

function y1=normalizesall(y)
ymin=min(y,[],'all');
ymax=max(y,[],'all');
y1=(y-ymin)./(ymax-ymin);
end