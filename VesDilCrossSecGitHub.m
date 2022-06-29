function [d,figh]=VesDilCrossSec(stack,pos,N,inputparams,basicparams,figure_master_print)

%%%%% Written primarily by Prakash Kara 2011-2012 %%%%%
% stack - 3D Matlab array of the images across the run (e.g. 512x512xNimages)  
% pos - XY coordinates of the ends of the line segment drawn across the vessel - diameter is computed across this line  
% N - cross section number - just used for labeling images when multiple cross sections are drawn per run  
% inputparams - structure with various parameters relevant for the diameter computation  
% basicparams - structure used across files having basic parameters used for the processing - folder names and such 
% figure_master_print - a parameter to determine if images will be drawn/saved  


%%% unpack parameter structures %%%%%%%%%%%% - most are only used in subfunctions
baseName=inputparams.baseName;      Nlines=inputparams.Nlines;
if isfield(inputparams,'lb')
    lb=inputparams.lb;                  ls=inputparams.ls;
    bkgframes=inputparams.bkgframes;    stimframes=inputparams.stimframes;
end
method=inputparams.method;          oversampl=inputparams.oversampl;
invertim=inputparams.invertim; % set to 1 if you want to invert the image, e.g., for blank vessel without dye
if invertim==1,    disp('inverting image luminance!!!!!!');  end

analysisDir=basicparams.analysisDir;

figure_crossection_print=1*figure_master_print;     figure_crossection_close=0;
figure_brightness_print=1*figure_master_print;      figure_brightness_close=0;

NaddLines=(Nlines-1)/2;  % number of additional lines; total Nlines= NaddLines*2+1

%calculate crossection line
le=sqrt((pos(1,1)-pos(2,1)).^2+(pos(1,2)-pos(2,2)).^2);
ler=round(le);
pos(2,1)=pos(1,1)+(pos(2,1)-pos(1,1))*ler/le;
pos(2,2)=pos(1,2)+(pos(2,2)-pos(1,2))*ler/le;
posX=linspace(pos(1,1),pos(2,1),ler*oversampl+1);
posY=linspace(pos(1,2),pos(2,2),ler*oversampl+1);
dY=(posX(1)-posX(2))*oversampl;
dX=(posY(1)-posY(2))*oversampl;
sh=(0:(NaddLines*2))-NaddLines;

% for each frame calculate profile, find max point, remove points around max point, find second max
% calculate distance between two maxs
figh(1)=figure;
hold on;
title([baseName '  ' num2str(N)],'interpreter','none');
xlabel('Vessel crossection');
sz=size(stack);
[x,y] = meshgrid(1:sz(2),1:sz(1));
% preallocate for vessel diameter
d=zeros(size(stack,3),1);
%preallocate for reponse intensity
maxI1=zeros(size(stack,3),1);
maxI2=zeros(size(stack,3),1);
midI=zeros(size(stack,3),1);
endI1=zeros(size(stack,3),1);
endI2=zeros(size(stack,3),1);

%**** get vessel profile and diameter for each frame *****
tmp_test=zeros(size(stack,3),length(sh));

for n=1:size(stack,3)
    if strcmp(method,'lumen')%chan==2
        ZI=zeros(length(sh),length(posX)-1);
    else
        ZI=zeros(length(sh),length(posX));
    end
    for n2=1:length(sh) %do for satellite crossections too
        %interp2_mod - interp2 accelerated for specific use here
        if strcmp(method,'lumen')%chan==2 %if channel 2 then use diff
            ZII = interp2_mod(x,y,stack(:,:,n),posX+dX*sh(n2),posY-dY*sh(n2));
            ZI(n2,:)=diff(ZII);
            spl=round(length(ZI(n2,:))/2);
            ZI(n2,spl:end)=-1*ZI(n2,spl:end);
        else
            ZI(n2,:) = interp2_mod(x,y,stack(:,:,n),posX+dX*sh(n2),posY-dY*sh(n2));
        end
    end

    if invertim
        ZI_orig=ZI;
        ZI=ZI_orig-(ZI_orig*2);
    end
    
    %save intensity
    midI(n)=mean(ZI(:,round(length(ZI)/2)));
    endI1(n)=mean(ZI(:,1));
    endI2(n)=mean(ZI(:,end));


    %plot profiles
    if exist('lb','var')
        if mod(n,lb+ls)==mod(bkgframes(end),lb+ls)
            %plot at the end of blanks
            plot((0:(length(ZI(NaddLines+1,1:end))-1))/oversampl,ZI(NaddLines+1,1:end),'-b');
        end
        if mod(n,lb+ls)==mod(stimframes(end),lb+ls)
            %plot at the end of stim
            plot((0:(length(ZI(NaddLines+1,1:end))-1))/oversampl,ZI(NaddLines+1,1:end),'-r');
        end
    end
    maxI1(n)=0;
    maxI2(n)=0;

    for n2=1:length(sh)
        [ma ind1]=max(ZI(n2,:));
        maxI1(n)=maxI1(n)+ma;

        done=0;
        removerange=1;
        while ~done
            removepoints1=max(ind1-removerange*oversampl,1);
            removepoints2=min(ind1+removerange*oversampl,length(ZI));
            ZI(n2,removepoints1:removepoints2)=0;
            %look for second max
            [ma ind2]=max(ZI(n2,:));
            %calculate diameter
            dtmp=abs(ind1-ind2)/oversampl;
            if dtmp>(removerange+2)
                done=1;
            else
                removerange=removerange+3;
                if sum(ZI(n2,:))==0
                    display(['Possible error in diameter calculation , frame - ' num2str(N)])
                    done=1;
                end
            end
        end

        maxI2(n)=maxI2(n)+ma;
        tmp_test(n,n2)=dtmp;
    end
    maxI1(n)=maxI1(n)/length(sh);
    maxI2(n)=maxI2(n)/length(sh);
end

saveas(figh(1),[analysisDir '\crossections_' num2str(N) '.eps'], 'psc2');
saveas(figh(1),[analysisDir '\crossections_' num2str(N) '.fig']);
saveas(figh(1),[analysisDir '\crossections_' num2str(N) '.jpg']);
if figure_crossection_print==1,        print(figh(1));    end
if figure_crossection_close==1,        close(figh(1));    end

meanSTD=mean(std(tmp_test,0,2));
for n=1:size(stack,3)
    dset=tmp_test(n,:);
    dset(abs(dset-median(dset))>5*meanSTD)=[]; % do not use outliers for diameter calculation
    d(n)=mean(dset);
end


% plot brightness
figh(2)=figure;
hold on;
title([baseName '  ' num2str(N)],'interpreter','none');
plot(maxI1,'r');
plot(maxI2,'r');
plot(midI,'b');
plot(endI1,'k');
plot(endI2,'k');
mmm=[maxI1 maxI2];
mean_wall=mean(mmm(:));
kk=[endI1 endI2];
mean_outside=mean(kk(:));
mean_inside=mean(midI(:));
xlabel(['mean wall brightness=' num2str(mean_wall) '  wall/inside=' num2str(mean_wall/mean_inside) '    wall/outside=' num2str(mean_wall/mean_outside)]);


saveas(figh(2),[analysisDir '\brightness_' num2str(N) '.eps'], 'psc2');
saveas(figh(2),[analysisDir '\brightness_' num2str(N) '.fig']);
saveas(figh(2),[analysisDir '\brightness_' num2str(N) '.jpg']);
if figure_brightness_print==1,        print(figh(2));    end
if figure_brightness_close==1,        close(figh(2));    end


%plot diameter all frames
% section is commented (needed for test only)
%     if NaddLines>0
%         figure;
%         plot(tmp_test);
%         xlabel('Frames');
%         ylabel('Diameter [pix]');
%         title('All sub-crossections');
%     end