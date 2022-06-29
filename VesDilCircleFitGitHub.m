function [ezt,eat,ebt,ealphat,d,cirdata]=VesDilCircleFit(stack,pos,inputparams,basicparams,reg_num)

%%%% Written primarily by Pratik Y. Chhatbar 2013-2014 %%%%%

% stack - 3D Matlab array of the images across the run (e.g. 512x512xNimages)  
% pos - XY coordinates of the region in which pixels will be thresholded and fit to the circle - drawn around each vessel in the calling function
% inputparams - structure with various parameters relevant for the cirle fitting
% basicparams - structure used across files having basic parameters used for the processing - folder names and such 
% reg_num - a parameter to save the fit images with the right region number if more than one region is fit per run.

showindiim = 2; % display fits of individual images: 1: color, 2: gray, 3: gray+black, 4: color+gray ROI in single fig [default], 5: gray+color ROI in single fig. x.yyy will give .yyy delay between the images
% the number value is delay between two consecutive display % for minimum display of 1 ms, say 0.001
showall=0; %set to 1 to plot the fit for each frame - you can close it once it starts and it will only draw the n2save frames if savefitimages=1 - but if you don't want to have to close it for every selection, set this to 0

%%% unpack parameter structures %%%%%%%%%%%%        
thresquant=inputparams.thresquant;  oversampl=inputparams.oversampl;    dosobel=inputparams.dosobel;
ellipsefit=inputparams.ellipsefit;  dotracefit=inputparams.dotracefit;
invertim=inputparams.invertim;% set to 1 if you want to invert the image, e.g., for blank vessel without dye
if invertim==1,    disp('inverting image luminance!!!!!!');  end

analysisDir=basicparams.analysisDir;    

if isfield(inputparams,'threshmeth'),   threshmeth=inputparams.threshmeth;  
else                                    threshmeth='new';
end
if isfield(basicparams,'savefitimages'),savefitimages=basicparams.savefitimages;
else                                    savefitimages=1;
end
n2save=10;%how many fit images to save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wb = waitbar(0,'calculate area');
%%%%%%%%%% create a mask %%%%%%%%%%%%%%%%
sz=size(stack);
szx = 1:sz(2); szy = 1:sz(1); % image pixel coordinates
[szxm,szym] = meshgrid(szx,szy); % mesh that matches actual image
% with oversampl included. set oversampl=1 if not needed.
maskorig = poly2mask(pos(:,1)*oversampl,pos(:,2)*oversampl,sz(1)*oversampl,sz(2)*oversampl);
[moy,mox] = find(maskorig);
moya = moy/oversampl; moxa = mox/oversampl; % actual image pixel coordinates
% maskmodi = zeros(1:max(moy)-min(moy),1:max(mox)-min(mox));
mmx = min(mox):max(mox); mmy = min(moy):max(moy);
mmxa = mmx/oversampl; mmya = mmy/oversampl;
maskmodi = maskorig(mmy,mmx);
osim = false.*maskmodi;
if dosobel %V7
    sobx = [1 2 1]'*[1 0 -1]; soby = sobx';
    mmsx = filter2(sobx,maskmodi);mmsy = filter2(soby,maskmodi);
    mmsxya = logical(abs(mmsx)+abs(mmsy));
    maskmodisi = maskmodi & mmsxya==0;
end
% preallocate for vessel CROSS SECTIONAL AREA
d=zeros(size(stack,3),1);%fit every frame in the stack and extract the ones you want for analysis further down
% preallocate for ellipse parameters
ezt = zeros(size(stack,3),2);  eat = zeros(size(stack,3),1);  ebt = zeros(size(stack,3),1);  ealphat = zeros(size(stack,3),1);%nFrames
if showindiim
    figure(333);%a waste figure so the FOV doesn't get overwritten when figure 451 gets closed 
    showtype = floor(abs(showindiim));
    switch showtype
        case 1
            figure(451); colormap default;
            ddl{1}='w--'; cfl{1}= 'w'; nfig = 1; spc=4;
        case 2
            figure(451); colormap gray;
            ddl{1}='y--'; cfl{1}= 'r';  nfig = 1; spc=4;
        case 3
            figure(451); colormap default; figure(452); colormap gray;
            ddl{1}='w--'; cfl{1}= 'w'; ddl{2}='y--'; cfl{2}= 'r'; nfig = 2; spc=4;
        case 5
            figure(451); colormap gray; nfig = 1; spc=5;
            ddl{1}='y--'; cfl{1}= 'r'; ddl{2}='w--'; cfl{2}= 'w';
        otherwise % 4 is default
            figure(451); colormap default; nfig = 1; spc=5;
            ddl{1}='w--'; cfl{1}= 'w'; ddl{2}='y--'; cfl{2}= 'r';
    end
end
saveclosed=0;
for n=1:size(stack,3)%nFrames
    waitbar(n/size(stack,3),wb,['calculate area : image ' num2str(n) ' of ' num2str(size(stack,3))]);%nFrames
    % threshold the image, and wipe out the mask
    curim = stack(:,:,n);
    if oversampl>1,     osimt = interp2_p(szxm,szym,curim,moxa,moya);
    else                osimt = curim(maskorig);
    end
    osim(maskmodi) = osimt;
    % patch for "horizontal" (old) vs "vertical" (new) thresholding
    if strcmp(threshmeth,'old')% old "horizontal" way
        thresil = quantile(osimt,thresquant);
    else % new "vertical" way
        thresil = min(osimt)+(max(osimt)-min(osimt))*thresquant;
    end
    %thresil = min(osimt)+(max(osimt)-min(osimt))*thresquant;
    imt = osim>thresil;
    if invertim
        imt(maskmodi) = logical(1-imt(maskmodi));
    end
    if dosobel %V7
        imsx = filter2(sobx,imt); imsy = filter2(soby,imt);
        imst = (abs(imsx)+abs(imsy)).*maskmodisi;
        % imst = filter2(sob2,imt).*maskmodi;
    else
        imst = imt.*maskmodi;
    end
    % find the vessel points
    [x2,x1] = find(imst);
    % give it original values
    x2 = mmya(x2);  x1 = mmxa(x1);
    try
        if ellipsefit
            % do ellipse fitting and estimate cross-sectional area
            if dotracefit,          [zt, at, bt, alphat] = fitellipse([x1;x2],'linear', 'constraint', 'trace');
            else                    [zt, at, bt, alphat] = fitellipse([x1;x2],'linear');
            end
            alphat = mod(alphat,2*pi);
            %8/29/14 - should use smaller radius of ellipse and compute vessel diameter instead of multiplying the two radii and
            %computing area - long radius is often unstable
            if at>=bt,      userad=bt;
            else            userad=at;
            end
        else
            % do circle fitting and estimate cross-sectional area
            [zt, at] = fitcircle([x1;x2]);   bt = at;  userad=at;  alphat = 0;
        end        
    catch %will just use previous values if not overwritten - won't work if it's the first frame though
        errstr='!! Fit failed for image %d . Using previous image value!!\n';
        fprintf(errstr,n);         
    end
    ezt(n,:)=zt;   eat(n,1)=at;   ebt(n,1)=bt;   ealphat(n,1)=alphat;
   
    d(n,1) = userad*2; cirdata='diam';
    if ~showall && n>5 && any(findobj('type','figure')==451)%let it draw a few just in case you want to see the first ones
        close(451);% gcf
    end
    %d(n,1) = pi*at*bt; cirdata='area';
    if (showindiim && any(findobj('type','figure')==451)) ...%ismember(451,findobj('type','figure')) % latter will ensure that figure is not drawn once closed 
            || savefitimages
        if ~any(findobj('type','figure')==451) && mod(n,floor(size(stack,3)/n2save))==0 % will save n2save images as samples even after closing the window for every frame fit 
            saveclosed=1;
            switch showtype % have to re-initialize figure 451 if it's been closed
                case 1
                    figure(451); colormap default;
                    ddl{1}='w--'; cfl{1}= 'w'; nfig = 1; spc=4;
                case 2
                    figure(451); colormap gray;
                    ddl{1}='y--'; cfl{1}= 'r';  nfig = 1; spc=4;
                case 3
                    figure(451); colormap default; figure(452); colormap gray;
                    ddl{1}='w--'; cfl{1}= 'w'; ddl{2}='y--'; cfl{2}= 'r'; nfig = 2; spc=4;
                case 5
                    figure(451); colormap gray; nfig = 1; spc=5;
                    ddl{1}='y--'; cfl{1}= 'r'; ddl{2}='w--'; cfl{2}= 'w';
                otherwise % 4 is default
                    figure(451); colormap default; nfig = 1; spc=5;
                    ddl{1}='w--'; cfl{1}= 'w'; ddl{2}='y--'; cfl{2}= 'r';
            end
        end
        tt = linspace(0,2*pi,100);
        for pp = 1:nfig
            if any(findobj('type','figure')==451) %ismember(451,findobj('type','figure'))
                figure(450+pp)
                subplot(1,spc,[1 2]);
                imagesc(curim);
                axis image; hold on; axis ij;
                plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],ddl{pp});
                if showtype>3, freezeColors; end
                if ellipsefit
                    plotellipse(zt,at,bt,alphat,cfl{pp}); hold off;
                else
                    plot(zt(1)+at*cos(tt),zt(2)+at*sin(tt),cfl{pp}); hold off;
                end
            end
            if any(findobj('type','figure')==451) %ismember(451,findobj('type','figure'))
                figure(450+pp)
                subplot(2,spc,spc-1);
                imagesc(mmxa,mmya,curim(mmya(1):mmya(end),mmxa(1):mmxa(end)));
                axis image; hold on; axis ij;
                plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],ddl{pp});
                xlabel(['diameter=' num2str(d(n)) ',a=' num2str(at) ',b=' num2str(bt) ',alpha=' num2str(alphat*180/pi)]);%
                title(analysisDir,'interpreter','none');
                if showtype>3, freezeColors; end
                if ellipsefit
                    plotellipse(zt,at,bt,alphat,cfl{pp}); hold off;
                else
                    plot(zt(1)+at*cos(tt),zt(2)+at*sin(tt),cfl{pp}); hold off;
                end
            end
            if any(findobj('type','figure')==451) %ismember(451,findobj('type','figure'))
                figure(450+pp)
                subplot(2,spc,2*spc-1);
                imagesc(mmxa,mmya,osim);
                axis image; hold on; axis ij;
                plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],ddl{pp});
                if showtype>3, freezeColors; end
                if ellipsefit
                    plotellipse(zt,at,bt,alphat,cfl{pp}); hold off;
                else
                    plot(zt(1)+at*cos(tt),zt(2)+at*sin(tt),cfl{pp}); hold off;
                end
                % colorbar;
            end
            if any(findobj('type','figure')==451) %ismember(451,findobj('type','figure'))
                figure(450+pp)
                subplot(2,spc,spc);
                imagesc(mmxa,mmya,imt);
                axis image; hold on; axis ij;
                plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],ddl{pp});
                if showtype>3, freezeColors; end
                if ellipsefit
                    plotellipse(zt,at,bt,alphat,cfl{pp}); hold off;
                else
                    plot(zt(1)+at*cos(tt),zt(2)+at*sin(tt),cfl{pp}); hold off;
                end
            end
            if any(findobj('type','figure')==451) %ismember(451,findobj('type','figure'))
                figure(450+pp)
                subplot(2,spc,2*spc);
                imagesc(mmxa,mmya,imst);
                axis image; hold on; axis ij;
                plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],ddl{pp});
                if showtype>3, freezeColors; end
                if showtype>4, colormap default; else colormap gray; end
                if ellipsefit
                    plotellipse(zt,at,bt,alphat,cfl{pp}); hold off;
                else
                    plot(zt(1)+at*cos(tt),zt(2)+at*sin(tt),cfl{pp}); hold off;
                end
            end
            if showtype>3
                if any(findobj('type','figure')==451) %ismember(451,findobj('type','figure'))
                    figure(450+pp)
                    subplot(2,spc,spc-2)
                    imagesc(mmxa,mmya,curim(mmya(1):mmya(end),mmxa(1):mmxa(end)));
                    axis image; hold on; axis ij; freezeColors;
                    plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],ddl{3-pp});
                    subplot(2,spc,2*spc-2);
                    imagesc(mmxa,mmya,osim);
                    axis image; hold on; axis ij; freezeColors;
                    plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)],ddl{3-pp});
                    if showtype>4, colormap default; else colormap gray; end
                    if ellipsefit
                        plotellipse(zt,at,bt,alphat,cfl{pp}); hold off;
                    else
                        plot(zt(1)+at*cos(tt),zt(2)+at*sin(tt),cfl{pp}); hold off;
                    end
                end
            end
        end
        if savefitimages
            
            if mod(n,floor(size(stack,3)/n2save))==0        
                try%if window is made small to speed it up this will crash
                    saveas(gcf,[analysisDir '\fitfig_reg' num2str(reg_num) 'frame' num2str(n) '.fig']);
                    saveas(gcf,[analysisDir '\fitfig_reg' num2str(reg_num) 'frame' num2str(n) '.jpg']);
                    saveas(gcf,[analysisDir '\fitfig_reg' num2str(reg_num) 'frame' num2str(n) '.eps'],'psc2');
                catch
                end
            end
        end
        pause(mod(showindiim,1));
        if saveclosed && any(findobj('type','figure')==451)
            close(451);% gcf
        end
    end
end
close(wb)




