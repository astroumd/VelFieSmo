function COBeamSmearingCorr(part,mtype,src)
%Apply a beam smearing correction to the EDGE CO mom2 maps
%Usage: COBeamSmearingCorr(part,mtype,*src)
%   part: 1=make the cubes, 2=after convolution in Miriad do BS correction,
%3=derive galaxy averages and write to file
%   mtype: optional parameter, which kind of making to use: 0: dilated,
%   1: smoothed, 2:gaussian; default is 2  
%   src: optional parameter, if given does only that galaxy, if not does all
%
% 1. Creates a rotation only data cube from the fitted model velocity field
% where at each pixel the CO line is a delta function with the velocity
% given by the model velocity field and the amplitude is from the cube's
% value at that pixel. Based on SimConvLinewidth.m
% BREAK: export to Miriad for convolution to appropriate beam size and make
% moment maps
% 2. From actual mom2 subtract model beam smeared mom2 in quadrature
% 3. If all sources run, write table of average velocity dispersions

cd ../EDGE

%load geometric parameters
srcdata=parseCsvTable('KinematicParameters.csv');
srclist=srcdata.Name;
MT=parseCsvTable('DETableFinal_2016Dec14.csv');

%src is an optional keyword
%if present, do only the specified source
%if absent, do all sources in srclist
if nargin==3
    istart=strmatch(src,srclist);
    iend=istart;
else
    istart=1;
    iend=length(srclist);
end

%pick the mask type to use (see header)
%if absent, use Gaussian
if nargin < 2
    mtype=2;
end


if part==1
    %create simulated data cube
    
    %create list of flags for whether there is a good RC and fit, will be used in part 3
    flaglist=zeros(size(srclist));     
        
    for i=istart:iend
        src=char(srclist(i));
        
        %get geometric parameters
        ra=MT.ledaRA(i)*15; %degrees
        dec=MT.ledaDE(i);
        pa=srcdata.PA(i);
        inc=srcdata.Inc(i);
        vsys=srcdata.coVsys(i);
        xoff=srcdata.coXoff(i); %arcsec
        yoff=srcdata.coYoff(i); %arcsec
        
        %load moment maps and coordinates
        [~, ~, rf, rr] = defvar(src,pa,inc,ra,dec,xoff,yoff);
        
        
        %import fitted CO rotation curve
        RC=importdata(['fits/data/',src,'.fitinfo.txt']);
        if isfield(RC,'data')==1
            R=RC.data(:,1);
            Vrot=RC.data(:,2);
            eVrot=RC.data(:,3);
       
            %make dummy Vrad and Vsys components
            Vrad=zeros(size(Vrot)); eVrad=Vrad;
            Vsys=Vrad; eVsys=Vrad;

            %make fit model RC using Persic+96
            [Vrotm,Vradm,Vsysm,~,R,Vrot,eVrot,~,~,~,~,~,~,~,~]=modelfit(src,2,R,Vrot,eVrot,Vrad,eVrad,Vsys,eVsys,rr);

            %check goodness of fit using reduced Chi2
            [~,chi2r]=chi2model(Vrot,eVrot,Vrotm,3);
            if round(chi2r) > 1E2 %model isn't a great fit
                vfmod=zeros(size(rf.data))+vsys; %make constant with V=Vsys
            else 
                %make the model velocity field
                if strcmp(src,'UGC10043')==1, inc=85; end %this galaxy is edge-on, reduce inc so the script doesn't break
                vfmod=mkvelfield(rf,pa,inc,R,ra,dec,xoff,yoff,Vrotm,Vradm,Vsysm);
                vfmod=vfmod+vsys; %add Vsys back in
                flaglist(i)=1;
            end
            
            %plot data and model RC
            figure(8); clf; hold on;
            errorbar(R,Vrot,eVrot,'bo-','MarkerFaceColor','b','Linewidth',1);
            plot(R,Vrotm,'k-','Linewidth',1);
            xlabel('Radius (")');
            ylabel('V_{rot} (km/s)');
            title(sprintf('%s: X^2_r=%.1f',src,chi2r));
            legend('Data','Persic Model','Location','SouthEast');
            set(gca,'box','on','FontSize',14,'XMinorTick','on','YMinorTick','on');
            print(['rotmodels/comp_Vrot/',src,'.modelcomp.vrot.eps'],'-depsc');
        
        else
            vfmod=zeros(size(rf.data))+vsys;
        end
        
        %convert to radio convention from relativistic since cubes are in
        %radio!!
        c=2.99792458E5; %km/s
        vfmod=c*(1-sqrt((1-vfmod/c)./(1+vfmod/c)));
        
        %use model velocity field to create a rotation only data cube
        %load actual datacube
        cube=rfits(['DataCubes/',src,'.co.cmmsk.fits']);
        rotcube=cube;
        rotcube.data=zeros(size(cube.data));
        
        %define new bmaj and bmin < 1 pixel < 1" for convolution
        rotcube.bmaj=0.5/3600; 
        bratio = cube.bmaj*3600/0.5;
        rotcube.bmin=rotcube.bmin/bratio; %keep beam aspect ratio
        
        %get velocity corresponding to each channel
        chanvel=rotcube.lstart:rotcube.lstep:rotcube.lstart+rotcube.lstep*rotcube.numpt(3)-1;
        
        %at each spaxel, put a delta function at the channel corresponding
        %to the model rotation velocity field
        sx=rotcube.numpt(1);
        sy=rotcube.numpt(2);
        for j=1:sx
            for k=1:sy
                [~,chan]=min(abs(vfmod(j,k)-chanvel));
                rotcube.data(j,k,chan)=cube.data(j,k,chan);
                %if rotcube.data(j,k,chan) < 0
                %    rotcube.data(j,k,chan)=0;
                %end
            end
        end
        
        %save cube
        wfits(rotcube,['DataCubes_SimInfRes/',src,'.co.infres.fits']);  
    end
    
    %save flaglist
        if nargin==1
            fp=fopen('EDGE_CO_VelDisp_flags.csv','w');
            fprintf(fp,'Name,flag\n');
            for j=1:length(srclist)
                fprintf(fp,'%s,%1.f\n',char(srclist(j)),flaglist(j));
            end
        end
end

%BREAK: port to Miriad to do the convolution and moment map generation

if part==2
    %remove instrumental and beam smearing from mom2 map   
    for i=istart:iend
        src=char(srclist(i));
        
        %load mom2
        if mtype==1
            mom2=rfits(['Smoothed/mom2/',src,'.co.de20_smo.mom2.fits']);
        elseif mtype==0
            if exist(['Dilated/mom2/',src,'.co.de20_dil.mom2.fits'],'file')==2
                mom2=rfits(['Dilated/mom2/',src,'.co.de20_dil.mom2.fits']);
            else
                mom2=rfits(['Smoothed/mom2/',src,'.co.de20_smo.mom2.fits']);
            end
        else
            mom2=rfits(['Gauss_20_radio/',src,'.co.de20_smo.fwhm.fits']);
        end
        mom2bscorr=mom2;
        m2=mom2bscorr.data;
        
        if mtype==2
            m2=m2/2.355; %convert gaussian FWHM to sigma
        end
        
        %load mom2 from beam smearing
        mom2bs=rfits(['Moments_BeamSmeared/',src,'.mom2.bs.fits']);
        m2bs=mom2bs.data;
        %m2bs(isnan(m2bs))=0;
        
        %remove beam smearing from mom2 in quadrature
        m2=real(sqrt(m2.^2-m2bs.^2));
        m2(m2==0)=NaN;
        mom2bs.data=m2;
        
        %write to file
        if mtype==1
            wfits(mom2bs,['Smoothed/mom2_bscorr/',src,'.co.de20_smo.bs.mom2.fits']);
        elseif mtype==0
            wfits(mom2bs,['Dilated/mom2_bscorr/',src,'.co.de20_dil.bs.mom2.fits']);
        else
            wfits(mom2bs,['Gauss_20_radio/mom2_bscorr/',src,'.co.de20_gauss.bs.mom2.fits']);
        end
        
        %mask to flat region
        [m2,~]=MaskToFlat(src,mom2bs,0,1);
        m2(m2==0)=NaN;
        mom2bs.data=m2;
        
        %write to file
        if mtype==1
            wfits(mom2bs,['Smoothed/mom2_bscorr_masked/',src,'.co.de20_smo.bs.mom2.fits']);  
        elseif mtype==0
            wfits(mom2bs,['Dilated/mom2_bscorr_masked/',src,'.co.de20_dil.bs.mom2.fits']);  
        else
            wfits(mom2bs,['Gauss_20_radio/mom2_bscorr_masked/',src,'.co.de20_gauss.bs.mom2.fits']);
        end
        
        %plot comparison of mom2 maps
        [sig,esig]=pltmom2bscorr(src,mtype);
    end  
end



%%%%%%%%% PART 3 %%%%%%%%%%%%%%
if part==3
    %derive galaxy averages for non-beam smearing corrected, beam smearing
    %corrected, and beam smearing corrected and inner 2*beam masked
    fp=fopen('EDGE_VelDisp.txt','w');
    fprintf(fp,'Version: %s\n',datestr(now));
    fprintf(fp,'COVelDisp:\n   Average CO velocity dispersion (sigma, not FWHM) from Gaussian fit to CO line mom2, beam-smearing corrected, excluding pixels with radii < 2*bmaj and > R_max,CO (km/s)\n');
    fprintf(fp,'eCOVelDisp:\n   Standard deviation of the CO velocity dispersions as above (km/s)\n');
    fprintf(fp,'HgVelDisp:\n   Average Hgamma velocity dispersion (sigma, not FWHM), beam-smearing corrected, instrumental linewidth removed, excluding pixels with radii < 2*bmaj and > R_max,CO (km/s)\n');
    fprintf(fp,'eHgVelDisp:\n   Standard deviation of the Hgamma velocity dispersions as above (km/s)\n');
    fclose(fp);

    fp=fopen('EDGE_VelDisp.csv','w');
    fprintf(fp,'Name,COVelDisp,eCOVelDisp,HgVelDisp,eHgVelDisp\n');
    
    %fpa=fopen('EDGE_CO_VelDisp_all.txt','w');
    %fprintf(fpa,'Version: %s\n',datestr(now));
    %fprintf(fpa,'COVelDispBSCorr:\n   Average velocity dispersion (sigma, not FWHM) from dilated-mask mom2, beam-smearing corrected, instrumental linewidth removed (km/s)\n');
    %fprintf(fpa,'eCOVelDispBSCorr:\n   Standard deviation of the velocity dispersions as above (km/s)\n');
    %fprintf(fpa,'COVelDispBSCorrMasked:\n   Average velocity dispersion (sigma, not FWHM) from dilated-mask mom2, beam-smearing corrected, instrumental linewidth removed, excluding pixel with radii < 2*bmaj and > R_max,CO (km/s)\n');
    %fprintf(fpa,'eCOVelDispBSCorrMasked:\n   Standard deviation of the velocity dispersions as above (km/s)\n');
    %fprintf(fpa,'flag:\n   Flag indicating whether the rotation smearing is corrected for (1) or not (0), galaxies without robust CO rotation curves or bad smooth model fits do not have the rotation smearing corrected\n');
    %fclose(fpa);

    %if mtype==0
    %   fpa=fopen('EDGE_CO_VelDisp_all_dilated.csv','w');
    %elseif mtype==1
    %    fpa=fopen('EDGE_CO_VelDisp_all_smoothed.csv','w');
    %else
    %    fpa=fopen('EDGE_CO_VelDisp_all_gauss.csv','w');
    %    fpdef=fopen('EDGE_CO_VelDisp_all.csv','w');
    %end
    %fprintf(fpa,'Name,COVelDisp,eCOVelDisp,COVelDispMasked,eCOVelDispMasked,flag\n');
   % 
    %if mtype==2
    %    fprintf(fpdef,'Name,COVelDisp,eCOVelDisp,COVelDispMasked,eCOVelDispMasked,flag\n');
    %end
    
    flaglist=parseCsvTable('EDGE_CO_VelDisp_flags.csv');
    flaglist=flaglist.flag;
    
    
    
    for i=1:length(srclist)
        %get source parameters
        src=char(srclist(i));

        %load original mom2
        %mom2=rfits(['Smoothed/mom2/',src,'.co.de20_smo.mom2.fits']);
        %m2=mom2.data;

        %remove instrumental linewidth, 20 km/s
        %instlw=20;
        %m2full=real(sqrt(m2.^2-instlw^2));
        
        %get average velocity dispersion of non-BS corr mom2
        %VelDispFull=mean(reshape(m2full,[1,numel(m2full)]),'omitnan');
        %eVelDispFull=std(reshape(m2full,[1,numel(m2full)]),'omitnan');
        
        %load BS corr mom2 masked 
        %if mtype==1
        %    mom2=rfits(['Smoothed/mom2_bscorr_masked/',src,'.co.de20_smo.bs.mom2.fits']);
        %elseif mtype==0
        %    mom2=rfits(['Dilated/mom2_bscorr_masked/',src,'.co.de20_dil.bs.mom2.fits']);
        %else
        %    mom2=rfits(['Gauss_20_radio/mom2_bscorr_masked/',src,'.co.de20_gauss.bs.mom2.fits']);
        %end
        %m2=mom2.data;
        
        %get average velocity dispersion of non-BS corr mom2
        %VelDisp=mean(reshape(m2,[1,numel(m2)]),'omitnan');
        %eVelDisp=std(reshape(m2,[1,numel(m2)]),'omitnan');
        
        %load BS corr mom2
        %if mtype==1
        %    mom2=rfits(['Smoothed/mom2_bscorr/',src,'.co.de20_smo.bs.mom2.fits']);
        %elseif mtype==0
        %    mom2=rfits(['Dilated/mom2_bscorr/',src,'.co.de20_dil.bs.mom2.fits']);
        %else
        %    mom2=rfits(['Gauss_20_radio/mom2_bscorr/',src,'.co.de20_gauss.bs.mom2.fits']);
        %end
        %m2=mom2.data;
        
        %get average velocity dispersion of non-BS corr mom2
        %VelDispa=mean(reshape(m2,[1,numel(m2)]),'omitnan');
        %eVelDispa=std(reshape(m2,[1,numel(m2)]),'omitnan')/sqrt(numel(m2(~isnan(m2))));
        
        [VelDisp,eVelDisp]=pltmom2bscorr(src,2);
        
        %also load in Hg velocity dispersions (sigma)
        vddat=parseCsvTable('../CALIFA/HgVelDisp_Gauss.csv');
        ix=strmatch(src,vddat.Name);
        sigHg=vddat.HgVelDisp(ix);
        esigHg=vddat.eHgVelDisp(ix);
        
        if isempty(strmatch(src,vddat.Name))==1
            sigHg=NaN;
            esigHg=NaN;
        end
        
        %print all to csv file
        %fprintf(fpa,'%s,%.1f,%.1f,%.1f,%.1f,%1.f\n',...
        %src,VelDispa,eVelDispa,VelDisp,eVelDisp,flaglist(i));
        
    %if mtype==2
        %print all to csv file
     %   fprintf(fpdef,'%s,%.1f,%.1f,%.1f,%.1f,%1.f\n',...
      %  src,VelDispa,eVelDispa,VelDisp,eVelDisp,flaglist(i));
    %end
    
        %now print to BS corr file if flag == 1
        if flaglist(i)==1 && isnan(VelDisp)==0
             fprintf(fp,'%s,%.1f,%.1f,%.1f,%.1f\n',...
            src,VelDisp,eVelDisp,sigHg,esigHg);
        end
    end
    
    fclose(fp);
    %fclose(fpa);
    %if mtype==2
    %    fclose(fpdef);
    %end
    
    !cp EDGE_VelDisp.txt ../EDGE-CALIFA/EDGE_VelDisp.txt
    !cp EDGE_VelDisp.csv ../EDGE-CALIFA/EDGE_VelDisp.csv

end

end
