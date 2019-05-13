function [sig,esig]=HgBeamSmearingCorr(part,src)
%Apply a beam smearing correction to the CALIFA Hg mom2 maps
%Usage: HgBeamSmearingCorr(part,*src)
%   part: 1=make the cubes, 2=after convolution in Miriad do BS correction
%   src: optional parameter, if given does only that galaxy, if not does all
%
% 1. Creates a rotation only data cube from the fitted model velocity field
% where at each pixel the Hg line is a delta function with the velocity
% given by the model velocity field and the amplitude is from the cube's
% value at that pixel. Based on COBeamSmearingCorr.m
% BREAK: export to Miriad for convolution to appropriate beam size and make
% moment maps
% 2. From actual mom2 subtract model beam smeared mom2 in quadrature
% 3. If all sources run, write table of average velocity dispersions

cd ../CALIFA

%load geometric parameters
srcdata=parseCsvTable('../EDGE/KinematicParameters.csv');
srclist=srcdata.Name;
MT=parseCsvTable('../EDGE/DETableFinal_2016Dec14.csv');

%src is an optional keyword
%if present, do only the specified source
%if absent, do all sources in srclist
if nargin==2
    istart=strmatch(src,srclist);
    iend=istart;
else
    istart=1;
    iend=length(srclist);
end

if part==1
    %create simulated data cube
    
    %create list of flags for whether there is a good RC and fit, will be used in part 3
    flaglist=zeros(size(srclist)); 
    sig=NaN;
    esig=NaN;
        
    for i=istart:iend
        src=char(srclist(i));
        disp(src);
        
        %get geometric parameters
        ra=MT.ledaRA(i)*15; %degrees
        dec=MT.ledaDE(i);
        pa=srcdata.PA(i);
        inc=srcdata.Inc(i);
        vsys=srcdata.haVsys(i);
        xoff=srcdata.haXoff(i); %arcsec
        yoff=srcdata.haYoff(i); %arcsec
        
        %load moment maps and coordinates
        [~, rf, rr] = Cdefvar(src,pa,inc,ra,dec,xoff,yoff);
        
        %get redshift
        z=MT.caZgas(i);
        lamHg=4340.47; %AA
        lamHgobs=lamHg*(1+z); %AA
        
        %import fitted CO rotation curve
        RC=importdata(['fits/data/',src,'.ifitinfo.txt']);
        if isfield(RC,'data')==1 & strcmp(src,'ARP220')==0 & strcmp(src,'UGC03253')==0 & strcmp(src,'UGC10043')==0
            %grab rotation curve points from RC file
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
            if round(chi2r) > 1e4 %model isn't a great fit
                vfmod=zeros(size(rf.data))+vsys; %model is constant with V=Vsys
            else 
                %make the model velocity field
                vfmod=Cmkvelfield(rf,pa,inc,ra,dec,xoff,yoff,R,Vrotm,Vradm,Vsysm);
                vfmod=vfmod+vsys; %Add Vsys back in
                flaglist(i)=1;
            end
            
            %plot data and model RC
            figure(8); clf; hold on;
            errorbar(R,Vrot,eVrot,'ro-','MarkerFaceColor','r','Linewidth',1);
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
        
        %convert to optical convention from relativistic since cubes are in
        %optical!!
        c=2.99792458E5; %km/s
        vfmod=c*(sqrt((1+vfmod/c)./(1-vfmod/c))-1);
        
        %use model velocity field to create a rotation only data cube
        
        %get channel velocities and nchans
        vstart=vsys-1000; %starting velocity for simcube
        vstep=0.7/lamHgobs*c; %km/s, spectral res
        %vstep=1;
        vend=vsys+1000;
        chans=vstart:vstep:vend;
        nchan=length(chans);
        
        %load V1200 FE to get final beam size
        hgfe=rfits(['flux_elines_V1200/flux_elines.',src,'.cube.fits']);
        if isfield(hgfe,'fwhm')==1
            finbeam=hgfe.fwhm;
        else
            finbeam=2.5;
        end
        
        %load co cube as spatial template
        cocube=rfits(['../EDGE/DataCubes/',src,'.co.cmmsk.fits']);
        sx=cocube.numpt(1);
        sy=cocube.numpt(2);
        
        %use cocube header for Miriad
        %CALIFA headers are not read well by Miriad...
        rotcube=cocube;
        rotcube.data=zeros(sx,sy,nchan);
        
        %update header accordingly
        rotcube.numpt(3)=nchan;
        rotcube.crval(3)=vstart*1000; %m/s
        rotcube.cdelt(3)=vstep*1000; %m/s
        rotcube.lstart=vstart;
        rotcube.lstep=vstep;
        rotcube.lwidth=vstep;
        rotcube.x{3}=chans*1000; %m/s       
        
        %define new bmaj and bmin < 1 pixel < 1" for convolution
        rotcube.bmaj=0.5/3600; 
        rotcube.bmin=rotcube.bmaj;
        rotcube.bpa=finbeam; %STORE FINBEAM IN BPA TO GRAB IN MIRIAD
        
        
        %at each spaxel, put a delta function at the channel corresponding
        %to the model rotation velocity field
        specres=2.3/2.355/lamHgobs*c;
        for j=1:sx
            for k=1:sy
                [~,chan]=min(abs(vfmod(j,k)-chans));
                %rotcube.data(j,k,chan)=max(cocube.data(j,k,:));
                %make a Gaussian with mean=chan, sigma=2.3/2.355
                gauss=exp(-0.5*(chans-vfmod(j,k)).^2./specres^2);
                rotcube.data(j,k,:)=gauss;
            end
        end
        
        %save cube
        wfits(rotcube,['DataCubes_SimInfRes/',src,'.hg.infres.fits']);  
    end
    
    %save flaglist
        if nargin==1
            fp=fopen('CALIFA_Hg_VelDisp_flags.csv','w');
            fprintf(fp,'Name,flag\n');
            for j=1:length(srclist)
                fprintf(fp,'%s,%1.f\n',char(srclist(j)),flaglist(j));
            end
        end
end

%BREAK: port to Miriad to do the convolution and moment map generation

if part==2
    
    if istart==1
        fp=fopen('HgVelDisp_all.csv','w');
        fprintf(fp,'Name,HgVelDisp,eHgVelDisp\n');
    else
        fp=1;
    end
    
    %remove instrumental and beam smearing from mom2 map   
    for i=istart:iend
        src=char(srclist(i));
        
        %load mom2
        mom2=rfits(['flux_elines_V1200/mom_remap/',src,'.Hg.mom2.fits']);
        mom2bscorr=mom2;
        m2=mom2bscorr.data;
        m2(m2<=0)=NaN;
        
        %remove instrumental LW, 2.3AA FWHM
        %m2=real(sqrt(m2.^2-2.3^2));
        
        %convert from FWHM to sigma
        m2=m2/2.355;
        
        %convert to velocity
        c=2.99792458E5; %km/s
        z=MT.caZgas(i);
        lamHg=4340.47; %AA
        lamHgobs=lamHg*(1+z); %AA
        m2=c*m2/lamHgobs;
        m2orig=m2;
        
        %load mom2 from beam smearing and remove its instrumental LW
        mom2bs=rfits(['Moments_BeamSmeared/',src,'.mom2.bs.fits']);
        m2bs=mom2bs.data;
        m2bs(isnan(m2bs))=0;
        %m2bs=real(sqrt(m2bs.^2-(1/2.355)^2));
        %m2bs=real(sqrt(m2bs.^2-(67.5)^2));
        
        %remove beam smearing from mom2
        m2(m2==0)=NaN;
        rmbs=sqrt(m2.^2-m2bs.^2);
        %m2=real(sqrt(m2.^2-m2bs.^2));
        %im2=imag(sqrt(m2.^2-m2bs.^2));
        m2=real(rmbs);
        im2=imag(rmbs);
        mom2bs.data=m2;
        m2nomask=m2;
        
        %write to file
        wfits(mom2bs,['flux_elines_V1200/mom2_bscorr/',src,'.Hg.bs.mom2.fits']);
        
        %mask to flat region
        [m2,~]=MaskToFlat(src,mom2bs,2,1);
        m2(m2==0)=NaN;
        mom2bs.data=m2;
        
        %write to file
        wfits(mom2bs,['flux_elines_V1200/mom2_bscorr_masked/',src,'.Hg.bs.mom2.fits']);
 
        
        %plot comparison of mom2 maps
        [sig,esig]=pltmom2bscorrHg(src);
        
        %plot BS corr mom2 with imaginary part as -velocity disp
        im2=-im2;
        im2(isnan(im2))=0;
        m2nomask(isnan(m2nomask))=0;
        m2im=m2nomask+im2;
        m2im(m2im==0)=NaN;
        medreal=median(m2im(m2im>0),'omitnan');
        medimag=median(m2im(m2im<0),'omitnan');
        
        %mask to flat region to compare
        mom2im=mom2bs; mom2im.data=m2im;
        [m2imm,~]=MaskToFlat(src,mom2im,2,1);
        medrealm=median(m2imm(m2imm>0),'omitnan');
        medimagm=median(m2imm(m2imm<0),'omitnan');
        
        figure(2); clf;
        cmap=[1 1 1; flipud(othercolor('RdYlBu8'))];
        colormap(cmap);
        imagesc(m2im');
        set(gca,'YDir','normal','FontSize',14); axis image;
        c=colorbar;
        ylabel(c,'Velocity Dispersion (km/s)');
        cax=median(abs(caxis));
        caxis([-cax cax]);
        print(['Moments_BeamSmeared/Imag/',src,'.m2.imag.png'],'-dpng');
        
        figure(3); clf; hold on;
        histogram(m2im,'FaceColor',nicegray,'binwidth',5);
        histogram(m2imm,'FaceColor','b','binwidth',5);
        set(gca,'box','on','fontsize',14);
        xlabel('Velocity Dispersion (km/s');
        ylabel('Number of pixels');
        title(src);
        legend(sprintf('All Pixels: Medians = %1.f, %1.f km/s',medreal,medimag),...
            sprintf('Masked Regions: Medians = %1.f, %1.f km/s',medrealm,medimagm));
        print(['Moments_BeamSmeared/Imag/',src,'.m2.imag.hist.png'],'-dpng');
        
        
        %write to file
        fprintf(fp,'%s,%.1f,%.1f\n',src,sig,esig);
    end  
    if istart==1
        fclose(fp);
    end
end

end