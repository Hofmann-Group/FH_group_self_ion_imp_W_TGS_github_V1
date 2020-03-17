close all
load('SRIM_20MEV_data_new.mat')

% data is imported from the text file manually or from the matlab data file

% Ang is the depth location in angstroms 
% RECOILS is vacancies created by recoils per ion per angstrom

%% plotting the damage profile

h1=figure;
plot(Ang*1e-4,RECOILS)
grid on
ylabel('Vac per ion per Angstrom')
xlabel('Depth (um)')
set(gcf,'color','w');
set(gca,'fontsize',14);
title('20MeV')



%converting all to SI
% taking a median value for the average damage level
% vac is then the vacancies per ion per m

vac=RECOILS(26)*1e10;
depth=Ang*1e-10;

% getting the fluence needed for the relevant dpa

dpa=[0.001 0.01 0.1 1 10];
atom_dens=6.258e28;

% fluence needed in ions/m2
fluence=(atom_dens*dpa)./(vac);



% plot(dpa,fluence)
% in ions mm-2
fluence2=fluence*1e-6;
% figure
% plot(dpa,fluence2)
% set(gcf,'color','w');
% set(gca,'fontsize',14);
% grid on
% xlabel('Dose (dpa)','FontSize',14)
% ylabel('Fluence (ions/mm2) ','FontSize',14)

fluence2
dpa


% filename = 'SRIM_data_sandia_v2.xlsx';
% xlswrite(filename,'30 MeV','Sheet1','A19')
% xlswrite(filename,dpa,'Sheet1','A20')
% xlswrite(filename,fluence2,'Sheet1','A21')
% % % 

%% dpa_depth_plot

% at index 20 is where the dpa is 0.1 , so we normalise to that value 
% making it for 1 dpa - normalising for the 1dpa location
RECOILS2=(RECOILS./RECOILS(26))*1;



%probed depth - 
TG_depth=2.7528/pi;


% close all
h1=figure;
plot(Ang*1e-4,RECOILS2,'b','LineWidth',2)
grid on
ylabel('Damage (dpa)')
xlabel('Depth ({\mu}m)')
set(gcf,'color','w');
set(gca,'fontsize',14);
axis([0 3 0 1.2])
h = vline(TG_depth,'r','Probed Depth');
% title('20MeV')



%% atomic densities of the injected ions for the analysis plot

% first - getting the fluence needed 
% load from the VACANCY file, depth as Ang and RECOILS 

load('SRIM_20MEV_data_new.mat')

% 
% h1=figure;
% plot(Ang*1e-4,RECOILS)
% grid on
% ylabel('Vac per ion per Angstrom')
% xlabel('Depth (um)')
% set(gcf,'color','w');
% set(gca,'fontsize',14);
% title('20MeV')



%converting all to SI
% taking a median value for the average damage level
% vac is then the vacancies per ion per m

vac=RECOILS(15)*1e10;
depth=Ang*1e-10;

% getting the fluence needed for the relevant dpa
load ('helsinki_summary_data_4_4.mat','dose')
dpa=dose;

% dpa=logspace(-4,1,50);
atom_dens=6.258e28;

% fluence needed in ions/m2
fluence=atom_dens*dpa./(vac);

%% load the RANGe - injected ion concentration data  - or it is loaded from the SRIM20MEVdata_new.mat
% load as depth2 and ion_conc - ion_conc is in atomsper cm3/ atoms per cm2
% converting to SI
ion_conc_2=ion_conc*1e2;  % this has units ions per m3 per ions per m2

% h2=figure;
% plot(Ang2*1e-4,ion_conc_2)
% grid on
% ylabel('Ions per m3 per ions per m2')
% xlabel('Depth (um)')
% set(gcf,'color','w');
% set(gca,'fontsize',14);
% title('20MeV')
% 
% 
% depth2=Ang2*1e-10;

% multiplying the ion_conc by the fluences for each dpa
% this gives the ion ranges / concentrations for each dpa implantation 
ion_conc_3= fluence'*ion_conc_2';

% figure
% plot(dpa,fluence)
% 
% 
% figure
% 
% plot(depth2,ion_conc_3(12,:))
% grid on 
% 
% % to summarise the implanted ion conc with dpa, getting the average of the
% % implanted ion concentration in the probed depth of lambdaTG/pi
% 
% lambdaTG= 2.7528e-06;
% pr_depth= lambdaTG/pi ; % probed depth 

% ind=find((depth2/pr_depth)>1,1);
% 
% % mean ion conc upto probed depth 
% 
% mean_ion=mean(ion_conc_3(1:14,1:ind-1),2);
% 
% % dividing by atomic density of tungsten 
% 
% mean_ion_at_fr= mean_ion/6.258e28;
% 
% close all
% figure 
% plot(dpa,mean_ion_at_fr,'bx')
% grid on
% xlabel('Dose (dpa)','FontSize',16)
% ylabel('Implanted Ion Concentration (at. fr.) ','FontSize',16)
% set(gcf,'color','w');
% set(gca,'fontsize',16);
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xticks([1e-6 1e-4 1e-3 1e-2 1e-1 1e-0 10])
% xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})

% save('implanted_ions_data.mat','dpa','mean_ion_at_fr')


% %% plotting ion conc at fr vs depth for 1 dpa - paper 
% 
% % ion concn at fr 
 ion_conc_4=ion_conc_3/6.258e28;
% %probed depth - 
% TG_depth=2.7528/pi;
% 
% 
% figure
% 
% plot(depth2*1e6,ion_conc_4(12,:),'b','LineWidth',2)
% grid on 
% ylabel('Implanted Ions (at. fr.)')
% xlabel('Depth ({\mu}m)')
% set(gcf,'color','w');
% set(gca,'fontsize',14);
% axis([0 3 0 .8e-4])
% h = vline(TG_depth,'r','Probed Depth');


%% combining the damage depth and ions vs depth plots 

fig = figure;
left_color = [ 0 0 1];
right_color = [1 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

% we plot the left axis 

yyaxis left
plot(Ang*1e-4,RECOILS2,'bd-','LineWidth',1)
grid on
ylabel('Damage (dpa)')
xlabel('Depth ({\mu}m)')
set(gcf,'color','w');
set(gca,'fontsize',14);
axis([0 3 0 1.2])

% give the PATCH

ax=gca  % setting the parent axis

% color

rightcolor = [1 1 1];
    leftcolor = [255/255 172/255 137/255];
    
    
    cdata(1,1,:) = leftcolor;
cdata(1,2,:) = rightcolor ;
cdata(1,3,:) = rightcolor;
cdata(1,4,:) =  leftcolor;

p =patch([0 1 1 0], [0 0 1.2 1.2],'k','Parent',ax);

set(p,'CData',cdata, ...
    'FaceColor','interp', ...
    'EdgeColor','none');

p.FaceVertexAlphaData = 0;    % Set constant transparency 
p.FaceAlpha = 'flat' ;          % Interpolate to find face transparency

uistack(p,'bottom') % Put gradient underneath everything else


% now plot the second graph 

hold on
yyaxis right 

plot(depth2*1e6,ion_conc_4(12,:),'mx-','LineWidth',1)
grid on 
ylabel('Implanted Ions (at. fr.)')
xlabel('Depth ({\mu}m)')
set(gcf,'color','w');
set(gca,'fontsize',14);
axis([0 3 0 7e-5])




