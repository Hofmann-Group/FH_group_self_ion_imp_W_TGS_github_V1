clc
clear all
close all

%% This is the analysis of the TC vs DOSE data 

% make sure it has the correct links to the thermal diffusivity vs dose
% data set and the MD +TEM defect density data set 

%% importing data
load ('helsinki_summary_data_4_4_archive.mat','TC','std_TC','dose')
%close all

% renamng TC as TD since it is thermal diffusivity to avoid confusion 
TD=TC;
std_TD=std_TC;

clear TC std_TC


% dose(1)=0.00001;  % setting the unimplanted dose to 0
% TD(1)=6.85e-5;

% this equation was derived for the atomic fraction of vacancies
% assuming all a re poitn defects 


% need to check the valuewe used for electronic heat capacity and
% scattering time 
%Cv=1.7e-7*(6.85./TD -1);

% from the new corrected calculation

% TD=6.5e-5;
% 
% te=21.4e-15; % felix calculation - get a ref 
%te=22.918e-15; % same calculation with 2 terms only

te=21.929e-15; % by inverting TC expression with our unimp data



Ce=26208; % from masons paper apendix and confirmed with the other phys rev paper 
vf=9.5e5; % mason paper 
rho=19.25e3;  %online/ref in transfer  
Cp=132; % ref in transfer 
% sigmav=7.6e15; % felix supplementary
% sigmai=21.6e15; 

% sigmav=2.39e15; % felix supplementary
% sigmai=2.99e15; 


% sigmav=7.64419e15; % latest - without regel limit and with temp adjustment
% sigmai=21.64419e15;     % these two are just sigma 0 adjsuted for temp effects 


% with our calculation and new lorensz number 
sigmav=6.05e15;
sigmai=17.3e15;


Cv=zeros(1,length(TD));

std_Cv= zeros(1,length(TD)); % using simple error theory 

for i=1:length(TD)
Cv(i)=((((te*Ce*(vf)^2)/(3*rho*Cp.*TD(i)))-1))/(te*(sigmav+sigmai)-2);

std_Cv(i)=(((((te*Ce*(vf)^2)/(3*rho*Cp.*(TD(i)).^2))))/(te*(sigmav+sigmai)-2)).*std_TD(i);

end


dose(1)=0.000001; 




% this removes the unimplanted density from the other points 
% this is needed since the other implanted points in TEM etc have only the
% irradiation induced defects, not the intrinsic ones 
% for i=2:length(Cv)
%     Cv(i)=Cv(i)-Cv(1);
% end
% %

% adjsuting since 0.0001 dpa gave higher value for TC , smaller value for
% Cv than unimp
for i=3:length(Cv)
    Cv(i)=Cv(i)-Cv(1);
end
%


% 
figure
errorbar(dose,Cv,std_Cv,'rx','LineWidth',1,'MarkerFaceColor','b','MarkerEdgeColor','b')

% hold on
% plot(dose,Cv,'--')
grid on
xlabel('Dose (dpa)','FontSize',16)
ylabel('No. of Frenkel Pairs (at. fr.) ','FontSize',16)
set(gcf,'color','w');
set(gca,'fontsize',16);
set(gca,'xscale','log')


% hold on 
% plot(dose,Cv,'bx')

% set(gca,'yscale','log')

xticks([1e-6 1e-4 1e-3 1e-2 1e-1 1e-0 10])
xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})

hold on
% plot(dose,dose,'r-')
% set(gca,'xscale','log')




% same one in number density 

% calc_num_dens=Cv*6.322e28  ; %  per m-3
% 
% figure 
% 
% plot(dose,calc_num_dens,'rx')
% % hold on
% % plot(dose,Cv,'--')
% grid on
% xlabel('Dose (dpa)','FontSize',16)
% ylabel('No. of Frenkel Pairs (Num density m-3) ','FontSize',16)
% set(gcf,'color','w');
% set(gca,'fontsize',16);
% set(gca,'xscale','log')
% 
% set(gca,'yscale','log')
% 
% xticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e-0 10])
% xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})

hold on



%% plotting with felix data  - run this after running the first section 
% % 
% load('Xiaoou_sand_TC_data.mat','dpa')
 
 % load the felix calculation dat a
 load('felix_calc_MDTEM.mat') % this has the dpa values as Xiaoo, and the estimated defect density visible and all, cirum and area
 % this was taken from the excel sheet which has the calculation -
%  felix_calc=felix_calc';
 
figure 
e=errorbar(dose,Cv,std_Cv,'rx','LineWidth',2,'MarkerFaceColor','b','MarkerEdgeColor','b');

% adjsuting fpr tehe negative error value in first two points - snding it
% to zero
e.YNegativeDelta(1:2)=e.YData(1:2)-1e-9;


% hold on
% plot(dose,Cv,'rx')

grid on
xlabel('Dose (dpa)','FontSize',16)
ylabel('Point defects (at. fr.) ','FontSize',16)
set(gcf,'color','w');
set(gca,'fontsize',16);
set(gca,'xscale','log')
set(gca,'yscale','log')
 hold on 

plot(dpa,felix_calc(:,1),'ks')
hold on
plot(dpa,felix_calc(:,2),'ms')
hold on
plot(dpa,felix_calc(:,3),'kd','MarkerFaceColor','k')
hold on
plot(dpa,felix_calc(:,4),'md','MarkerFaceColor','m')


xticks([1e-6 1e-4 1e-3 1e-2 1e-1 1e-0 10])
xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})

 yticks([1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2])
 yticklabels({'0','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}'})


hold on 
% adding the implanted ions data 

% load('implanted_ions_data.mat')
% 
% mean_ion_at_fr(1)=0;
% plot(dpa,mean_ion_at_fr,'b^')
% 
% legend({'Estimated from TC (TGS + KT)','In visible loops (TEM) (Area)','In visible loops (TEM) (Circum)','Estimated from MD + TEM (area)', 'Estimated from MD + TEM (Circum)','Implanted Ions' },'Location','northwest','FontSize',12)
% legend boxoff

% adding implanted ions data - actual values 
% this is the actual  data we got from the implantation people 
fluence=[0,2.7e10,8.13e10,2.42e11,8.03e11,2.55e12,4.61e12,8.2e12,1.42e13,2.54e13,8.11e13,2.53e14,8.1e14,2.53e15]; % ions/cm2

fluence2=fluence*1e4;
ion_dens=fluence2/2e-6; % 2e-6 is the implanted thickness

ion_dens_at_fr=ion_dens/6.322e28;

plot(dose,ion_dens_at_fr,'b^','MarkerFaceColor','b')

legend({'Estimated from TGS + KT','In TEM visible loops (Area)','In TEM visible loops (Circum)','Estimated from MD + TEM (Area)', 'Estimated from MD + TEM (Circum)','Implanted Ions' },'Location','northwest','FontSize',12)
% legend boxoff

xlim manual

xticks([1e-6 1e-4 1e-3 1e-2 1e-1 1e-0 10])
xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})

% ylim manual
 yticks([1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2])
 yticklabels({'0','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}'})

 
 