clc
clear all
close all

%% This is the analysis of the TC vs DOSE data 
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
te=21.4e-15; % felix calculation - get a ref 
te=22.918e-15; % same calculation with 2 terms only

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
% 
% 
% Cv=zeros(1,length(TD));
% 
% std_Cv= zeros(1,length(TD)); % using simple error theory 
% 
% for i=1:length(TD)
% Cv(i)=((((te*Ce*(vf)^2)/(3*rho*Cp.*TD(i)))-1))/(te*(sigmav+sigmai)-2);
% 
% std_Cv(i)=(((((te*Ce*(vf)^2)/(3*rho*Cp.*(TD(i)).^2))))/(te*(sigmav+sigmai)-2)).*std_TD(i);
% 
% end


dose(1)=0.000001; 




% this removes the unimplanted density from the other points 
% this is needed since the other implanted points in TEM etc have only the
% irradiation induced defects, not the intrinsic ones 
% for i=2:length(Cv)
%     Cv(i)=Cv(i)-Cv(1);
% end
% %

% adjsuting since 0.0001 dpa gave higher value for TC , smaller value for
% % Cv than unimp
% for i=3:length(Cv)
%     Cv(i)=Cv(i)-Cv(1);
% end
% %



  load('felix_calc_MDTEM.mat') % this has the dpa values as Xiaoo, and the estimated defect density visible and all, cirum and area


  %%  TC back calculations from defects 
 
 % taking the defect density from the required calculated value - here
 % felix's calculation ,can use my calcualtion too
%  defect_dens=felix_calc(:,1);

 % giving zero as unimp dpa and hten the values from xiaooou 

  dpa2=[1e-6 ;dpa];  % for plotting purposes wioth ref value 
  dpa2=dpa2';

  
  % getting the different defect densities 
  
  defect_dens_vis_area=[0;felix_calc(:,1)]; % adding unimplanted value and rest is from the calculated 
  defect_dens_vis_circum=[0;felix_calc(:,2)];
    defect_dens_total_area=[0;felix_calc(:,3)];
  defect_dens_total_circum=[0;felix_calc(:,4)];
    
 % just the expression from scirep 
 % for each defect density we calculate what the back calculated TC would
 % be 
 
for i=1:length(defect_dens_vis_area)
 TC_vis_area(i)= ((te*Ce*(vf)^2)/(3*rho*Cp)).*(1./(defect_dens_vis_area(i)*(te*(sigmav+sigmai)-2) +1));
 
  TC_vis_circum(i)= ((te*Ce*(vf)^2)/(3*rho*Cp)).*(1./(defect_dens_vis_circum(i)*(te*(sigmav+sigmai)-2) +1));
  
   TC_total_area(i)= ((te*Ce*(vf)^2)/(3*rho*Cp)).*(1./(defect_dens_total_area(i)*(te*(sigmav+sigmai)-2) +1));
 
   
   TC_total_circum(i)= ((te*Ce*(vf)^2)/(3*rho*Cp)).*(1./(defect_dens_total_circum(i)*(te*(sigmav+sigmai)-2) +1));
 
end


figure 

errorbar(dose,TD,std_TD,'rx','LineWidth',2,'MarkerFaceColor','b','MarkerEdgeColor','b')

xlabel('Dose (dpa)','FontSize',16)
ylabel('Thermal Diffusivity (m^{2}s^{-1}) ','FontSize',16)
grid on
set(gcf,'color','w');
set(gca,'fontsize',16);
set(gca,'xscale','log')

hold on 
plot(dpa2,TC_vis_area,'kd','MarkerFaceColor','k')
hold on
plot(dpa2,TC_vis_circum,'md','MarkerFaceColor','m')
hold on 
plot(dpa2,TC_total_area,'k^','MarkerFaceColor','k')
hold on
plot(dpa2,TC_total_circum,'m^','MarkerFaceColor','m')
xticks([1e-6 1e-4 1e-3 1e-2 1e-1 1e-0 10])
xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})

axis([1e-6 10 1e-5 7.5e-5])

% adding the as received point 
hold on 
load('tungsten_plansee_as_received_1_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','map_vel','std_vel','ph','pv','phh')

TC_as=mean(mean(map_diffuse));
std_TC_as=mean(mean(std_diffuse));
sawv_as=mean(mean(map_vel));
std_sawv_as=mean(mean(std_vel));


hold on 
errorbar(dose(1),TC_as,std_TC_as,'gx','LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','r')


hold off
legend({'TGS measurement','TEM visible defects (Area)','TEM visible defects (Circum)','Total defects (Area)', 'Total defects (Circum)', 'As Rolled' },'Location','northwest','FontSize',12)



 