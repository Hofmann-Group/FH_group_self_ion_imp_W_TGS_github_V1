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

%% adding line if there was no annhilation of defects 
dose2=logspace(-6,1,100);

% this is in dp - now in number density 

% dose2_num_dens=dose2.*6.322e28;


% plot(dose2,dose2_num_dens,'r')
% axis([0 10 0 3e26])

plot(dose2,dose2,'m')
axis([0 10 0 3e-3])


hold on

%% fitting the curve from 0.00032 to 0.01 dpa before the saturation  - doing a Y=Ax type for 0.00032 upto 0.032 dpa

dose3=dose(2:5);
Cv3=Cv(2:5);
 Cv3=Cv(2:5);
%   cftool

% % % this was from a 2 term exponential
% % 
% %        a =   8.567e+25 ;
% %        b =       3.334 ;
% %        c =  -7.696e+25 ;
% %        d =      -135.6;
% %  fit_exp_no_sat = a*exp(b.*dose2) + c*exp(d.*dose2);
% %  
% % plot(dose2, fit_exp_no_sat,'k')
% 
% % cftool
    

%        a =   2.666e+26 ;
%        b =      0.3218  ;
% 
%  fit_power_1_no_satur = a*dose2.^b;
% plot(dose2,fit_power_1_no_satur,'k')

% from the Y=Ax fitting 

% plot(dose2, 0.08416*dose2,'k-')
% plot(dose2, 0.1375*dose2,'k-')

plot(dose2, 0.2274*dose2,'k-') % upto 0.0032

l=legend ('Obtained Density','Without Recombination','Without Saturation');
l.FontSize=12;
legend('boxoff')
% % y=ax + c fitting 
%        p1 =     0.06288 ;
%        p2 =   0.0002848 ;
% 
%     no_sat = p1*dose2 + p2;

% plot(dose2, no_sat,'m-')

