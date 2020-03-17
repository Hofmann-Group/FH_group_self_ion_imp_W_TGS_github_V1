clc
clear all
close all

% this is the summary plot - for diffusiviyt / saw vs dose 

% import the map data and average , with the errors 

n=13;
dose=zeros(1,n);
TC=zeros(1,n);
std_TC=zeros(1,n);
sawf=zeros(1,n);
std_sawf=zeros(1,n);




i=1;
load('Processed Data/helsinki_unimp_scan2_17delay/helsinki_unimp_scan2_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.000001;
TC(i)=mean(mean(map_diffuse));
std_TC(i)=mean(mean(std_diffuse));
sawf(i)=mean(mean(map_saw));
std_sawf(i)=mean(mean(std_saw));




i=2;
load('Processed Data/helsinki_imp_0_0001_scan2/helsinki_imp_0_0001_scan2_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.0001;
TC(i)=mean(mean(map_diffuse));
std_TC(i)=mean(mean(std_diffuse));
sawf(i)=mean(mean(map_saw));
std_sawf(i)=mean(mean(std_saw));


i=3;
load('Processed Data/helsinki_imp_0_00032_dpa_scan4_17delay/helsinki_imp_0_00032_dpa_scan4_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.00032;
TC(i)=mean(mean(map_diffuse));
std_TC(i)=mean(mean(std_diffuse));
sawf(i)=mean(mean(map_saw));
std_sawf(i)=mean(mean(std_saw));



i=4;
load('Processed Data/helsinki_imp_0_001dpa_scan2_17delay/helsinki_imp_0_001dpa_scan2_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.001;
TC(i)=mean(mean(map_diffuse));
std_TC(i)=mean(mean(std_diffuse));
sawf(i)=mean(mean(map_saw));
std_sawf(i)=mean(mean(std_saw));

clear map_saw std_saw map_diffuse std_diffuse ph pv 


i=5;
load('Processed Data/helsinki_imp_0_0032dpa_scan1_17delay/helsinki_imp_0_0032dpa_scan1_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.0032;
TC(i)=mean(mean(map_diffuse));
std_TC(i)=mean(mean(std_diffuse));
sawf(i)=mean(mean(map_saw));
std_sawf(i)=mean(mean(std_saw));

clear map_saw std_saw map_diffuse std_diffuse ph pv 





i=6;
load('Processed Data/helsinki_imp_0_01dpa_scan2_17delay/helsinki_imp_0_01dpa_scan2_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.01;
TC(i)=mean(mean(map_diffuse));
std_TC(i)=mean(mean(std_diffuse));
sawf(i)=mean(mean(map_saw));
std_sawf(i)=mean(mean(std_saw));

clear map_saw std_saw map_diffuse std_diffuse ph pv 


i=7;
load('Processed Data/helsinki_imp_0_018dpa_scan2_17delay/helsinki_imp_0_018dpa_scan2_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.018;
TC(i)=mean(mean(map_diffuse(:,7:end)));
std_TC(i)=mean(mean(std_diffuse(:,7:end)));
sawf(i)=mean(mean(map_saw(:,7:end)));
std_sawf(i)=mean(mean(std_saw(:,7:end)));

clear map_saw std_saw map_diffuse std_diffuse ph pv 


i=8;
load('Processed Data/helsinki_imp_0_032dpa_scan2_17delay/helsinki_imp_0_032dpa_scan2_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.032;
TC(i)=mean(mean(map_diffuse(2:end,:)));
std_TC(i)=mean(mean(std_diffuse(2:end,:)));
sawf(i)=mean(mean(map_saw(2:end,:)));
std_sawf(i)=mean(mean(std_saw(2:end,:)));

clear map_saw std_saw map_diffuse std_diffuse ph pv 

i=9;
load('Processed Data/helsinki_imp_0_056dpa_scan2_17delay/helsinki_imp_0_056dpa_scan2_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.056;
TC(i)=mean(mean(map_diffuse));
std_TC(i)=mean(mean(std_diffuse));
sawf(i)=mean(mean(map_saw));
std_sawf(i)=mean(mean(std_saw));


i=10;
load('Processed Data/helsinki_imp_new0_1_dpa_scan2/helsinki_imp_new0_1_dpa_scan2_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.1;
TC(i)=mean(mean(map_diffuse));
std_TC(i)=mean(mean(std_diffuse));
sawf(i)=mean(mean(map_saw));
std_sawf(i)=mean(mean(std_saw));


 clear map_saw std_saw map_diffuse std_diffuse ph pv 
 
 i=11;
load('Processed Data/helsinki_imp_0_32dpa_scan4_17delay/helsinki_imp_0_32dpa_scan4_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=0.32;
TC(i)=mean(mean(map_diffuse(1:end-6,1:end-6)));
std_TC(i)=mean(mean(std_diffuse(1:end-6,1:end-6)));
sawf(i)=mean(mean(map_saw(1:end-6,1:end-6)));
std_sawf(i)=mean(mean(std_saw(1:end-6,1:end-6)));

 clear map_saw std_saw map_diffuse std_diffuse ph pv 

 
 
  i=12;
load('Processed Data/helsinki_imp_1dpa_scan4_17delay/helsinki_imp_1dpa_scan4_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=1;
TC(i)=mean(mean(map_diffuse(:,3:end-3)));
std_TC(i)=mean(mean(std_diffuse(:,3:end-3)));
sawf(i)=mean(mean(map_saw(:,3:end-3)));
std_sawf(i)=mean(mean(std_saw(:,3:end-3)));

 clear map_saw std_saw map_diffuse std_diffuse ph pv 
 
  i=13;
load('Processed Data/helsinki_imp_3_2_dpa_scan1_17delay/helsinki_imp_3_2_dpa_scan1_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=3.2;
TC(i)=mean(mean(map_diffuse(7,:)));
std_TC(i)=mean(mean(std_diffuse(7,:)));
sawf(i)=mean(mean(map_saw(7,:)));
std_sawf(i)=mean(mean(std_saw(7,:)));
 
 clear map_saw std_saw map_diffuse std_diffuse ph pv 
 
  i=14;
load('Processed Data/helsinki_imp_10dpa_scan3_17delay/helsinki_imp_10dpa_scan3_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','ph','pv')
dose(i)=10;
TC(i)=mean(mean(map_diffuse(6:7)));
std_TC(i)=mean(mean(std_diffuse));
sawf(i)=mean(mean(map_saw));
std_sawf(i)=mean(mean(std_saw));


h1=figure;
errorbar(dose,TC,std_TC)

hold on
plot(dose,TC, 'r','LineWidth',1)
grid on 
xlabel('Dose (dpa)','FontSize',16)
ylabel('Thermal Diffusivity (m^{2}s^{-1})','FontSize',16)
set(gcf,'color','w');
set(gca,'fontsize',16);
 set(gca,'xscale','log')
xticks([1e-6 1e-4 1e-3 1e-2 1e-1 1e-0,10])
xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})

% h2=figure;
% 
% errorbar(dose,sawf,std_sawf)
% 
% hold on
% plot(dose,sawf, 'r','LineWidth',1)
% grid on 
% xlabel('Dose (dpa)','FontSize',16)
% ylabel('Peak SAW Frequency (Hz)','FontSize',16)
% set(gcf,'color','w');
% set(gca,'fontsize',16);
%  set(gca,'xscale','log')
% 
% xticks([1e-6 1e-4 1e-3 1e-2 1e-1 1e-0 10])
% xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})



% sawv=sawf*2.7528e-06;
% std_sawv=std_sawf*2.7528e-06;
% 
% h3=figure;
% 
% errorbar(dose,sawv,std_sawv)
% 
% hold on
% plot(dose,sawv, 'r','LineWidth',1)
% grid on 
% xlabel('Dose (dpa)','FontSize',16)
% ylabel('Peak SAW Velocity (ms^{-1})','FontSize',16)
% set(gcf,'color','w');
% set(gca,'fontsize',16);
%  set(gca,'xscale','log')
% 
% xticks([1e-6 1e-4 1e-3 1e-2 1e-1 1e-0 10])
% xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})






saving=0;
%% saving


if saving==1
  %cd Analysis
 save('helsinki_summary_data_4_4.mat','sawv','std_sawv','TC','std_TC','dose','sawf','std_sawf')
end

if saving==1

    savefig(h1,'helsinki_summary_4_4_TC.fig')
saveas(h1,'helsinki_summary_4_4_TC.png')

    savefig(h3,'helsinki_summary_4_4_SAWV.fig')
saveas(h3,'helsinki_summary_4_4_SAWV.png')

end

%% adding the as received point 
load('Processed Data/tungsten_plansee_as_received_1_17delay/tungsten_plansee_as_received_1_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','map_vel','std_vel','ph','pv','phh')

TC_as=mean(mean(map_diffuse));
std_TC_as=mean(mean(std_diffuse));
sawv_as=mean(mean(map_vel));
std_sawv_as=mean(mean(std_vel));

figure 
errorbar(dose,TC,std_TC)

hold on
plot(dose,TC, 'r','LineWidth',1)
grid on 
xlabel('Dose (dpa)','FontSize',16)
ylabel('Thermal Diffusivity (m^{2}s^{-1})','FontSize',16)
set(gcf,'color','w');
set(gca,'fontsize',16);
 set(gca,'xscale','log')
xticks([1e-6 1e-4 1e-3 1e-2 1e-1 1e-0,10])
xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})

hold on 
errorbar(dose(1),TC_as,std_TC_as)

hold on
plot(dose(1),TC_as, 'r','LineWidth',1)

