clc
clear all
close all
 
 
load('summary_with_refdata.mat')
% load('helsinki_summary_data_25_3.mat','sawv','std_sawv','TC','std_TC','dose')
load('helsinki_summary_data_4_4.mat','TC','std_TC','dose')
dose(1)=1e-5;

h3=figure;

% h2=figure;
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
% xticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e-0 10])
% xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})

%% plotting the unimplanted refs


plot(fukuda_dose(1),fukuda_TD(1),'ks','MarkerFaceColor','k','MarkerSize',10)

grid on 
xlabel('Dose (dpa)','FontSize',16)
ylabel('Thermal Diffusivity (m^{2}s^{-1})','FontSize',16)
set(gcf,'color','w');
set(gca,'fontsize',16);
 set(gca,'xscale','log')
 axis([1e-5 100 1e-5 8e-5])
xticks([ 1e-5 1e-4 1e-3 1e-2 1e-1 1e-0,10])
xticklabels({'Ref.','0.0001','0.001','0.01','0.1','1','10'})


hold on
plot(fukuda_dose(2),fukuda_TD(2),'kd','MarkerFaceColor','k','MarkerSize',10)

hold on
plot(fukuda_dose(3),fukuda_TD(3),'ko','MarkerFaceColor','k','MarkerSize',10)

hold on
plot(fukuda_dose(4),fukuda_TD(4),'kp','MarkerFaceColor','k','MarkerSize',10)


hold on
plot(fukuda_dose(5),fukuda_TD(5),'k^','MarkerFaceColor','k','MarkerSize',10)

hold on

% % adding felix TGS unimp value after review 
 plot(1e-5,6.75e-5,'k>','MarkerFaceColor','g','MarkerSize',10)
 


hold on


%% plotting the implanted ones


% plot(cui_dose,cui_TD,'bs')
errorbar(cui_dose,cui_TD,cui_std*2,'bs','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2)




 hold on
%  plot(deuch_dose(1),deuch_TD(1),'bx','MarkerFaceColor','b')
% hold on

% %% removing since his is cu ions
% errorbar(deuch_dose(2),deuch_TD(2),deuch_std*2,'ms','MarkerFaceColor','m','MarkerSize',10,'LineWidth',2)
% 
% hold on
plot(ESS_dose,ESS_TD,'mo','MarkerFaceColor','m','MarkerSize',10)
hold on
plot(peackock_dose,peackock_TD,'cd','MarkerFaceColor','c','MarkerSize',10)

hold on


errorbar(dose,TC,std_TC,'rx','LineWidth',2,'MarkerFaceColor','b','MarkerEdgeColor','b')

% hold on
% plot(dose,TC,'b.','MarkerFaceColor','b')

% h=legend('[8]',' [33] Poly. Cryst.','[33] Sin. Cryt.','[34]','[11]',',[12], Self-ion','[35], T_{irr}= 150 ^{o}C ','[36], T_{irr}= 200 ^{o}C ','This work');
% h.Location='northeast';
% h.FontSize=11;

% newOrder = [2 3 5 4 1 6 7 8 9];
% [~,~,plot_h,text_strings] = legend;
% % legend(plot_h(newOrder),text_strings{newOrder})
% 
% set(h,'String',text_strings(newOrder),'Data',plot_h(newOrder))


%      savefig(h3,'helsinki_summary_with_refs_25_3_TC.fig')
% saveas(h3,'helsinki_summary_with_refs_25_3_TC.png')

% adding the as received point 
hold on 
load('tungsten_plansee_as_received_1_17delay_analysis.mat','map_saw','std_saw','map_diffuse','std_diffuse','map_vel','std_vel','ph','pv','phh')

TC_as=mean(mean(map_diffuse));
std_TC_as=mean(mean(std_diffuse));
% 


hold on 
errorbar(dose(1),TC_as,std_TC_as,'gx','LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','r')

% h=legend('Fujitsuka, Pure W [15]','Fukuda, Pure W [41]','Fukuda, Sin. Cryst. [41]','Touloukian, Pure W [42]','Tanabe, Pure W [18]','Hofmann, Pure W; TGS [16]','Cui, Self-ion [19]','Habainy, ''n''; T_{irr}= 150 ^{o}C [43]','Linke, ''n''; T_{irr}= 200 ^{o}C [44]','This work (Prior Annealed)', 'This work (As Received) ');

h=legend('Pure W, Fujitsuka [15]','Pure W, Fukuda [41]','Pure W, Sin. Cryst., [41]','Pure W, Touloukian [42]','Pure W, Tanabe [18]','Pure W, TGS, Hofmann [16]','Self-ion impl., Cui [19]','''n'' T_{irr}= 150 ^{o}C, Habainy [43]','''n'' T_{irr}= 200 ^{o}C, Linke [44]','This work (Prior Annealed)', 'This work (As Received) ');

h.Location='northeastoutside';
legend('boxoff')

% h.Location='northeast';
h.FontSize=11;

axis([0.00001 10 2e-5 7.5e-5])


