% This script plots ground motion predictions using empirical Gorund Motion 
% Prediction Models developed for Iran. 
% The predictions correspond to horizontal, vertical and vertical-to-horizontal components of ground motions. 
% 
% Author: Atefe Darzi, Last Edit: 2019
% 
% References:
% 1- Darzi A, Zolfaghari MR, Cauzzi C, Fäh D (2019). An Empirical Ground Motion Model 
% for Horizontal PGV, PGA and 5%-Damped Elastic Response Spectra (0.01-10 s) in Iran. 
% Bulletin of the Seismological Society of America. 10.1785/0120180196 

% 2- Zolfaghari MR, Darzi A (2019). Ground-Motion Models for Predicting Vertical 
% Components of PGA, PGV and 5%-damped Spectral Acceleration (0.01-10 s) in Iran. 
% Bulletin of Earthquake Engineering, 10.1007/s10518-019-00623-2  

% 3- Zolfaghari MR, Darzi A (2019). A Prediction Model for Vertical-to-Horizontal 
% Ratios of PGA, PGV and 5%-damped Response Spectra (0.01-10 s) for Iran. 
% Journal of Seismology. 10.1007/s10950-019-09836-z 
    
clear all;close all;clc

%% predict V/H 
Tnga = load('Periods-nga.txt'); T=Tnga(1:60)';%SA(0.01-...-10)
per = [-2 -1 0 T]; % PGD,PGV,PGA, SA (1*63)

% Select an approriate period-dependent cofficient (pdc) file based on
% your choice of distance type (Rjb, Repi, Rhyp, Rrup) and functional form (fun1, fun3)
pdc_H2V=load('IrV2H-Rjb-fun3.txt')'; % Rjb-fun3
%Coef_fun3={'fRV';'fSS';'C1';'m1';'h';'r1';'S(II)';'S(III-IV)';'phi';'tau';'sigma'};

%%%% plot GM predictions(V/H) vs distance
gcf1=figure(1);

subplot(1,2,1); ord=find(per==0)
RR=logspace(log10(1),log10(250),50); 
[Sa1,sig]=V2H_GMM_ZD19('Rjb','fun3',[5 6 7],RR,ord,pdc_H2V,'II','R',0); 
[Sa11,sig]=V2H_GMM_ZD19('Rjb','fun3',[5 6 7],RR,ord,pdc_H2V,'II','R',1);
[Sa12,sig]=V2H_GMM_ZD19('Rjb','fun3',[5 6 7],RR,ord,pdc_H2V,'II','R',-1);
r1=loglog(RR,Sa1,'k','linewidth',1.2); hold on; grid on;
r2=loglog(RR,Sa11,':k','linewidth',1); loglog(RR,Sa12,':k','linewidth',1);

ylabel('\bf \fontsize{20} V/H Ratio');xlabel('\bf \fontsize{20} Distance (km)'); 
ax=gca;set(ax,'FontSize',22); xlim([10 200]);axis square;
set(ax,'XTick',[10 100],'XTickLabel',{'10','100'});
set(ax,'YTick',[0.1 1 10],'YTickLabel',{'0.1','1','10'});
ylim([0.05 10]); title('\bf \fontsize{22} PGA')
legend({'Mw 5','Mw 6','Mw 7'},'Location','southeast');
%
subplot(1,2,2); ord=find(per==1)
RR=logspace(log10(1),log10(250),50); 
[Sa1,sig]=V2H_GMM_ZD19('Rjb','fun3',[5 6 7],RR,ord,pdc_H2V,'II','R',0); 
[Sa11,sig]=V2H_GMM_ZD19('Rjb','fun3',[5 6 7],RR,ord,pdc_H2V,'II','R',1);
[Sa12,sig]=V2H_GMM_ZD19('Rjb','fun3',[5 6 7],RR,ord,pdc_H2V,'II','R',-1);
r1=loglog(RR,Sa1,'k','linewidth',1.2); hold on; grid on;
r2=loglog(RR,Sa11,':k','linewidth',1); loglog(RR,Sa12,':k','linewidth',1);

ylabel('\bf \fontsize{20} V/H Ratio');xlabel('\bf \fontsize{20} Distance (km)'); 
ax=gca;set(ax,'FontSize',22); xlim([10 200]);axis square;
set(ax,'XTick',[10 100],'XTickLabel',{'10','100'});
set(ax,'YTick',[0.1 1 10],'YTickLabel',{'0.1','1','10'});
ylim([0.05 10]); title('\bf \fontsize{22} T=1 s')

saveFigure(gcf1,'V2H-distance.jpg'); 
   
%%%% plot GM predictions(V/H) vs magnitude
gcf2=figure(2);
subplot(1,2,1); ord=find(per==0)
Mag=4.5:0.1:7.8; 
[Sa1,sig]=V2H_GMM_ZD19('Rjb','fun3',Mag,[30 60 90],ord,pdc_H2V,'II','R',0); 
[Sa11,sig]=V2H_GMM_ZD19('Rjb','fun3',Mag,[30 60 90],ord,pdc_H2V,'II','R',1);
[Sa12,sig]=V2H_GMM_ZD19('Rjb','fun3',Mag,[30 60 90],ord,pdc_H2V,'II','R',-1);
r1=semilogy(Mag,Sa1,'k','linewidth',1.5); hold on; grid on;
r2=semilogy(Mag,Sa11,':k','linewidth',1);semilogy(Mag,Sa12,':k','linewidth',1);

ylabel('\bf \fontsize{20} V/H Ratio'); xlabel('\bf \fontsize{20} Mw');
ax=gca;set(ax,'FontSize',22); xlim([5 7.5]);axis square;
set(ax,'XTick',[5 6 7],'XTickLabel',{'5','6','7'});
set(ax,'YTick',[0.1 1 10],'YTickLabel',{'0.1','1','10'});
ylim([0.05 10]); title('\bf \fontsize{22} PGA')
legend({'R 30 km','R 60 km','R 90 km'},'Location','southeast');
%
subplot(1,2,2); ord=find(per==1)
Mag=4.5:0.1:7.8; 
[Sa1,sig]=V2H_GMM_ZD19('Rjb','fun3',Mag,[30 60 90],ord,pdc_H2V,'II','R',0); 
[Sa11,sig]=V2H_GMM_ZD19('Rjb','fun3',Mag,[30 60 90],ord,pdc_H2V,'II','R',1);
[Sa12,sig]=V2H_GMM_ZD19('Rjb','fun3',Mag,[30 60 90],ord,pdc_H2V,'II','R',-1);
r1=semilogy(Mag,Sa1,'k','linewidth',1.5); hold on; grid on;
r2=semilogy(Mag,Sa11,':k','linewidth',1);semilogy(Mag,Sa12,':k','linewidth',1);

ylabel('\bf \fontsize{20} V/H Ratio'); xlabel('\bf \fontsize{20} Mw');
ax=gca;set(ax,'FontSize',22); xlim([5 7.5]);axis square;
set(ax,'XTick',[5 6 7],'XTickLabel',{'5','6','7'});
set(ax,'YTick',[0.1 1 10],'YTickLabel',{'0.1','1','10'});
ylim([0.05 10]); title('\bf \fontsize{22}  T=1 s')
saveFigure(gcf2,'V2H-magnitude.jpg'); 
 
%% predict Horizontal accelerations 
Tnga = load('Periods-nga.txt'); % T=0.01-20 sec
T1=Tnga(1:60)'; per = [-2 -1 0 T1]; % PGD,PGV,PGA, SA

% Select an approriate period-dependent cofficient (pdc) file based on
% your choice of distance type (Rjb, Repi, Rhyp, Rrup) 
pdc_H=load('IrH-jb.txt')'; % Rjb

gcf3= figure(3)
R= [5 80];% mech=SS Vs30=800
Mw=5;
ord=1:length(per);
[Sa1,sigmat]=H_GMM_Darzi19(Mw,R,ord,pdc_H,'I','S',0);
[Sa11,sigmat1]=H_GMM_Darzi19(Mw,R,ord,pdc_H,'I','S',1);
[Sa12,sigmat2]=H_GMM_Darzi19(Mw,R,ord,pdc_H,'I','S',-1);
l1=loglog(per,Sa1,'k','linewidth',3);  hold on; grid on
l2=loglog(per,Sa11,'--k','linewidth',2.2); loglog(per,Sa12,'--k','linewidth',2.2);

xlabel('\bf \fontsize{24} period (s)'); ylabel('\bf \fontsize{24} Horizontal SA (cm/s^{2})');
xlim([0.01 10]);title('\bf\fontsize{24} Mw 5');axis square
ax=gca;set(ax,'FontSize',22); set(ax,'XTick',[0.1 1 10],'XTickLabel',{'0.1','1','10'}); 
ylim([1 1500]); set(ax,'YTick',[10 100 1000],'YTickLabel',{'10','100','1000'}); 

saveFigure(gcf3,'H-PSA.jpg'); 


%% predict Vertical accelerations 
Tnga = load('Periods-nga.txt'); % T=0.01-20 sec
T1=Tnga(1:60)'; per = [-2 -1 0 T1]; % PGD,PGV,PGA, SA

% Select an approriate period-dependent cofficient (pdc) file based on
% your choice of distance type (Rjb, Repi, Rhyp, Rrup) 
pdc_V=load('IrV-jb.txt')'; % Rjb

gcf4= figure(4)
R= [5 80];% mech=SS Vs30=800
Mw1=5;
ord=1:length(per); 
[Sa1,sig]=V_GMM_ZD19(Mw1, R, ord, pdc_V,'I','S',0);
l1=loglog(per,Sa1(:,1)./981,'-k','linewidth',2.5); hold on; grid on;
l2=loglog(per,Sa1(:,2)./981,'-k','linewidth',2.5);

ylabel('\bf \fontsize{15} Vertical SA (g)');%xlabel('\bf \fontsize{15} period (s)');
ax=gca;set(ax,'FontSize',15); xlim([0.01 10]); ylim([0.0001 2]); %axis square;
set(ax,'XTick',[ 0.1 1 10 ],'XTickLabel',{'0.1','1','10'});  
set(ax,'YTick',[0.001 0.01 0.1 1],'YTickLabel',{'0.001','0.01','0.1','1'});
text(0.1,0.003,'\bf \fontsize{13} RJB=80 km');text(.5,.05,'\bf \fontsize{11} RJB=5 km');
title('\bf \fontsize{15} Mw=5, site class I'); axis square

saveFigure(gcf4,'V-PSA.jpg'); 

