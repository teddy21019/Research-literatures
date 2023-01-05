clear all

%% load the workspaces containing impulse responses 
load unanticipated
load anticipated_run_does_not_happen
load anticipated_run_happens_at_time_3

%% extract and compute relevant variables

kh_prun=x_path(1:T);
q_prun=x_path(T+1:2*T);
d_prun=x_path(2*T+1:3*T);

n_prun=P.sigma*( (q_prun(2:end)+z(2:end)-mc_val*q_prun(1:end-1)).*(1-kh_prun(1:end-1)) - d_prun(1:end-1)   ) + P.w;

n_prun=[N_ss;n_prun];

ch_prun=y(2:end)-(1-P.sigma)*(n_prun(2:end)-P.w)./P.sigma - P.alpha/2*kh_prun(2:end).^2-mc_val*q_prun(1:end-1).*(1-kh_prun(1:end-1));
ch_prun=[ch_ss;ch_prun];
phi_prun=q_prun.*(1-kh_prun)./n_prun;
q_bar_prun=gamma_val*d_prun(1:end-2)./(1-kh_prun(1:end-2))-z(2:end-1)-mc_val*q_prun(1:end-2);

pai_check=P.pbar*(1-(q_run(3:end)+z(3:end)).*(1-kh_prun(2:end-1))./d_prun(2:end-1)).^P.delta_p;
pai=max(P.pbar*(1-(q_run(3:end)+z(3:end)).*(1-kh_prun(2:end-1))./d_prun(2:end-1)).^P.delta_p,0);
pai=[0;pai;pai(end)];


Rk_z=[1/beta+spread_val; (z(3:end) + q_z(3:end))./q_z(2:end-1) ;1/beta+spread_val ];
Rh_z=[1/beta; (z(3:end) + q_z(3:end))./(q_z(2:end-1)+alpha*kh_z(2:end-1)) ;1/beta+spread_val ];

Rk_zrun=[1/beta+spread_val; Rk_z(2:t_run); (z(t_run+2:end) + q_zrun(t_run+2:end))./q_zrun(t_run+1:end-1) ;1/beta+spread_val ];

Spreadk_z=Rk_z-R_z;
Spreadk_zrun=Rk_zrun-R_zrun;
Spreadk_prun=Rk_prun-R_prun;
Spreadd_prun=R_prun-R_rf_prun;
Run_z=1-(q_run_plot(3:end)+z(3:end)).*(1-kh_z(2:end-1))./(gamma_val*d_z(2:end-1));
Run_p=1-(q_run(3:end)+z(3:end)).*(1-kh_prun(2:end-1))./d_prun(2:end-1);


phi_after_run=q_path_after_run(3:end).*(1-kh_path_after_run(2:end-1))./((1+P.sigma)*P.w);
pai_dev=pai(2:end-3);
deviation_loss=(1-pai_dev).*(1-P.sigma+P.sigma*P.theta*phi_prun(3:end-2)).*( (z(3:end-2) + q_prun(3:end-2))./q_prun(2:end-3) - R_prun(2:end-3)  ).*phi_prun(2:end-3);
deviation_gain=pai_dev.*P.beta.*(1-P.sigma+P.sigma*P.theta*phi_after_run).*( (z(3:end-2) + q_run(3:end-2))./q_prun(2:end-3) );
daviation_value=deviation_gain-deviation_loss;

%% GENERATE PLOTS


%% Z SHOCK
SS=zeros(size(z(2:1+N)));
figure 
 
subplot(3,3,1)
plot((z(2:1+N)-Z_ss)/Z_ss  ,'LineWidth',2   )     %plot(oo_.irfs.z_eps_z(1:PERIODS)) 
hold on
plot(SS,'r');
title('z')
ylabel('%\Delta from ss') 


subplot(3,3,2)
plot((y_z(2:1+N)-y_ss)/y_ss,'LineWidth',2)  
hold on
plot(SS,'r');
title('y')
ylabel('%\Delta from ss') 


subplot(3,3,3)
plot((1-kh_z(2:1+N)-kb_ss)/kb_ss,'LineWidth',2) 
hold on
plot(SS,'r');
title('kb')
ylabel('%\Delta from ss') 


subplot(3,3,4)
plot((q_z(2:1+N)-Q_ss)/Q_ss,'LineWidth',2) 
hold on
plot(SS,'r');
title('Q')
ylabel('%\Delta from ss') 



subplot(3,3,5)
plot(4*(Spreadk_z(2:1+N)-spread_val),'LineWidth',2) 
hold on
plot(SS,'r');
title('ER^{b}-R')
ylabel('Ann. \Delta from ss') 


subplot(3,3,6)
plot((n_z(2:1+N)-N_ss)/N_ss,'LineWidth',2) 
hold on
plot(SS,'r');
title('n')
ylabel('%\Delta from ss') 



subplot(3,3,7)
plot(4*(R_z(2:1+N)-R_ss)/R_ss,'LineWidth',2) 
hold on
plot(SS,'r');
title('R')
xlabel('Quarters')
ylabel('Ann. %\Delta from ss') 


subplot(3,3,8)
plot((ch_z(2:1+N)-ch_ss)/ch_ss,'LineWidth',2) 
hold on
%legend('Recession','Location','SouthOutside')
plot(SS,'r');

title('ch')

xlabel('Quarters')
ylabel('%\Delta from ss') 
legend('z impulse','Location','EastOutside','Orientation', 'Horizontal')

subplot(3,3,9)
plot((cb_z(2:1+N)-cb_ss)/cb_ss,'LineWidth',2) 
hold on
plot(SS,'r');
title('cb')
xlabel('Quarters')
ylabel('%\Delta from ss') 














%% Unanticipated Run

figure 

subplot(3,3,1)
plot((y_z(2:1+N)-y_ss)/y_ss,'r--','LineWidth',2)  
hold on
plot((y(2:1+N)- P.alpha/2*kh_zrun(2:1+N).^2-y_ss)/y_ss,'LineWidth',2)
plot(SS,'r');
title('y')
ylabel('%\Delta from ss') 


subplot(3,3,2)
plot((1-kh_z(2:1+N)-kb_ss)/kb_ss,'r--','LineWidth',2) 
hold on
plot((1-kh_zrun(2:1+N)-kb_ss)/kb_ss,'LineWidth',2)
plot(SS,'r');
title('kb')
ylabel('%\Delta from ss') 


subplot(3,3,3)
plot((q_z(2:1+N)-Q_ss)/Q_ss,'r--','LineWidth',2) 
hold on
plot((q_zrun(2:1+N)-Q_ss)/Q_ss,'LineWidth',2) 
plot(SS,'r');
title('Q')
ylabel('%\Delta from ss') 

subplot(3,3,4)

plot((Run_z(1:N)),'r--','LineWidth',2) 
hold on
plot(SS(1:end),'r');
title('RUN')



subplot(3,3,5)

%plot((q_bar_norun(2:N+1)-(d_z(1)./(1-kh_z(1))-Z_ss))/(d_z(1)./(1-kh_z(1))-Z_ss),'LineWidth',2) 
plot( (q_run_plot(3:N+2)-q_run_plot(end))/q_run_plot(end),'r--','LineWidth',2) 
hold on
plot(SS,'r');
title('Q^{*}')

ylabel('%\Delta from ss') 



% 
subplot(3,3,6)
plot((phi_z(2:1+N)-phi_ss)/phi_ss,'r--','LineWidth',2) 
%hold on
plot((phi_zrun(2:1+N)-phi_ss)/phi_ss,'LineWidth',2)
hold on
plot(SS,'r');
title('\phi^{*}')
ylabel('%\Delta from ss') 



subplot(3,3,7)
plot((ch_z(2:1+N)-ch_ss)/ch_ss,'r--','LineWidth',2) 
hold on
plot((ch_zrun(2:1+N)-ch_ss)/ch_ss,'LineWidth',2) 
plot(SS,'r');
title('ch')
xlabel('Quarters')
ylabel('%\Delta from ss') 


subplot(3,3,8)
plot((cb_z(2:1+N)-cb_ss)/cb_ss,'r--','LineWidth',2) 
hold on
plot(((1-P.sigma)*(n_zrun(2:1+N)-P.w)./P.sigma-cb_ss)/cb_ss,'LineWidth',2) 
plot(SS,'r');

title('cb')
xlabel('Quarters')
ylabel('%\Delta from ss') 
legend('z impulse','unanticipated run','Location','EastOutside','Orientation', 'Horizontal')

subplot(3,3,9)
plot(4*(Spreadk_zrun(2:1+N)-spread_val),'LineWidth',2) 
hold on
plot(4*(Spreadk_z(2:1+N)-spread_val),'r--','LineWidth',2) 
hold on
plot(SS,'r');
title('ER^{b}-R')
xlabel('Quarters')
ylabel('Ann. \Delta from ss') 


%% Anticipated run

figure 
 
subplot(3,3,1)
plot((pai(2:1+N))  ,'LineWidth',2     ) 
hold on
plot(SS,'r');
title('p')
ylabel('\Delta from ss') 


subplot(3,3,2)
plot((y_z(2:1+N)-y_ss)/y_ss,'r--','LineWidth',2)  
hold on
plot((y(2:1+N)- P.alpha/2*kh_prun(2:1+N).^2-mc_val*q_prun(1:N).*(1-kh_prun(1:N))-y_ss)/y_ss,'LineWidth',2)
plot(SS,'r');
title('y')
ylabel('%\Delta from ss') 


subplot(3,3,3)
plot((1-kh_z(2:1+N)-kb_ss)/kb_ss,'r--','LineWidth',2) 
hold on
plot((1-kh_prun(2:1+N)-kb_ss)/kb_ss,'LineWidth',2)
plot(SS,'r');
title('kb')
ylabel('%\Delta from ss') 


subplot(3,3,4)
plot((q_z(2:1+N)-Q_ss)/Q_ss,'r--','LineWidth',2) 
hold on
plot((q_prun(2:1+N)-Q_ss)/Q_ss,'LineWidth',2) 
plot(SS,'r');
title('Q')
ylabel('%\Delta from ss') 

subplot(3,3,5)
plot((phi_z(2:1+N)-phi_ss)/phi_ss,'r--','LineWidth',2) 
hold on
plot((phi_prun(2:1+N)-phi_ss)/phi_ss,'LineWidth',2)
plot(SS,'r');
title('\phi')
ylabel('%\Delta from ss') 

 

subplot(3,3,6)
plot((n_z(2:1+N)-N_ss)/N_ss,'r--','LineWidth',2) 
hold on
plot((n_prun(2:1+N)-N_ss)/N_ss,'LineWidth',2) 
plot(SS,'r');
title('n')
ylabel('%\Delta from ss') 


subplot(3,3,7)
plot(4*(Spreadk_prun(2:1+N)-spread_val),'LineWidth',2) 
hold on
plot(4*(Spreadk_z(2:1+N)-spread_val),'r--','LineWidth',2)
plot(SS,'r');
title('ER^{b}-R^{d}')
xlabel('Quarters')
ylabel('Ann. \Delta from ss') 



subplot(3,3,8)
plot(4*Spreadd_prun(2:1+N),'LineWidth',2) 
hold on
plot(SS,'r--','LineWidth',2);
title('R^{d}-R^{free}')
xlabel('Quarters')
ylabel('Ann. \Delta from ss')
legend('p and z impulse','z impulse','Location','EastOutside','Orientation', 'Horizontal')

subplot(3,3,9)
plot(4*(R_z(2:1+N)-R_ss)/R_ss,'r--','LineWidth',2) 
hold on
plot(4*(R_rf_prun(2:1+N)-R_ss)/R_ss,'LineWidth',2) 
plot(SS,'r');
title('R^{free}')
%legend('z impulse','p and z impulse','Location','SouthOutside')
xlabel('Quarters')
ylabel('Ann. %\Delta from ss') 



