

%% Produces Figure 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Excess Bond Premia and compute quarterly averages 
[xff txt_xff xff_raw ]=xlsread('GZ.xls','Foglio1');
EBP_quarter=mean(reshape(xff(61:end,2),3,length(xff(61:end,1))/3));

%% Load financial index data 
[xff txt_xff xff_raw ]=xlsread('GZ.xls','financial index');
Fin_nw=xff(21:end);
Fin_nw=Fin_nw(1:length(EBP_quarter));


%% Load impulse responses associated to a recession with positive probability of a run and a run that actually happens one year after the initial shock
load anticipated_run_happens_at_time_4


%% Pick Model counterparts and let the steady state value correpond to the quarter before Bear Stern (10th entry)
Model_Spread=4*(Spreadk_prun-spread_val);

Model_Spread_EBP_07_on=Model_Spread;
Model_Spread_EBP_07_on(1)=EBP_quarter(10);
Model_Spread_EBP_07_on(2:end)=Model_Spread_EBP_07_on(1)+Model_Spread_EBP_07_on(2:end)*100;

Model_v=(theta*phi_prun.*n_prun-theta*phi_ss.*N_ss)/(theta*phi_ss.*N_ss);
Model_v_07_on=Model_v;
Model_v_07_on(1)=Fin_nw(10);
Model_v_07_on(2:end)=Model_v_07_on(1)*(1+Model_v_07_on(2:end));


% BS=13;
% LB=16;



figure
subplot(2,1,2)

plot(Fin_nw(10:end));
hold on
plot([NaN*ones(2,1) ; Model_v_07_on(1:length(EBP_quarter)-11)],'r')
subplot(2,1,1)
plot(EBP_quarter(10:end));
hold on
plot([NaN*ones(2,1) ; Model_Spread_EBP_07_on(1:length(EBP_quarter)-11)],'r')

quarters={'1st quarter','2nd quarter', '3rd quarter','4th quarter'};
years={'2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010'};
rep_quarters=repmat(quarters,1,length(years));
rep_years=reshape(repmat(years,4,1),1,length(rep_quarters));
rep_ticks=strcat(rep_quarters,rep_years);


set(gca,'XTickLabel',years);
