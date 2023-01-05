clear all
global P y z N_ss pai ch_run q_run T kh_ss d_ss mc_val ch_ss phi_val gamma_val
tic;


%% calibration 
pbar=1;
delta_p=1;
alpha_val=.00797;
phi_val=10;
spread_val=.0025;
mc_val=0;
sigma_val=.95;
inno_z=.05;
inno_p=.01;
rho_z_val=.95;
delta_rho=0;
gamma_val=1; %0.985; 


%% parameters targets and steadyt state
rho_z=rho_z_val;
beta =.99;
totmc=.1*.096;
alpha =alpha_val;
k=1;
wh=.045;
sigma= sigma_val ;
Q_ss=1;
phi_ss=phi_val;    
s=spread_val;
mc=mc_val; % disregard this

R_ss=1/beta;
Rk_ss=R_ss+s;
Rh_ss=1/beta;
Z_ss=Q_ss*(Rk_ss-1);
kh_ss=(Z_ss+Q_ss-Rh_ss*Q_ss)/(alpha*Rh_ss);
kb_ss=k-kh_ss; 
N_ss=Q_ss*kb_ss/phi_ss;
w=N_ss*( 1- sigma*((Rk_ss-mc-R_ss)*phi_ss+R_ss));
y_gross_ss=Z_ss+wh+w;
y_ss=y_gross_ss-(alpha/2)*(kh_ss^(2))-mc*kb_ss; 
cb_ss=(1-sigma)*((Rk_ss-mc-R_ss)*phi_ss+R_ss)*N_ss;
ch_ss=y_ss-cb_ss;
theta=( (1-sigma)*beta*(R_ss+phi_ss*(s-mc)) )/(phi_ss*(1-sigma*beta*(R_ss+phi_ss*(s-mc))) );

P=struct('theta',theta,'sigma',sigma,'alpha',alpha,'beta',beta,'Z_ss',Z_ss,'w',w,'wh',wh,'k',k,'rho_z',rho_z,'delta_p',delta_p,'pbar',pbar);


actual_run=3;
antirun=1;
Periods=120;
Periods_onlyhh=1; %Periods after ther run in which banks don't operate
T=Periods+2;
N=40; %Plot periods





%% Anticipated Run


%% Steady State: the long run value of a system in which a run is possible but never happens
%%set the path for y z and ch_run to their ss value
y=y_gross_ss*ones(T,1);
z=Z_ss*ones(T,1);
ch_run=(P.Z_ss+P.wh-P.alpha/2)*ones(T,1);

%%loop over q_run_ss to find the price of capital when a run happens in
%%steady state. 
q_run_ss=.9;
x_path_ss_after_run=[kh_ss*ones(T,1) ; ones(T,1) ; (1-kh_ss-N_ss)/beta*ones(T,1)];
for i=1:100
    %%compute the mopdified steady state a
    [P.kh_ss_prun P.d_ss_prun P.q_ss_prun]=prun_ss_endo( [kh_ss ; 1 ; (1-kh_ss-N_ss)/beta] , q_run_ss);
    q_run=q_run_ss*ones(T,1);
    
    [x_path_ss_after_run,residual_value,convergence_myrun]=fsolve(@residual_endorun_afterrun,[kh_ss*ones(T,1) ; ones(T,1) ; (1-kh_ss-N_ss)/beta*ones(T,1)]);
    
        
    kh_ss_after_run=x_path_ss_after_run(2);
    q_ss_after_run=x_path_ss_after_run(T+2);
    ch_ss_after_run=y(2)-(1-P.sigma)*P.w - P.alpha/2*kh_ss_after_run.^2;
    
    q_run_ss_new=-P.alpha + P.beta*(  ch_run(1)/ch_ss_after_run*( z(2) + q_ss_after_run  )  );
    
    if abs(q_run_ss_new-q_run_ss)<.00000001
        break
    else
        q_run_ss=0*q_run_ss+1*q_run_ss_new;
    end
    i
    
    
    
    
end
    
    



%% Z path
%%create the path for z from t=1 to t=T
z_path=exp(-inno_z*rho_z.^(0:T-3))*P.Z_ss;

%%initialize output
q_path_after_run=ones(Periods,1);
kh_path_after_run=ones(Periods,1);
q_after_run=q_ss_after_run;
ch_after_run=ch_ss_after_run;

%%Inductively create a path for the economy after a run happened at time t, 
%%x_path_ss_after_run, and a path of run values for q and ch from t onwards
%%(See Section 7.3 of Appendix)
 
for i=1:T-2
    t=T-2-i+1;
    z=[P.Z_ss;z_path(t);z(2:end-1)];
    y=[y_gross_ss;z_path(t)+z_path(t)./P.Z_ss*P.wh+w;y(2:end-1)];
    
    
    ch_run_new=z_path(t)+z_path(t)./P.Z_ss*P.wh-P.alpha/2;
    q_run_new=-P.alpha + P.beta*(  ch_run_new/ch_after_run*( z(3) + q_after_run  )  );
    
    if t==actual_run
        q_actual_run=q_run_new;
        ch_actual_run=ch_run_new;
        x_path_run=x_path_ss_after_run;
    end
    
    
    
    [x_path_ss_after_run,residual_value,convergence_myrun]=fsolve(@residual_endorun_afterrun,x_path_ss_after_run);
    
    
    
    kh_after_run=x_path_ss_after_run(2);
    q_after_run=x_path_ss_after_run(T+2);
    ch_after_run=y(2)-(1-P.sigma)*P.w - P.alpha/2*kh_after_run.^2;
    q_path_after_run(t)=q_after_run;
    kh_path_after_run(t)=kh_after_run;
    
    ch_run=[ch_run(1:2);ch_run_new;ch_run(3:end-1)];
    q_run=[q_run(1:2);q_run_new;q_run(3:end-1)];
    
    
    
    
    
end

q_run=[q_run(2:end);q_run_ss];
ch_run=[ch_run(2:end);(P.Z_ss+P.wh-P.alpha/2)];


[x_path,residual_value,convergence_myrun]=fsolve(@residual_endorun,x_path_ss_after_run);
if ~convergence_myrun
   load anticipated_run_root
   [x_path,residual_value,convergence_myrun]=fsolve(@residual_endorun,x_path_original);
      
end
x_path_original=x_path;


if actual_run<120
%%run


x_path(T+actual_run+1)=q_actual_run;



%%after run
x_path(actual_run+1:T)=x_path_run(1:T-actual_run);
x_path(T+actual_run+2:2*T)=x_path_run(T+2:2*T-actual_run);
x_path(2*T+actual_run+1:3*T)=x_path_run(2*T+1:3*T-actual_run);

end








kh_prun=x_path(1:T);
q_prun=x_path(T+1:2*T);
d_prun=x_path(2*T+1:3*T);

n_prun=P.sigma*( (q_prun(2:end)+z(2:end)).*(1-kh_prun(1:end-1)) - d_prun(1:end-1)   ) + P.w;

n_prun=[N_ss;n_prun];

if actual_run<120
n_prun(actual_run+1)=0;
n_prun(actual_run+2)=(1+P.sigma)*P.w;
end
ch_prun=y(2:end)-(1-P.sigma)*(n_prun(2:end)-P.w)./P.sigma - P.alpha/2*kh_prun(2:end).^2;
ch_prun=[ch_ss;ch_prun];
if actual_run<120
ch_prun(actual_run+1)=ch_run(actual_run+1);
q_norun=x_path_original(T+actual_run+1);
kh_norun=x_path_original(actual_run+1);
n_norun=P.sigma*( (q_norun+z(actual_run+1)).*(1-kh_prun(actual_run)) - d_prun(actual_run)   ) + P.w;
ch_norun=y(actual_run+1)-(1-P.sigma)*(n_norun-P.w)./P.sigma - P.alpha/2*kh_norun.^2;
end


phi_prun=q_prun.*(1-kh_prun)./max(n_prun,.0000000000001);
q_bar_prun=gamma_val*d_prun(1:end-2)./(1-kh_prun(1:end-2))-z(2:end-1);

pai_check=P.pbar*(1-min((q_run(3:end)+z(3:end)).*(1-kh_prun(2:end-1))./max(d_prun(2:end-1),.000000001),1)).^P.delta_p;
pai=P.pbar*(1-min((q_run(3:end)+z(3:end)).*max(1-kh_prun(2:end-1),.00001)./max(d_prun(2:end-1),.000000001),1)).^P.delta_p;
pai=[0;pai;pai(end)];
 
reco_r= pai(2:end-1).*P.beta.*ch_prun(2:end-1)./ch_run(3:end).*( z(3:end) + q_run(3:end)   ).*(1-kh_prun(2:end-1))./max((  phi_prun(2:end-1) -1 ).*n_prun(2:end-1),.00000000001  );
reco_rk=pai(2:end-1).*P.beta.*ch_prun(2:end-1)./ch_run(3:end).*( z(3:end) + q_run(3:end)   )./(q_prun(2:end-1)+ P.alpha*kh_prun(2:end-1));


R_prun =(1-reco_r)./ (  (1-pai(2:end-1)) .* P.beta.*ch_prun(2:end-1)./ch_prun(3:end)  );
R_prun=[1/beta;R_prun;1/beta];
if actual_run<120
R_prun(actual_run) =(1-reco_r(actual_run-1))./ (  (1-pai(actual_run)) .* P.beta.*ch_prun(actual_run)./ch_norun  );
end


Rk_prun=[1/beta+spread_val; (z(3:end) + pai(2:end-1).*q_run(3:end)+(1-pai(2:end-1)).*q_prun(3:end))./q_prun(2:end-1) ;1/beta+spread_val ];
if actual_run<120
Rk_prun(actual_run)=(z(actual_run+1) + pai(actual_run).*q_run(actual_run+1)+(1-pai(actual_run)).*x_path_original(T+actual_run+1))./q_prun(actual_run);
end
R_rf_prun=[ 1/beta; 1./ (  (1-pai(2:end-1))* P.beta.*ch_prun(2:end-1)./ch_prun(3:end) + pai(2:end-1)*P.beta.*ch_prun(2:end-1)./ch_run(3:end) ) ; 1/beta  ];
if actual_run<120
R_rf_prun(actual_run)= 1./ (  (1-pai(actual_run))* P.beta.*ch_prun(actual_run)./ch_norun + pai(actual_run)*P.beta.*ch_prun(actual_run)./ch_run(actual_run+1) );
end

Spreadk_prun=Rk_prun-R_prun;
Spreadd_prun=R_prun-R_rf_prun;
Run_p=1-(q_run(3:end)+z(3:end)).*(1-kh_prun(2:end-1))./d_prun(2:end-1);


if actual_run<120
    
    run_time=num2str(actual_run);
    save (strcat('anticipated_run_happens_at_time_',run_time))
    
else
    save anticipated_run_does_not_happen
end
    
