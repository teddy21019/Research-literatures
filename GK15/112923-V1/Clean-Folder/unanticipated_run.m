clear all
global P y z N_ss pai ch_run q_run T kh_ss d_ss mc_val ch_ss phi_val gamma_val
tic;


%% calibration 

% parameters, targets and steady state
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
gamma_val=1; 
rho_z=rho_z_val;
beta =.99;
totmc=.1*.096;
alpha =alpha_val;
k=1;
wh=.045;
sigma= sigma_val ;
Q_ss=1;
phi_ss=phi_val;     %In the paper it is 8
s=spread_val;
mc=mc_val; %disregard this

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


% Experiment parameters 
t_run=3; % Period when the run happens
t_run_dev=t_run-1;% Period in which deviation is checked 
antirun=1;
Periods=120;% Simulation Periods
Periods_onlyhh=1; %Periods after ther run during which banks don't operate
T=Periods+2; 
N=40; %Plot periods

%% Initialize output 

q_run=zeros(Periods,1);
ch_run=zeros(Periods,1);
eps_t=zeros(Periods,1);


%% Recession

i=1;
%%Get Impulse responses to a shock to Z from dynare
dynare PF_level_loop_zshock noclearall

q_z=Q;
kh_z=kh;
d_z=d;
phi_z=phi;
y_z=y;
n_z=n;
R_z=R;
ch_z=ch;
cb_z=cb;
x_path_norun=[kh_z;q_z;d_z];

%%check convergence
convergence_zshock=oo_.deterministic_simulation.status;

%%Find the threshold price at each time
q_bar_norun=gamma_val*d_z(1:end-2)./(1-kh_z(1:end-2))-exp(z(2:end-1))*Z_ss;

%% Unanticipated Run


%%initialize output
q_zrun=ones(Periods,1);

%%Compute paths for the price of capital and consumption upon a run (q_run and ch_run) and the
%%ex-post run impulse responses
for j=1:Periods
    
    %%Starting from the end of the simulation and working backwards compute
    %%the path of the economy after a run happens back to steady state
     i=(Periods+Periods_onlyhh-j+1);
     dynare PF_level_loop_convex noclearall
    %%Use the household euler equation and the resource constraint to
    %%derive consumption and asset prices when bankers are inactive
     z_ini=exp(-rho_z.^(i-Periods_onlyhh-1:i-1)'*inno_z)*Z_ss;
     ch_only_hh=[exp(-rho_z.^(i-Periods_onlyhh-1:i-2)'*inno_z)*(Z_ss+wh)-alpha_val/2;ch(2)];
     q_zrun(Periods_onlyhh+1:end)=Q(2:Periods-Periods_onlyhh+1);
     
     for hh_only=1:Periods_onlyhh
         t=Periods_onlyhh+1-hh_only;
         q_zrun(t)=beta*ch_only_hh(t)/ch_only_hh(t+1)*(z_ini(t+1)+q_zrun(t+1))-alpha_val;
     end
         
         
     eps_t(i-Periods_onlyhh)=rho_z.^(i-Periods_onlyhh-1)*inno_z;
     q_run(i-Periods_onlyhh)=q_zrun(1);
     ch_run(i-Periods_onlyhh)=ch_only_hh(1);
     convergence_dynare_run(i)=oo_.deterministic_simulation.status;
     
     %%Once we have worked backwards to time of the run we construct the ex-post run paths    
     if j==Periods-(t_run-1)
         q_zrunn=[1;q_z(2:t_run);q_zrun;1]; 
         
kh_zrun=[kh_ss;kh_z(2:t_run);ones(Periods_onlyhh,1);kh(2:Periods-Periods_onlyhh+2)];
d_zrun=[0;d_z(2:t_run);zeros(Periods_onlyhh,1);d(2:Periods-Periods_onlyhh+2)];
phi_zrun=[phi_ss;phi_z(2:t_run);zeros(Periods_onlyhh,1);phi(2:Periods-Periods_onlyhh+2)];
y_zrun=[y_ss;y_z(2:t_run);exp(-rho_z.^(0:Periods_onlyhh-1)'*inno_z)*(Z_ss+wh)+w-alpha_val/2;y(2:Periods-Periods_onlyhh+2)];
n_zrun=[N_ss;n_z(2:t_run);(1:Periods_onlyhh)'*w;n(2:Periods-Periods_onlyhh+2)];
ch_zrun=[ch_ss;ch_z(2:t_run);ch_only_hh;ch(3:Periods-Periods_onlyhh+2)];
R_zrun=[1/beta;R_z(2:t_run);1/beta*ch_zrun(t_run+2:end)./ch_zrun(t_run+1:end-1);1/beta];
cb_zrun=[cb_ss;cb_z(2:t_run);zeros(Periods_onlyhh,1);cb(2:Periods-Periods_onlyhh+2)];


     end
 
     j
end


%%% Adjust Unanticipated Run Variables
q_zrun=q_zrunn(1:end-(t_run-1));
kh_zrun=kh_zrun(1:end-(t_run-1));
d_zrun=d_zrun(1:end-(t_run-1));
phi_zrun=phi_zrun(1:end-(t_run-1));
y_zrun=y_zrun(1:end-(t_run-1));
n_zrun=n_zrun(1:end-(t_run-1));
ch_zrun=ch_zrun(1:end-(t_run-1));
R_zrun=R_zrun(1:end-(t_run-1));
cb_zrun=cb_zrun(1:end-(t_run-1));
q_run_plot=q_run;
q_run_plot=[q_run_plot(1);q_run_plot;q_run_plot(end)];

save unanticipated



