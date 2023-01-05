


%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

%var kb kh y cb ch Q Rh Rk n phi nu mu omega spread R Qr Qbar ch_run Run con w_run w_norun
var kb kb_rc kh y cb ch Q n phi spread R y_gross d 
z Z;                         % Shocks Processes

%Rk_run Rk_bar Qr run

varexo eps_z eps_mc eps_nw;  

parameters beta alpha theta gamma Z_ss k w sigma wh s f rho_z sigma_eps mc nw_start
phi_ss R_ss Rk_ss Rh_ss kh_ss kb_ss y_ss N_ss Q_ss ch_ss cb_ss mu_ss nu_ss Q_bar_ss Q_run_ss Kbar ch_run_ss y_run_ss y_gross_ss
Rk_run_ss;

%----------------------------------------------------------------
% 2. Calibration of Parameters
%----------------------------------------------------------------
%SS_Basic_Step([.0075 .965 .99 4 .005 1 1 0 .0])

rho_z=rho_z_val;
sigma_eps=rho_z^(i-1)*inno_z;
%sigma_eps=.05;
sigma_eps_p_run=.01;
%sigma_eps=.0000000001;
beta =.99;
totmc=.1*.096;
alpha =alpha_val;
gamma=.99;
k=1;
wh=.045;
sigma= sigma_val ;
Q_ss=1;
phi_ss=phi_val;     %In the paper it is 8
s=spread_val;
mc=mc_val;
f=1; %.75;
Kbar=totmc/alpha;


R_ss=1/beta;
Rk_ss=R_ss+s;
Rh_ss=1/beta;

Z_ss=Q_ss*(Rk_ss-1);
kh_ss=(Z_ss+Q_ss-Rh_ss*Q_ss)/(alpha*Rh_ss);
kb_ss=k-kh_ss;
  
N_ss=Q_ss*kb_ss/phi_ss;
w=N_ss*( 1- sigma*((Rk_ss-mc-R_ss)*phi_ss+R_ss));
nw_start=Periods_onlyhh*w*sigma;
y_gross_ss=Z_ss+wh+w;
y_ss=y_gross_ss-(alpha/2)*(kh_ss^(2))-mc*kb_ss; 

cb_ss=(1-sigma)*((Rk_ss-mc-R_ss)*phi_ss+R_ss)*N_ss;
ch_ss=y_ss-cb_ss;

theta=( (1-sigma)*beta*(R_ss+phi_ss*(s-mc)) )/(phi_ss*(1-sigma*beta*(R_ss+phi_ss*(s-mc))) );
omega_ss=1-sigma+sigma*theta*phi_ss;

nu_ss=beta*omega_ss*R_ss;
mu_ss=beta*omega_ss*s;

Q_run_ss=(beta/(1-beta))*(Z_ss-alpha*Kbar)-alpha*Kbar;
Q_bar_ss=f*Q_ss*(1-1/phi_ss)*R_ss-Z_ss;
ch_run_ss=y_gross_ss-w+(alpha/2)*(Kbar)^(2)-alpha*Kbar;
sp_run_ss=(Z_ss+Q_run_ss) /(Q_run_ss)-R_ss;

y_run_ss=y_gross_ss+(alpha/2)*(Kbar)^(2)-alpha*Kbar;
Rk_run_ss=(Z_ss+Q_run_ss)/Q_run_ss;




%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;

%%%%%%%%% PRODUCTION SIDE and MARKET CLEARING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=(kb)+(kh);

kb(-1)=kb_rc-eps_mc;

(y_gross)=Z_ss*exp(z)+w+wh*exp(z);

(y)=(y_gross)-(alpha/2)*((kh))^(2)-mc*Q(-1)*kb_rc;                 
                               
(y)=(ch)+(cb);                                                           


%%%%%%%%%% HOUSEHOLD PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%
beta*((ch)/(ch(+1)))*(R)=1;

beta*((ch)/(ch(+1)))*(Z_ss*exp(z(+1))+(Q(+1)))/( (Q) + alpha*(kh) )=1;
 
%%%%%%%%%%%% BANKER'S PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(Q)*(kb)=(phi)*(n);

(Q)*(kb)=(n)+(d)/(R);

(phi)=  beta*( (1-sigma) + sigma*theta*(phi(+1))  )*(R)  /  (  theta -   beta*( (1-sigma) + sigma*theta*(phi(+1))  )*(( Z_ss*exp(z(+1))+(Q(+1)) )/(Q)-mc-(R))  )  ;

(n)=   sigma*(  ( Z_ss*exp(z)+(Q)-mc*Q(-1) )*(kb(-1)) - (d(-1)) )  + w +eps_nw;

(cb)=(1-sigma)*(  ( Z_ss*exp(z)+(Q)-mc*Q(-1) )*(kb(-1)) - (d(-1)) );

%Qbar=(f*d/kb(-1)) -exp(z)*Z_ss; 

%%%%%%%%%%%%%%%%%% SHOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Exp Shocks

z=rho_z*z(-1)-eps_z;                            %Productivity shock  

Z=Z_ss*exp(z);

spread=4*(( Z_ss*exp(z(+1))+(Q(+1)) )/(Q)-R);
 
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

%%%%%%%

initval;

kb_rc=0;
kb=0; 
d=0;
z=0;

kh=k;
n=N_ss;

y=y_ss;
cb=cb_ss;
ch=ch_ss;
Q=Q_ss;
phi=phi_ss; 
spread=s;
R=1/beta;
y_gross=y_gross_ss;
Z=Z_ss;

end;
%steady;
%%%%%%%%


endval;


z=0;

Q=Q_ss;

kh =  (kh_ss );  

kb_rc=kb_ss;

kb =   (kb_ss); 

y =    (y_ss);

ch =    ( ch_ss   );

cb =   (cb_ss  );

n =    (N_ss );

d=((kb_ss-N_ss)/beta);

R =  (R_ss);

phi   =	log(phi_ss);

Z=Z_ss;

%Qbar=Q_bar_ss;

spread=s;

end;
steady(solve_algo=0);




%steady(solve_algo=0);
%check(solve_algo=0);

%shocks;
%var eps_z;
%stderr sigma_eps;
%end;

shocks;
var eps_z;
periods 1;
values (sigma_eps);
var eps_nw;
periods 1;
values (nw_start);


end;


simul(periods=120, stack_solve_algo=0);
%stoch_simul(order=1,dr_algo=0,periods=1000,irf=60) z y cb ch kh kb n Q R spread phi ;





