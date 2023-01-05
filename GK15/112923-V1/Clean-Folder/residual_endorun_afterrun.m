function res=residual_endorun_afterrun(x)

global P y z q_run ch_run T mc_val ch_ss phi_val

kh=x(1:T);
q=x(T+1:2*T);
d=x(2*T+1:3*T);

n=P.sigma*( (q(2:end)+z(2:T)-mc_val*q(1:end-1)).*(1-kh(1:end-1)) - d(1:end-1)   ) + P.w;

n=[P.w;n];
n(2)=n(2)+(P.sigma)*P.w;

ch=y(2:T)-(1-P.sigma)*(n(2:end)-P.w)./P.sigma - P.alpha/2*kh(2:end).^2-mc_val*q(1:end-1).*(1-kh(1:end-1));
ch=[ch_run(1);ch];
phi=q(2:end).*(1-kh(2:end))./n(2:end);
phi=[phi_val;phi];
%q_bar=gamma_val*d(1:end-1)./(1-kh(1:end-1))-z(2:end)+mc_val*q(1:end-1);

%p=max(P.pbar*(1-(q_run(3:end)+z(3:T)).*(1-kh(2:end-1))./d(2:end-1)).^P.delta_p,0);
p=P.pbar*(1-min((q_run(3:end)+z(3:T)).*(1-kh(2:end-1))./d(2:end-1),1)).^P.delta_p;
p=[0;p;p(end)];
reco_r= p(2:end-1).*P.beta.*ch(2:end-1)./ch_run(3:end).*( z(3:T) + q_run(3:end)   ).*(1-kh(2:end-1))./((  phi(2:end-1) -1 ).*n(2:end-1)  );
reco_rk=p(2:end-1).*P.beta.*ch(2:end-1)./ch_run(3:end).*( z(3:T) + q_run(3:end)   )./(q(2:end-1)+ P.alpha*kh(2:end-1));


R =(1-reco_r)./ (  (1-p(2:end-1)) .* P.beta.*ch(2:end-1)./ch(3:end)  );

res=[ kh(end)-P.kh_ss_prun;
    
     q(end)-P.q_ss_prun;
    
     d(end)-P.d_ss_prun;
     
     kh(1)-1;
    
     q(1)-1;
    
     d(1)-0;


      (P.sigma*n(2:end-1).*(  phi(2:end-1).*(   (z(3:T) + q(3:end))./q(2:end-1)-mc_val - R ) + R  )+ P.w-n(3:end)  ) ;   
    
      ((1-p(2:end-1)) .* P.beta.*ch(2:end-1)./ch(3:end).*( z(3:T) + q(3:end)   )./(q(2:end-1)+ P.alpha*kh(2:end-1))+reco_rk-1);
      
      (phi(2:end-1)-(1-p(2:end-1))./P.theta.*P.beta.*(  phi(2:end-1).*(   (z(3:T) + q(3:end))./q(2:end-1) - R -mc_val) + R  ).*(1-P.sigma+P.sigma*P.theta.*phi(3:end)))

        
          ];

end



