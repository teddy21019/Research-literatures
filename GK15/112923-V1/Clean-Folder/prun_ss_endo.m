function [kh_ss_prun d_ss_prun q_ss_prun]=prun_ss_endo(xo,q_run)
global P ch_run



function res=residual_ss(x)



kh=x(1);
q=x(2);
d=x(3);

y=P.Z_ss+P.wh+P.w;
n=P.sigma*( (q+P.Z_ss).*(1-kh) - d  ) + P.w;
ch=y-(1-P.sigma)*(n-P.w)./P.sigma - P.alpha/2*kh.^2;
phi=q.*(1-kh)./n;
q_bar=d./(1-kh)-P.Z_ss;

%p=max(P.pbar*(1-(q_run+P.Z_ss).*(1-kh)./d).^P.delta_p,0);
p=P.pbar*(1-min((q_run+P.Z_ss).*(1-kh)./d,1)).^P.delta_p;



reco_r= p.*P.beta.*ch./ch_run(1).*( P.Z_ss + q_run   ).*(1-kh)./((  phi -1 ).*n  );
reco_rk=p.*P.beta.*ch./ch_run(1).*( P.Z_ss + q_run   )./(q+ P.alpha*kh);


R =(1-reco_r)./ (  (1-p) .* P.beta.*ch./ch  );

res=[ 
      (P.sigma*n.*(  phi.*(   (P.Z_ss + q)./q- R ) + R  )+ P.w-n  ) ;   
    
      ((1-p) .* P.beta.*ch./ch.*( P.Z_ss+ q   )./(q+ P.alpha*kh)+reco_rk-1);
      
      (phi-(1-p)./P.theta.*P.beta.*(  phi.*(   (P.Z_ss + q)./q - R) + R  ).*(1-P.sigma+P.sigma*P.theta.*phi))

        
          ];

end

x=fsolve(@residual_ss,xo);
kh_ss_prun=x(1);
q_ss_prun=x(2);
d_ss_prun=x(3);

end