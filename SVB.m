function [hsparse,tauGrid_next,betaGrid_next,iter]=SVB(y,tol,maxiter,noisevar,a,b,mu_old,tauGrid,betaGrid,q_alpha,K,resolu_tau,resolu_beta,s_tx,freq,time,Q_rx)
[A,AD_tau,ADD_tau,AD_beta,ADD_beta]=gen_dictionary_mat_data(s_tx,freq,time,Q_rx,tauGrid,betaGrid,q_alpha);
Phi=A;PhiD_tau=AD_tau;PhiDD_tau=ADD_tau;
PhiD_beta=AD_beta;PhiDD_beta=ADD_beta;
row_len=size(A,1);col_len=size(A,2);tauMax=max(tauGrid);betaMax=max(betaGrid);
c=1e-6;d=1e-6;
%gamma=c/d;
% mu_old=zeros(col_len,1);
Sigma=mu_old*mu_old';
% alpha_mean=1./(A'*y);%(a/b)*ones(col_len,1);
alpha_mean=(a+1)./(b+real(diag(Sigma))+(abs(mu_old).^2));
D_alpha=diag(alpha_mean);
gamma=(row_len+c)/(d+norm(y-Phi*mu_old,2)^2+real(trace(Phi*Sigma*Phi')));
% gamma=1/noisevar;
converged = false;
iter=1;tauGrid_next=tauGrid;
betaGrid_next=betaGrid;
while ~converged
    Sigma=(D_alpha+Phi'*Phi*gamma)^(-1);
    mu=Sigma*Phi'*y*gamma;%zeros(col_len,1);
    if norm(mu - mu_old)/norm(mu) < tol || iter >= maxiter
        converged = true;
    end
    mu_old=mu;iter=iter+1;
    [~,pos]=maxk(mu,K,'ComparisonMethod','abs');
    pos=sort(pos);
    %% tau refinement
    temp=tauGrid_next;tauGrid_next=tauGrid;tauGrid_next(pos)=temp(pos);
    temp1=Phi;Phi=A;Phi(:,pos)=temp1(:,pos);temp2=PhiD_tau;PhiD_tau=AD_tau;PhiD_tau(:,pos)=temp2(:,pos);temp3=PhiDD_tau;PhiDD_tau=ADD_tau;PhiDD_tau(:,pos)=temp3(:,pos);
    for i = 1:1
        for n = 1:K
            index=pos(n);
            der1 = -2*real(mu(index) * (y-Phi*mu)'*PhiD_tau(:,index))+2*real(PhiD_tau(:,index).'*conj(Phi)*Sigma(index,:).');%+PhiD(:,index)'*Phi*Sigma(:,index));
            der2 = -2*real(mu(index) * (y-Phi*mu)'*PhiDD_tau(:,index)) +2*abs(mu(index))^2*(PhiD_tau(:,index)'*PhiD_tau(:,index))+...
                real(PhiDD_tau(:,index).'*conj(Phi)*Sigma(index,:).'+Sigma(index,index)*PhiD_tau(:,index).'*conj(PhiD_tau(:,index))+...
                PhiDD_tau(:,index)'*Phi*Sigma(:,index)+Sigma(index,index)*PhiD_tau(:,index)'*PhiD_tau(:,index));
            if der2 > 0
                tau_next = tauGrid_next(index) - der1/der2;
            else
                tau_next = tauGrid_next(index)-sign(der1)*(tauMax/1e3)*rand(1);%tau;%sign(der1)*(tauMax/20)*rand(1);
            end
            if tau_next>=0 && tau_next<=tauMax
                if tau_next>=tauGrid(index)-resolu_tau/2 && tau_next<=tauGrid(index)+resolu_tau/2
                    [a_next,ad_next_tau,add_next_tau,~,~]=gen_dictionary_mat_data(s_tx,freq,time,Q_rx,tau_next,betaGrid_next(index),q_alpha);

                                        if norm(y-Phi*mu)>norm(y-Phi*mu+Phi(:,index)*mu(index)-a_next*mu(index))
                    tauGrid_next(index)=tau_next;
                    Phi(:,index)=a_next;PhiD_tau(:,index)=ad_next_tau;PhiDD_tau(:,index)=add_next_tau;

                                        end
                end
            end
        end
       
    end

   %% beta refinement
   temp=betaGrid_next;betaGrid_next=betaGrid;betaGrid_next(pos)=temp(pos);
    temp1=Phi;Phi=A;Phi(:,pos)=temp1(:,pos);temp2=PhiD_beta;PhiD_beta=AD_beta;PhiD_beta(:,pos)=temp2(:,pos);temp3=PhiDD_beta;PhiDD_beta=ADD_beta;PhiDD_beta(:,pos)=temp3(:,pos);
    for i = 1:1
        for n = 1:K
            index=pos(n);         
            der1 = -2*real(mu(index) * (y-Phi*mu)'*PhiD_beta(:,index))+2*real(PhiD_beta(:,index).'*conj(Phi)*Sigma(index,:).');%+PhiD(:,index)'*Phi*Sigma(:,index));
            der2 = -2*real(mu(index) * (y-Phi*mu)'*PhiDD_beta(:,index)) +2*abs(mu(index))^2*(PhiD_beta(:,index)'*PhiD_beta(:,index))+...
                real(PhiDD_beta(:,index).'*conj(Phi)*Sigma(index,:).'+Sigma(index,index)*PhiD_beta(:,index).'*conj(PhiD_beta(:,index))+...
                PhiDD_beta(:,index)'*Phi*Sigma(:,index)+Sigma(index,index)*PhiD_beta(:,index)'*PhiD_beta(:,index));
            if der2 > 0
                beta_next = betaGrid_next(index) - der1/der2;
            else
                beta_next = betaGrid_next(index)-sign(der1)*(2*betaMax/1e3)*rand(1);
            end
            if beta_next>=-betaMax &&   beta_next<=betaMax
                if beta_next>=betaGrid(index)-resolu_beta/2 && beta_next<=betaGrid(index)+resolu_beta/2
                    [a_next,~,~,ad_next_beta,add_next_beta]=gen_dictionary_mat_data(s_tx,freq,time,Q_rx,tauGrid_next(index),beta_next,q_alpha);

                                        if norm(y-Phi*mu)>norm(y-Phi*mu+Phi(:,index)*mu(index)-a_next*mu(index))
                    betaGrid_next(index)=beta_next;
                    Phi(:,index)=a_next;PhiD_beta(:,index)=ad_next_beta;PhiDD_beta(:,index)=add_next_beta;

                                        end
                end
            end
        end
       
    end

    gamma=(row_len+c)/(d+norm(y-Phi*mu,2)^2+real(trace(Phi*Sigma*Phi')));
    alpha_mean=(a+1)./(b+real(diag(Sigma))+(abs(mu).^2));
    D_alpha=diag(alpha_mean);
    Sigma=(D_alpha+Phi'*Phi*gamma)^(-1);
    mu_new=Sigma*Phi'*y*gamma;

end

hsparse=mu_new;
end
