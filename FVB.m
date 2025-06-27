function [hsparse,tauGrid_next,betaGrid_next,iter]=FVB(y,tol,maxiter,noisevar,a,b,mu_old,tauGrid,betaGrid,K,resolu,s_tx,freq,time,Q_rx,q_alpha)

[A,AD_tau,~,AD_beta,~] = gen_dictionary_mat_data(s_tx,freq,time,Q_rx,tauGrid,betaGrid,q_alpha);
row_len=size(A,1);col_len=size(A,2);Phi=A;tauMax=max(tauGrid);betaMax=max(betaGrid);
c=1e-6;d=1e-6;
% gamma=c/d;
% mu_old=zeros(col_len,1);
Sigma=mu_old*mu_old';
% alpha_mean=1./abs(A'*y);%(a/b)*ones(col_len,1);
alpha_mean=(a+1)./(b+real(diag(Sigma))+(abs(mu_old).^2));
D_alpha=diag(alpha_mean);
% gamma=1/noisevar;
gamma=(row_len+c)/(d+norm(y-Phi*mu_old,2)^2+real(trace(Phi*Sigma*Phi')));
converged = false;
iter=1;del_tau=zeros(size(tauGrid));tauGrid_next=tauGrid;
del_beta=zeros(size(betaGrid));betaGrid_next=betaGrid;
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
    temp=del_tau; del_tau=zeros(size(tauGrid));del_tau(pos)=temp(pos);
    P_ =real((AD_tau'*AD_tau) .*((mu*mu').')+Sigma.'.*(AD_tau'*AD_tau));%real((PhiD'*PhiD).*(mu*(mu.')+Sigma.'));%
    P=P_(pos,pos);
    v_=(real((y-Phi*mu)'*AD_tau*diag(mu))-real(diag(Sigma*Phi'*AD_tau).')).';%real(diag(mu)*PhiD.'*conj(y)-diag((mu*mu'+Sigma)*Phi'*PhiD));
    v=v_(pos);
    temp1 = P \ v;
    if any(abs(temp1)>resolu/2) || any(diag(P)==0)
        for i = 1:5
            for n = 1:K
                temp_del = del_tau(pos);
                temp_del(n) = 0;
                del_tau(pos(n)) = (v(n) - (P(n,:) * temp_del.')) / P(n,n);
                if del_tau(pos(n)) > resolu/2
                    del_tau(pos(n)) = resolu/2;
                end
                if del_tau(pos(n)) < -resolu/2
                    del_tau(pos(n)) = -resolu/2;
                end
                if P(n,n) == 0
                    del_tau(pos(n)) = 0;
                end

            end
        end
    else
        del_tau=zeros(size(tauGrid));
        del_tau(pos)=temp1;
    end
    tauGrid_next1=tauGrid_next+del_tau;
    %% beta refinement
    temp=del_beta; del_beta=zeros(size(betaGrid));del_beta(pos)=temp(pos);
    P_ =real((AD_beta'*AD_beta) .*((mu*mu').')+Sigma.'.*(AD_beta'*AD_beta));%real((PhiD'*PhiD).*(mu*(mu.')+Sigma.'));%
    P=P_(pos,pos);
    v_=(real((y-Phi*mu)'*AD_beta*diag(mu))-real(diag(Sigma*Phi'*AD_beta).')).';%real(diag(mu)*PhiD.'*conj(y)-diag((mu*mu'+Sigma)*Phi'*PhiD));
    v=v_(pos);
    temp1 = P \ v;
    if any(abs(temp1)>resolu/2) || any(diag(P)==0)
        for i = 1:5
            for n = 1:K
                temp_beta = del_beta(pos);
                temp_beta(n) = 0;
                del_beta(pos(n)) = (v(n) - (P(n,:) * temp_beta.')) / P(n,n);
                if del_beta(pos(n)) > resolu/2
                    del_beta(pos(n)) = resolu/2;
                end
                if del_beta(pos(n)) < -resolu/2
                    del_beta(pos(n)) = -resolu/2;
                end
                if P(n,n) == 0
                    del_beta(pos(n)) = 0;
                end

            end
        end
    else
        del_beta=zeros(size(betaGrid));
        del_beta(pos)=temp1;
    end
    betaGrid_next1=betaGrid_next+del_beta;
    tauGrid_next=tauGrid_next+del_tau;
    tauGrid_next(tauGrid_next<=0)=0;tauGrid_next(tauGrid_next>=tauMax)=tauMax;
    betaGrid_next=betaGrid_next+del_beta;
    betaGrid_next(betaGrid_next<=-betaMax)=-betaMax;betaGrid_next(betaGrid_next>=betaMax)=betaMax;
    for ii=1:length(pos)
        %         if (tauGrid(pos(ii))-resolu/2) <=tauGrid_next1(pos(ii)) && tauGrid_next1(pos(ii))<= (tauGrid(pos(ii))+resolu/2)
        [a_next,~,~,~,~]=gen_dictionary_mat_data(s_tx,freq,time,Q_rx,tauGrid_next1(pos(ii)),betaGrid_next1(pos(ii)),q_alpha);

        Phi(:,pos(ii))=a_next;
        %         else
        %             del_tau(pos(ii))=temp(pos(ii));
        %         end
    end



    
    gamma=(row_len+c)/(d+norm(y-Phi*mu,2)^2+real(trace(Phi*Sigma*Phi')));
    alpha_mean=(a+1)./(b+real(diag(Sigma))+(abs(mu).^2));
    D_alpha=diag(alpha_mean);
    Sigma=(D_alpha+Phi'*Phi*gamma)^(-1);
    mu_new=Sigma*Phi'*y*gamma;
end

hsparse=mu_new;
end