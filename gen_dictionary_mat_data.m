function [A_dic,AD_dic_tau,ADD_dic_tau,AD_dic_beta,ADD_dic_beta] = gen_dictionary_mat_data(s_tx,freq,time,Q_rx,taup,betap,q_alpha)
alphap=q_alpha.^betap;
taps=length(taup);
A_dic=zeros(size(Q_rx,1),taps); AD_dic_tau=zeros(size(A_dic));ADD_dic_tau=zeros(size(A_dic));
AD_dic_beta=zeros(size(A_dic));ADD_dic_beta=zeros(size(A_dic));
for ip=1:taps
    Q0=(time(2)-time(1))*exp(-1i*2*pi*(freq/alphap(ip)).'.*time);
    Q1=sqrt(1/alphap(ip))*diag(exp(-1i*2*pi*freq*taup(ip)));
    Q1D_tau=sqrt(1/alphap(ip))*diag(exp(-1i*2*pi*freq*taup(ip)).*(-1i*2*pi*freq));
    Q1DD_tau=sqrt(1/alphap(ip))*diag(exp(-1i*2*pi*freq*taup(ip)).*(-1i*2*pi*freq).^2);
    
    A_dic(:,ip)=Q_rx*Q1*Q0*s_tx;
    AD_dic_tau(:,ip)=Q_rx*Q1D_tau*Q0*s_tx;
    ADD_dic_tau(:,ip)=Q_rx*Q1DD_tau*Q0*s_tx;
    Q0D_alpha=(time(2)-time(1))*exp(-1i*2*pi*(freq/alphap(ip)).'.*time).*(-1i*2*pi*freq.'.*time*(-1/alphap(ip)^2));
    Q0DD_alpha=(time(2)-time(1))*exp(-1i*2*pi*(freq/alphap(ip)).'.*time).*(-1i*2*pi*freq.'.*time*(-1/alphap(ip)^2)).^2+...
        (time(2)-time(1))*exp(-1i*2*pi*(freq/alphap(ip)).'.*time).*(-1i*2*pi*freq.'.*time*(2/alphap(ip)^3));
    Q1D_alpha=-0.5*alphap(ip)^(-3/2)*diag(exp(-1i*2*pi*freq*taup(ip)));
    Q1DD_alpha=0.5*3/2*alphap(ip)^(-5/2)*diag(exp(-1i*2*pi*freq*taup(ip)));
    Q0D_beta=(q_alpha^betap(ip))*log(q_alpha)*Q0D_alpha;
    Q1D_beta=(q_alpha^betap(ip))*log(q_alpha)*Q1D_alpha;
    Q0DD_beta=(q_alpha^betap(ip))*log(q_alpha)*Q0DD_alpha;
    Q1DD_beta=(q_alpha^betap(ip))*log(q_alpha)*Q1DD_alpha;
    AD_dic_beta(:,ip)=Q_rx*Q1D_beta*Q0*s_tx+Q_rx*Q1*Q0D_beta*s_tx;
    ADD_dic_beta(:,ip)=Q_rx*Q1DD_beta*Q0*s_tx+Q_rx*Q1D_beta*Q0D_beta*s_tx+Q_rx*Q1D_beta*Q0D_beta*s_tx+Q_rx*Q1*Q0DD_beta*s_tx;
end
