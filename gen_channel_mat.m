function H=gen_channel_mat(hp,alphap,taup,freq,time, Q_tx,Q_rx)
H=zeros(size(Q_rx,1),size(Q_tx,2));
for ip=1:length(hp)
    Q0=(time(2)-time(1))*exp(-1i*2*pi*(freq/alphap(ip)).'.*time);
    Q1=diag(exp(-1i*2*pi*freq*taup(ip)));
    H=H+hp(ip)*sqrt(1/alphap(ip))*Q_rx*Q1*Q0*Q_tx;
end
end