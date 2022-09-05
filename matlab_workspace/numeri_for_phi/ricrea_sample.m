%funzione che ricrea con interpolazione lineare i sample perduti;
function [pos,vel,tau]=ricrea_sample(pos_temp, vel_temp, tau_temp, hz)
M=pos_temp(end,1)*hz;
pos=zeros(M,8);
vel=zeros(M,8);
tau=zeros(M,8);
pos(1,:)=pos_temp(1,:);
vel(1,:)=vel_temp(1,:);
tau(1,:)=tau_temp(1,:);
j=2;
v=length(pos_temp);
for i=2:v
    k=floor(pos_temp(i,1)*hz)-floor(hz*pos_temp(i-1,1));
    for m=1:k
        pos(j,:)=pos_temp(i-1,:)+(pos_temp(i,:)-pos_temp(i-1,:))/k*m;
        vel(j,:)=vel_temp(i-1,:)+(vel_temp(i,:)-vel_temp(i-1,:))/k*m;
        tau(j,:)=tau_temp(i-1,:)+(tau_temp(i,:)-tau_temp(i-1,:))/k*m;
        j=j+1;
        100/v*i
    end
end
end