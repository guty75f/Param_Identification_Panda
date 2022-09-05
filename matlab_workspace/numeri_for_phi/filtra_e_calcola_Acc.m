%calcola l'accellerazione con differenzazione numerica e filtra tutti i
%dati esclusi i primi e gli ultimi per eliminare effetti di inizio  o di
%fine
function [pos,vel,tau,acc]=filtra_e_calcola_Acc(d, pos_temp, vel_temp, tau_temp, Dt)
pos=filtfilt(d,pos_temp(100:end,2:8)); 
vel=filtfilt(d,vel_temp(100:end,2:8));
tau=filtfilt(d,tau_temp(100:end,2:8));
pos=pos(100:end-200,:);
vel=vel(100:end-200,:);
tau=tau(100:end-200,:);
acc(1,:)=zeros(1,7);
acc(2:length(vel),:)=diff(vel)/Dt;
end