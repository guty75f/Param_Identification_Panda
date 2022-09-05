%crea il regressore numerico classico e risolve il OLS problem per il
%metodo classico. Richiede come dati dal cobot le posizioni, le velocità e
%le coppie in ogni istante.

clear all
%l'accellerazione può essere calcolata con la velocità filtrata o con la
%velocità non filtrata per poi filtrare il risultato
Acc_filtrata_direttamente=true; 

% addpath("dati_jointxjoint\");
% pos=importdata('dati_jointxjoint\a_posizione_giunto_test_ID_speriamo_in_DIO1.txt');
% vel=importdata('dati_jointxjoint\a_velocità_giunto_test_ID_speriamo_in_DIO1.txt');
% torque=importdata('dati_jointxjoint\a_tau_J_test_ID_speriamo_in_DIO1.txt');

f=1000; %frequenza di campionamento
HPF=10*0.15*pi*5/(2*pi*f/2); %Half power frequency del Butterworth filter
D1=designfilt('lowpassiir','FIlterOrder',4,...
    'HalfPowerFrequency',HPF,'DesignMethod','butter');
[pos,vel,tau_class]=ricrea_sample(pos,vel,torque,f); 
v=length(pos);
% % modalità che calcola prima l'acc e poi la filtra: 
if Acc_filtrata_direttamente
acc1(1,:)=zeros(1,7);
for i=1:v-1
        acc1(i+1,:)=(vel(i+1,2:8)-vel(i,2:8))*f;
end
acc1=filtfilt(D1,acc1);
acc1=acc1(199:end-200,:);
[pos,vel,tau_class,acc]=filtra_e_calcola_Acc(D1,pos,vel,tau_class,1/f); %filtro e calcolo l'acc 
v=length(pos);
Ytot_class=zeros(v*7,64);
for i=1:v   
    Yr_temp=compute_Y_classical_simplified(acc1(i,:),vel(i,:),pos(i,:),9.81);
    Ytot_class(7*i-6:7*i,:)=Yr_temp; %salvo il regressore ottenuto dai dati
    creazione_YR=i*100/v        %debug percentuale
end
% prima filtro la velocità poi calcolo l'accellerazione
else
[pos,vel,tau_class,acc]=filtra_e_calcola_Acc(D1,pos,vel,tau_class,1/f); %filtro e calcolo l'acc 
v=length(pos);
Ytot_class=zeros(v*7,64);
for i=1:v   
    Yr_temp=compute_Y_classical_simplified(acc(i,:),vel(i,:),pos(i,:),9.81);
    Ytot_class(7*i-6:7*i,:)=Yr_temp; %salvo il regressore ottenuto dai dati
    creazione_YR=i*100/v        %debug percentuale
end    
end
 for i=1:v 
    tau_list_class(i*7-6:i*7,:)=tau_class(i,:)'; %rendo tau_class un vettore 
 end
PI_classic=inv(Ytot_class'*Ytot_class)*Ytot_class'*tau_list_class;
%calcolo i valori non realistici ottenuti con il classical method

Y_stack_LI=Ytot_class;
tau_stack=tau_list_class';
save("Y_e_tau_class_2.mat","Y_stack_LI","tau_stack");

[perc_class]=calcolo_variazione_percentuale_dati(PI_classic,tau_list_class, Ytot_class);
save("PI_class_2.mat", "PI_classic", "perc_class");