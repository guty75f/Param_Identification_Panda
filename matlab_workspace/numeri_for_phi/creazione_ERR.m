%calcola l'errore tra il modello dinamico senza attrito e i valori misurati
%del cobot.
clear all
addpath("Dati")
f=1000; %frequenza di campionamento
dT=10*0.15*pi*5/(2*pi*(f/2)); %half power frequency per filtro di tipo butter
d1 = designfilt('lowpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',dT,'DesignMethod','butter');
%creazione filtro di BUTTERWORTH
%% ----------------------------
%   lettura e filtraggio dati 
%  ---------------------------
%in caso di aggiunta di nuovi dati oltre al primo aggiungere le variabili
%chiamandole 2 e aggiungerle con ; a q,dq,ddq e tau_stack;
pos1_temp=importdata('Dati/a_posizione_giunto_test_ID_friction_test_4.txt');
vel1_temp=importdata('Dati/a_velocità_giunto_test_ID_friction_test_4.txt');
tau1_temp=importdata("Dati/a_tau_J_test_ID_friction_test_4.txt");

[pos1,vel1,tau1]=ricrea_sample(pos1_temp,vel1_temp,tau1_temp,f); 
%ricrea i sample perduti dai dati presi

[pos1,vel1,tau1,acc1]=filtra_e_calcola_Acc(d1,pos1,vel1,tau1,1/f);
%filtra e calcola con differenzazione numerica l'accellerazione dei dati presi
 q=[pos1];
 dq=[vel1];
 ddq=[acc1];
 tau_stack=[tau1];

v=length(dq);
EE=zeros(v,7);
taux=zeros(v,7);
%% calcolo dell'errore
for i=1:v
   C_temp=Vettore_coriolis(dq(i,:),q(i,:)); %calcolo della componente di Coriolis
   G_temp=Vettore_Gravita(q(i,:),9.81); %calcolo della componente gravitazionale
   B=Matrice_inerzia(q(i,:)); %matrice inerziale
   M_temp=B*ddq(i,:)'; %calcolo componente inerziale
   tau_temp=G_temp+M_temp+C_temp; %calcolo delle coppia date dal modello dinamico senza attrito
   EE_temp=tau_stack(i,:)'-tau_temp; %calcolo dell'errore tra il modello dinamico senza attrito e i valori misurati
   EE(i,:)=EE_temp';
   taux(i,:)=tau_temp'; 
   100/v*i %debug con percentuale di avanzamento
end
save('friction_test_x', 'ddq', 'dq', 'EE', 'q', 'tau_stack', 'taux'); %salvataggio dati