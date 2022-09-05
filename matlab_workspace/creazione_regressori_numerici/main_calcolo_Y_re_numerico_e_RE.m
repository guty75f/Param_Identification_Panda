%crea il regressore numerico per il reverse eng e risolve il OLS problem per il
%metodo reverse. Richiede come dati dal cobot le posizioni, le velocità 
%le coppie, i valori di inerzia e di gravità in ogni istante.

clear all

%% caricamento dati

% pos=importdata('dati_jointxjoint\a_posizione_giunto_test_ID_3.txt');
% grav=importdata('dati_jointxjoint\a_gravity_test_ID_3.txt');
% inertia=importdata('dati_jointxjoint\a_inertia_test_ID_3.txt');
% torque=importdata('dati_jointxjoint\a_tau_J_test_ID_3.txt');
% vel=importdata('dati_jointxjoint\a_velocità_giunto_test_ID_3.txt');

%% calcolo regressore e PI
v=length(vel);
Yre_temp=zeros(35,43);
Yre_tot=zeros(35*v,43);
%creazione del regressore per il reverse eng
for i=1:v
    q=pos(i,2:8);
    Yre_temp=compute_Y_re_simplified(q,9.81);
    Yre_tot(35*i-34:35*i,:)=Yre_temp; 
    creazione_Y_RE=(i/v)*100 %debug
end
Inertia= inertia(:,2:50);
gravity=grav(1:v,2:8);

%creo il vettore d'inerzia in modo da aver solo i valori della matrice
%superiore di inerzia in quanto simmetri
for i=1:v
    XX=reshape(Inertia(i,:),[7,7]);
    XX=triu(XX);
    YY=nonzeros(XX);
    Iner(i,1:28)=YY';
end
% in questo caso il nostro tau in tau_re=Yre*PI_re è formato dai valori
% della matrice superiore di inerzia (Iner) e i valori del vettore di gravità
% (gravity)
tau=[Iner,gravity];

%il "tau" viene messo in un vettore
for i=1:v
   tau_list_re(i*35-34:i*35)=tau(i,:)';
end
%calcolo dei parametri attraverso LSM
PI_re=inv(Yre_tot'*Yre_tot)*Yre_tot'*tau_list_re';


%% deviazione standard
[perc,sigma]=calcolo_variazione_percentuale_dati(PI_re,tau_list_re', Yre_tot);
Y_stack_LI=Yre_tot;
tau_stack=tau_list_re';
% tau_stack=tau_list_re;
save("Y_e_tau_re.mat", "Y_stack_LI","tau_stack", "PI_re");
save("risultati_RE", "perc", "PI_re");

