%funzione di inizializzazione per poter far partire i simulatori. partirà
%il simulatore inserendo il numero richiesto nella command window.
%in \funzioni son presenti diversi file per creare le traiettorie.
%è necessario RTB by peter_corke. Funzioni create con matlab 2021a 
clear all;
%addpath(genpath('RTB'))
addpath(genpath('funzioni'))            
dq0=zeros(1,7); 
dq0=ones(7,1)*10^(-7); %initial velocities
%q0= [1, -pi/4, 0, -3 * pi/4, 0, pi/2, pi/4]; %initial position
%q0=[-1.1119,0.278665, -1.5173, -2.01865, -2.44908, 3.03067,
%-0.366567]; corretto con Wf, a, b corretti
q0=zeros(1,7); %con Wf e a,b mezzi
Ts=1/100; %step time
%% limit of franka
qmax=[ 2.8973 	1.7628 	2.8973 	-0.0698 	2.8973 	3.7525 	2.8973]; %max joint position
qmin=[-2.8973 	-1.7628 	-2.8973 	-3.0718 	-2.8973 	-0.0175 	-2.8973]; %min joint position
dqmax=[2.1750 	2.1750 	2.1750 	2.1750 	2.6100 	2.6100 	2.6100] ; %max vel (min vel=-max vel)
ddqmax=[15 7.5 10 12.5 15 20 20];%max acc (min acc=-max acc)
dddqmax=[7500 3750 5000 6250 7500 10000 10000]; %max jerk
taumax=[87, 87, 87, 87, 12, 12, 12]; %max tau joint side
dtaumax=[1000,1000, 1000, 1000, 1000, 1000, 1000]; % max derivate of tau in joint side
%% internal joint control
Kgain=[ 600.0, 600.0, 600.0, 600.0, 250.0, 150.0, 50.0]; %internal stiffness gains
%% simulation input parameters 
% omega_max=0.2; %amplitude of prova_velocità_sinusoidale
% angle=pi/4; %parameter of prova_velocità_cartesiana
% vel_max=0.1; %parameter of prova_velocità_cartesianaend

fprintf('choose the desired simulator by entering the corresponding number \n-11 for simulation with position-controlled lagrange euler model with RE data \n')
fprintf( '-12 for simulation with speed-controlled lagrange euler model with RE data \n' )
fprintf('-21 for simulation with position-controlled lagrange euler model with CLS data \n' )
fprintf('-22 for simulation with speed-controlled lagrange euler model with CLS data \n' )
fprintf('-31 for simulation with position-controlled RNE model with RE data \n' )
fprintf('-32 for simulation with speed-controlled RNE model with RE data \n' )
fprintf('-41 for simulation with position-controlled RNE model with OSI data \n')
fprintf('-42 for simulation with position-controlled RNE model with OSI data \n')
prompt='\n';
Num = input(prompt)
if (Num==11)
    mdl_panda_RE(); %modello robot con EE
    run_time=8; %run time
    q_fin=[1.63724, 0.240464, -1.38321, -2.08673, -1.15816, 1.45767, -1.47525]; %per prova polinomiale
    q0= [1, -pi/4, 0, -3 * pi/4, 0, pi/2, pi/4];
    Controllo_inverso_LE_RE_posizione
end
if (Num==31)
    mdl_panda_RE(); %modello robot con EE
    run_time=8; %run time
    q_fin=[1.63724, 0.240464, -1.38321, -2.08673, -1.15816, 1.45767, -1.47525]; %per prova polinomiale
    q0= [1, -pi/4, 0, -3 * pi/4, 0, pi/2, pi/4];
    Controllo_inverso_RNE_RE_posizione
end
if (Num==21)
    mdl_panda_OSI(); %modello robot con EE
    run_time=8; %run time
    q_fin=[1.63724, 0.240464, -1.38321, -2.08673, -1.15816, 1.45767, -1.47525]; %per prova polinomiale
    q0= [1, -pi/4, 0, -3 * pi/4, 0, pi/2, pi/4];
    Controllo_inverso_LE_CLS_posizione
end
if (Num==41)
    mdl_panda_OSI(); %modello robot con EE
    run_time=8; %run time
    q_fin=[1.63724, 0.240464, -1.38321, -2.08673, -1.15816, 1.45767, -1.47525]; %per prova polinomiale
    q0= [1, -pi/4, 0, -3 * pi/4, 0, pi/2, pi/4];
    Controllo_inverso_RNE_OSI_posizione
end
if (Num==12)
    mdl_panda_RE(); %modello robot con EE
    A=[0.75 -0.75 0.75 -0.75 -0.75 0.75 -0.75];
    T=[3.5 2.0 3.0 1.75 4.0 2.0 4.0];
    run_time=43;
    q0=[0 0 0 -0.1 0.6 0.71 0.6 ];
    Controllo_inverso_LE_RE_velocita
end
if (Num==32)
    mdl_panda_RE(); %modello robot con EE
    A=[0.75 -0.75 0.75 -0.75 -0.75 0.75 -0.75];
    T=[3.5 2.0 3.0 1.75 4.0 2.0 4.0];
    run_time=43;
    q0=[0 0 0 -0.1 0.6 0.71 0.6 ];
    Controllo_inverso_RNE_RE_velocita
end
if (Num==22)
    mdl_panda_OSI(); %modello robot con EE
    A=[0.75 -0.75 0.75 -0.75 -0.75 0.75 -0.75];
    T=[3.5 2.0 3.0 1.75 4.0 2.0 4.0];
    run_time=43;
    q0=[0 0 0 -0.1 0.6 0.71 0.6 ];
    Controllo_inverso_LE_CLS_velocita
end
if (Num==42)
    mdl_panda_OSI(); %modello robot con EE
    A=[0.75 -0.75 0.75 -0.75 -0.75 0.75 -0.75];
    T=[3.5 2.0 3.0 1.75 4.0 2.0 4.0];
    run_time=43;
    q0=[0 0 0 -0.1 0.6 0.71 0.6 ];
    Controllo_inverso_RNE_OSI_velocita
end