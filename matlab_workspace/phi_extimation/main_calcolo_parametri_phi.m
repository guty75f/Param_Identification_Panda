%funzione per calcolare l'attito con il modello sigmoideo
% author: Riccardo Bisogni , Gatti Gabriele
% date: 26_07_2022
%  deve contenere i dati EE e dq estratti da creazione_ERR

close all
clc
clear all
global EE dq SA_step TESTNUM 
% inserire qui il test da cui estrarre i parametri, deve contenere i dati
% EE e dq estratti da creazione_ERR


%% ----------------------------------
% inizializzazione
%----------------------------------
num_x=3; %numero di variabili del modello
LI= [-1;-100;-0.1;]; %limite inferiore
LS=[+1;100;0.1;]; %limite superiore
N=7;    %numero di giunti
N_prove_indipendenti = 2; % numero di prove indipendenti con X0 diversi
num_of_SA_steps = 2; % dalla prima soluzione quante volte cercare di migliorarlo(non davvero necessario con fminsearch)
    
%% ----------------------------------
% algoritmo
%----------------------------------

for T=1:N
TESTNUM=T;
%vettore che contiene le perdite di tutti i test fatti giunto per giunto,
% inizializzato a infinito
LOSSES = Inf*ones(N_prove_indipendenti,1);

% SOL contiene i valori i valori che ottimizzano al meglio il modello
SOL = zeros(num_x,N_prove_indipendenti);

% OUTPUTS and EXITFLAGS are variables related to 'fminsearch' Matlab
% function
OUTPUTS = cell(N_prove_indipendenti,1);
EXITFLAGS = zeros(N_prove_indipendenti,1);

%----------------------------------
% optimizator: fminsearch con Nelder_Mead optimization
%----------------------------------

for i=1:N_prove_indipendenti
    for SA_step=1:num_of_SA_steps
        stringtodisp = sprintf('prova %d di %d, step %d of %d del giunto %d su %d',i,N_prove_indipendenti,SA_step,num_of_SA_steps, TESTNUM, N);
        disp(stringtodisp);
        if SA_step==1
            X0 = rand(num_x,1).*(LS-LI) + LI;
        end
      
            [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@try_this_things,X0); %funzione di ottimizzazione.
        X0 = X;
    end
    disp('---------------------------');
    disp(sprintf('.........LOSS = %f',FVAL)); %debugga la perdite per il test
    disp('---------------------------');
    LOSSES(i) = FVAL;
    SOL(:,i) = X;
    OUTPUTS{i} = OUTPUT;
    EXITFLAGS(i) = EXITFLAG;
  
end

% Estraggo il miglior test: quello con il loss minore
min_idx = find(LOSSES==min(LOSSES));
optimal_solution = SOL(:,min_idx);

% salvo i dati in un vettore per ogni parametro del modello
phi1(TESTNUM) = optimal_solution(1);
phi2(TESTNUM) = optimal_solution(2);
phi3(TESTNUM) = optimal_solution(3);

end
% salvo i dati del test con i valori phi calcolati
save('phi_con_test_x', 'phi1', 'phi2', 'phi3');
%----------------------------------
% Auto Validazione
%----------------------------------

% global EE dq TESTNUM
[MSE]=plot_risultati(phi1,phi2,phi3,tau_stack,taux);
%mostra i risulati con delle figure per ogni giunto
% e permette di calcolare il MSE (mean square error)
