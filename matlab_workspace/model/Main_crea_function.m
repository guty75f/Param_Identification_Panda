%% crea le funzioni per B, G e vettore di coriolis
% una volta calcolati i risultati grazie ai dati far partire questa parte
% inserendo i parametri base calcolati non semplificati per ottenere le funzioni di Lagrange
% Eulero con quei parametri
%
clear all

funzione_coriolis(JJ,N,C2,m,com,tipogiunto,theta,dq,PI);
funzione_inerzia(N,tipogiunto,theta,JJ, B2, com, m,PI);
funzione_gravity(theta, tipogiunto, G2, N, com, m,g,PI);
