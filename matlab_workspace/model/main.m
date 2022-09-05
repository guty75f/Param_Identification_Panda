clear all
%% inizializzazione. far partire solo questa se si hanno già matrici, regressori salvati
DH=1; %DH=0: parametri classici; DH=1 parametri modificati
%definizione incognite
N=7; %numero di giunti
syms q1 q2 q3 q4 q5 q6 q7 real;
tipogiunto=[1, 1, 1, 1, 1, 1, 1]; % 1 se rotazionale, 0 se prismatico


Tool=[1 0 0 0;
    0 1 0 0;
    0 0 1 0.107;
    0 0 0 1];         %tool parameters
% classical DH Parameters 
if(DH==0)
    a=[0,0,0.0825,-0.0825,0,0.088,0,0]; %valori DH a (spostamento su asse x)
    d=[0.333,0,0.316,0,0.384,0,0.107,0.107]; % valori DH d (spostamento su asse z)
    alpha=[-pi/2,pi/2,pi/2,-pi/2,pi/2,pi/2,0,0]; % valori DH di alpha (rotazione su asse x)
    theta=[q1,q2,q3,q4,q5,q6,q7]; %valori DH di theta (rotazione su asse z)

end


% Modified Dh Parameters
if (DH==1)
    a=[0,0,0,0.0825,-0.0825,0,0.088]; %valori MDH a (spostamento su asse x)
    d=[0.333,0,0.316,0,0.384,0,0]; % valori MDH d (spostamento su asse z)
    alpha=[0,-pi/2,pi/2,pi/2,-pi/2,pi/2,pi/2]; % valori MDH di alpha (rotazione su asse x)
    theta=[q1,q2,q3,q4,q5,q6,q7]; %valori MDH di theta (rotazione su asse z)
end
alphad=rad2deg(alpha); %per poter mettere avere sin(0)=0 e sin(pi/2)=1 nei calcoli vengono usati i deg al posto dei rad

% velocità accelerazioni
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real;  %velocità di ogni giunto
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real; % accellerazione di ogni giunto
dq= [dq1,dq2,dq3,dq4,dq5,dq6,dq7]; %vettore velocità
ddq= [ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,ddq7]; %vettore accellerazione
% massa
syms m1 m2 m3 m4 m5 m6 m7 positive % masse dei link

m=[m1,m2,m3,m4,m5,m6,m7]; %vettore delle masse

% inertia tensor
syms I1xx I1xy I1xz I1yy I1yz I1zz real

I1=[ I1xx I1xy I1xz; I1xy I1yy I1yz; I1xz I1yz I1zz];

syms I2xx I2xy I2xz I2yy I2yz I2zz real

I2=[ I2xx I2xy I2xz; I2xy I2yy I2yz; I2xz I2yz I2zz];

syms I3xx I3xy I3xz I3yy I3yz I3zz real

I3=[ I3xx I3xy I3xz; I3xy I3yy I3yz; I3xz I3yz I3zz];

syms I4xx I4xy I4xz I4yy I4yz I4zz real

I4=[ I4xx I4xy I4xz; I4xy I4yy I4yz; I4xz I4yz I4zz];

syms I5xx I5xy I5xz I5yy I5yz I5zz real

I5=[ I5xx I5xy I5xz; I5xy I5yy I5yz; I5xz I5yz I5zz];

syms I6xx I6xy I6xz I6yy I6yz I6zz real

I6=[ I6xx I6xy I6xz; I6xy I6yy I6yz; I6xz I6yz I6zz];

syms I7xx I7xy I7xz I7yy I7yz I7zz real 
 
I7=[ I7xx I7xy I7xz; I7xy I7yy I7yz; I7xz I7yz I7zz];


I=[I1,I2,I3,I4,I5,I6,I7]; % matrice delle inerzie 

% incognite centro di massa
syms c1x c1y c1z c2x c2y c2z c3x c3y c3z c4x c4y c4z c5x c5y c5z c6x c6y c6z c7x c7y c7z real; 
com=[[c1x c1y c1z]',[ c2x c2y c2z]', [c3x c3y c3z ]',[c4x c4y c4z ]', [c5x c5y c5z]', [c6x c6y c6z]',[c7x c7y c7z]']; 
% gravity
%g=9.81;
syms g positive
g0=[0; 0; g]; %vettore della gravità classica(z0 è impostata verticale)


% creazione parametri J 
% definiti come :
% Jxx=Ixx+m*cy^2+m*cz^2
% Jzz=Izz+m*cy^2+m*cx^2
% Jyy=Iyy +m*cx^2+m*cz^2
% Jxy= Ixy -m*cx*cy
% Jxz= Ixz -m*cx*cz
% Jyz= Iyz -m*cy*cz

syms J1xx J1xy J1xz J1yy J1yz J1zz J2xx J2xy J2xz J2yy J2yz J2zz J3xx J3xy J3xz J3yy J3yz J3zz J4xx J4xy J4xz J4yy J4yz J4zz J5xx J5xy J5xz J5yy J5yz J5zz J6xx J6xy J6xz J6yy J6yz J6zz J7xx J7xy J7xz J7yy J7yz J7zz real 

Jxx=[J1xx,J2xx,J3xx,J4xx,J5xx,J6xx,J7xx];
Jyy=[J1yy,J2yy,J3yy,J4yy,J5yy,J6yy,J7yy];
Jzz=[J1zz,J2zz,J3zz,J4zz,J5zz,J6zz,J7zz];
Jxy=[J1xy,J2xy,J3xy,J4xy,J5xy,J6xy,J7xy];
Jxz=[J1xz,J2xz,J3xz,J4xz,J5xz,J6xz,J7xz];
Jyz=[J1yz,J2yz,J3yz,J4yz,J5yz,J6yz,J7yz];
JJ=[Jxx', Jyy', Jzz' Jxy', Jxz', Jyz'];
% %% crezione parametri mc
% % equivalgono a dire mcx=
% syms mc1x mc1y mc1z mc2x mc2y mc2z mc3x mc3y mc3z mc4x mc4y mc4z mc5x mc5y mc5z mc6x mc6y mc6z mc7x mc7y mc7z real; 
% mc=[[mc1x mc1y mc1z]',[ mc2x mc2y mc2z]', [mc3x mc3y mc3z ]',[mc4x mc4y mc4z ]', [mc5x mc5y mc5z]', [mc6x mc6y mc6z]',[mc7x mc7y mc7z]']; 
% %% parametri di attrito modello coulomb+viscoso
syms f01 f02 f03 f04 f05 f06 f06 f07 real
f0=[f01 f02 f03 f04 f05 f06 f06 f07];

syms fc1 fc2 fc3 fc4 fc5 fc6 fc7 real
fc=[fc1 fc2 fc3 fc4 fc5 fc6 fc7];

syms fv1 fv2 fv3 fv4 fv5 fv6 fv7 real
fv=[fv1 fv2 fv3 fv4 fv5 fv6 fv7];
%% identification 

[A,T,pl,J]=param_id(a,d,alphad,theta,com,DH,N,Tool,tipogiunto); %funzione che calcola i dati indicati

% A(:,4*i-3:4*i)        trasformazione omogenea tra giunti consecutivi: tra il giunto i-1 e i
% T(:,4*i-3:4*i)        trasormazione omogena tra il giunto 0 e il giunto i: 
% pl(:,i)               centro di massa del link i pl
% J(:,N*i-(N-1):N*i)    Jacobiano del link i

%%  matrice Inerzia 

B=calcolo_Inerzia1(m,J,I,T,N); %calcolo matrice di inerzia
[B2,B3]=simplify_inertia(B,JJ,N,com,I,a,alphad,theta,d,m,tipogiunto);%semplifica la matrice di inerzia
%B2 vettore della matrice superiore di B semplificata
%B3 matrice semplificata di B

%% Gravità

[G,G2]=calcoloCoppiaGrav1(g0,J,m,N,com); %calcolo matrice di gravità
% G è il vettore di gravità
% G2 è il vettore semplificato di G

%% coriolis 
[C,C2]=calcolo_C(B3,N,JJ,com,a,alphad,theta,d,m,tipogiunto,dq); %calcolo matrice di coriolis
% C è la matrice di coriolis non semplificata
% C2 è il vettore C2=C*dq semplificato

%% Calcolo regressore per il calcolo dei parametri per S (dai valori di inerzia e gravità)

[Yr,PIr,Yrs,PIrs]=calcolo_regr_e_PI_reverse(B2,G2,JJ,N,com,a,alphad,theta,d,m,tipogiunto);
%B2 valori della matrice superiore di inerzia messi in vettore
% B3 matrice B semplificata
% G2 valori di G semplificati
% 
%% regressore classico con l'attrito modello 
[Yc,Pic,Ycs,PIcs]=regressore_classico_friction(JJ,N,B3,ddq,C2,G2,m,com,tipogiunto,theta,alphad,a,d,fc, fv ,f0 ,dq);
%% creazione funzioni per i regressori
matlabFunction(Yrs,'File','Regressore_reverse_semplificato.m','Vars',{theta,g});
matlabFunction(Ycs, 'file', 'Regressore_classico_semplificato.m', 'Vars', {ddq, dq, theta, g});
