%funzione che crea una traiettoria che impone pari a 0 le velocità e le
%accelerazioni all'inizio e alla fine.
%se le velocità o le accellerazioni sono troppo grandi ritorna un errore
function dq=prova_polinomiale_vel (t,q0,q_fin,run_time)
% clear all
% t=9.8;
    %q0= [0, -pi/4, 0, -3 * pi/4, 0, pi/2, pi/4,-pi/4];
   %q_fin=[0   -1.5708   -1.5708    1.5708         0   -1.5708    1.5708 -pi/4];
    %run_time=10;
    dqmax=[2.1750 	2.1750 	2.1750 	2.1750 	2.6100 	2.6100 	2.6100 0]; % max absolute velocity for every joint
    ddqmax=[15 	7.5 	10 	12.5 	15 	20 	20 0]; % max acceleration for every joint
    n=8;
    q_temp=zeros(8,1);
    dq_temp=zeros(8,1);
    ddq=zeros(8,1);
    a0=q0;
    e=q_fin-q0;
    a3=zeros(8,1);
    a4=zeros(8,1);
    a5=zeros(8,1);
    if t>run_time
        q_temp=q_fin;
        dq_temp=zeros(8,1);
    else
    for i=1:n
        a3(i)=10*e(i)/(run_time^3);
        a4(i)=-15*e(i)/(run_time^4);
        a5(i)=6*e(i)/(run_time^5);
        q_temp(i)=a5(i)*(t^5)+a4(i)*(t^4)+a3(i)*(t^3)+q0(i);
        dq_temp(i)=5*a5(i)*(t^4)+4*a4(i)*(t^3)+3*a3(i)*(t^2);
        ddq(i)=20*a5(i)*(t^3)+12*a4(i)*(t^2)+6*a3(i)*t;
        if abs(ddq(i))>ddqmax(i) || abs(dq_temp(i))>dqmax(i)
            error('provapoinomiale','velocità o acellerazione troppo elevata');
        end;
            
    end
    end
   dq=dq_temp;
    
end