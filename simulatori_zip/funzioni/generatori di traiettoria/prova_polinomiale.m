%funzione che crea una traiettoria che impone pari a 0 le velocità e le
%accelerazioni all'inizio e alla fine.
%se le velocità o le accellerazioni sono troppo grandi ritorna un errore
function q=prova_polinomiale(q0 ,q_fin ,run_time,qmax,qmin,dqmax,ddqmax, u)
n=7;
t=u;
q_temp=zeros(n,1);
dq=zeros(n,1);
ddq=zeros(n,1);
a0=q0;
e=q_fin-q0;
a3=zeros(n,1);
a4=zeros(n,1);
a5=zeros(n,1);
if t>run_time
    q_temp=q_fin;
else
for i=1:n
    a3(i)=10*e(i)/(run_time^3);
    a4(i)=-15*e(i)/(run_time^4);
    a5(i)=6*e(i)/(run_time^5);
    q_temp(i)=a5(i)*(t^5)+a4(i)*(t^4)+a3(i)*(t^3)+q0(i);
    dq(i)=5*a5(i)*(t^4)+4*a4(i)*(t^3)+3*a3(i)*(t^2);
    ddq(i)=20*a5(i)*(t^3)+12*a4(i)*(t^2)+6*a3(i)*t;
    if abs(ddq(i))>ddqmax(i) || abs(dq(i))>dqmax(i)
        error('provapoinomiale','velocità o acellerazione troppo elevata');
    end
    if q_temp(i)>qmax(i) || q_temp(i)<qmin(i)
         error('provapoinomiale','posizione fuori dai limiti');
    end
    if q_fin(i)>qmax(i) || q_fin(i)<qmin(i)
     error('provapoinomiale','posizione di arrivo fuori dai limiti');
    end

    
end
q=q_temp;
end
end