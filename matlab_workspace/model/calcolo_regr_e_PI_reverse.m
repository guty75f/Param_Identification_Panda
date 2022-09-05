% crea il regressore e i PI per il calcolo con reverse_eng. Yrs e PIrs sono senza gli zeri
% è del tipo [Yr,PIr,Yrs,PIrs]=calcolo_regr_e_PI_reverse(B2,G2,JJ,N,com,a,alphad,theta,d,m,tipogiunto)
function [Yr,PIr,Yrs,PIrs]=calcolo_regr_e_PI_reverse(B2,G2,JJ,N,com,a,alphad,theta,d,m,tipogiunto)
%% Ridefinizione JJ in Jxx ed ecc

Jxx=JJ(:,1)';
Jyy=JJ(:,2)';
Jzz=JJ(:,3)';
Jxy=JJ(:,4)';
Jxz=JJ(:,5)';
Jyz=JJ(:,6)';
c=sym(zeros(N,1));
s=sym(zeros(N,1));

for i=1:N
    if (tipogiunto(i)==1) % se rotazionale creo i coseni e i seni come incognite
   
        c(i)=cos(theta(i));
    
        s(i)=sin(theta(i));
    end
    if (tipogiunto(i)==0) % se prismatico i coseni e i seni son calcolati in deg per evitare che sin(0)!=0 e sin(pi/2)!=1
        thetad=rad2deg(theta(i));
        c(i)=cosd(thetad);
        s(i)=sind(thetad);
    end

end
V=N*(N+1)/2; %valore dei dati della matrice B utilizzabili
Vs=V+N;% lunghezza del vettore s=(B2;G2);
%% creazione matrice s=(B2;G2)
S=sym(zeros(1,Vs));
S(1,1:V)=B2;
S(1,V+1:Vs)=G2;
%% Calcolo Derivate utili di S per la costruzione del regressore di S
    % calcolo le derivate di S. i vari debug son utili per capire lo stato
    % della funzione. tutti i valori son già stati semplificati
    % precedentemente
    sm=sym(zeros(N,Vs));
    smcx=sym(zeros(N,Vs));
    smcy=sym(zeros(N,Vs));
    smcz=sym(zeros(N,Vs));
    SM=sym(zeros(N,Vs));
    sIxx=sym(zeros(N,Vs));
    sIxy=sym(zeros(N,Vs));
    sIxz=sym(zeros(N,Vs));
    sIyy=sym(zeros(N,Vs));
    sIyz=sym(zeros(N,Vs));;
    sIzz=sym(zeros(N,Vs));
    n_c=11; %numero di derivate da calcolare, serve per il debug
    tot_int=N*n_c;
    for i=1:N;
                cunt=0;
        sm(i,:)=diff(S,m(i)); %derivata su mi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        smcx(i,:)=diff(sm(i,1:Vs),com(1,i)); %derivata su mi e cxi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        smcy(i,:)=diff(sm(i,1:Vs),com(2,i)); %derivata su mi e cyi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        smcz(i,:)=diff(sm(i,1:Vs),com(3,i)); %derivata su mi e czi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        SM(i,:)=sm(i,1:Vs)-smcx(i,1:Vs)*com(1,i)-smcy(i,1:Vs)*com(2,i)-smcz(i,1:Vs)*com(3,i); %moltiplicatore di mi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        sIxx(i,:)=diff(S,Jxx(i)); %moltiplicatore per Ixxi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        sIyy(i,:)=diff(S,Jyy(i)); %moltiplicatore per Iyyi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        sIzz(i,:)=diff(S,Jzz(i)); %moltiplicatore di Izzi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        sIxy(i,:)=diff(S,Jxy(i)); %moltiplicatore di Ixyi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        sIxz(i,:)=diff(S,Jxz(i)); %moltiplicatore di Ixzi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        sIyz(i,:)=diff(S,Jyz(i)); %moltiplicatore di Iyzi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

    end

%% Accorpazione dei parametri secondo le regole di??

mz=sym(zeros(1,N)); %cx*m
my=sym(zeros(1,N)); %cy*m
mx=sym(zeros(1,N)); %cmz*m
mm=m;
f=alphad;
for i=1:N
 mz(i)=mm(i)*com(3,i);
    my(i)=mm(i)*com(2,i);
    mx(i)=mm(i)*com(1,i);
end

for i=N:-1:2
   
    Jxx(i)=Jxx(i)-Jyy(i);
    Jxx(i-1)=Jxx(i-1)+Jyy(i)+2*d(i)*mz(i)+d(i)^2*mm(i);
    Jxy(i-1)=Jxy(i-1)+a(i)*sind(f(i))*mz(i)+a(i)*d(i)*sind(f(i))*mm(i);
    Jxz(i-1)=Jxz(i-1)-a(i)*cosd(f(i))*mz(i)-a(i)*d(i)*cosd(f(i))*mm(i);
    Jyy(i-1)=Jyy(i-1)+cosd(f(i))^2*Jyy(i)+2*d(i)*cosd(f(i))^2*mz(i)+mm(i)*(a(i)^2+d(i)^2*cosd(f(i))^2);
    Jyz(i-1)=Jyz(i-1)+cosd(f(i))*sind(f(i))*Jyy(i)+2*d(i)*cosd(f(i))*sind(f(i))*mz(i)+d(i)^2*cosd(f(i))*sind(f(i))*mm(i);
    Jzz(i-1)=Jzz(i-1)+sind(f(i))^2*Jyy(i)+2*d(i)*sind(f(i))^2*mz(i)+mm(i)*(a(i)^2+d(i)^2*sind(f(i))^2);
    mx(i-1)=mx(i-1)+a(i)*mm(i);
    my(i-1)=my(i-1)-sind(f(i))*mz(i)-d(i)*sind(f(i))*mm(i);
    mz(i-1)=mz(i-1)+cosd(f(i))*mz(i)+d(i)*cosd(f(i))*mm(i);
    mm(i-1)=mm(i-1)+mm(i);
end

%% Calcolo di PI e Y
PI_Ixx=[sIxx];
PI_Iyz=[sIyz];
PI_Izz=[sIzz];
PI_Ixz=[sIxz];
PI_Ixy=[sIxy];
PI_mcx=[smcx];
PI_mcy=[smcy];
PIr=sym(zeros(1,N*N));
Yr=sym(zeros(Vs,N*N));
for i=1:N
        PIr(1,(N*i-(N-1)):N*i)=[Jxx(i), Jxy(i), Jxz(i), Jyz(i), Jzz(i), mx(i), my(i)];
         Yr(:,(N*i-(N-1)):N*i)=[PI_Ixx(i,:)',PI_Ixy(i,:)', PI_Ixz(i,:)', PI_Iyz(i,:)', PI_Izz(i,:)',PI_mcx(i,:)',PI_mcy(i,:)']; 
end
% tutti i valori Iyy, mz e mm non sono linearmente indipendenti e quindi
% vengono esclusi di partenza dal regressore.
[Yrs,PIrs]=delete_zero(Yr,PIr,N); %vengono cancellati i parametri nulli

end
