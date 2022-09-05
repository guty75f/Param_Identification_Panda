% crea la funzione di coriolis vettore
% necessita di JJ, N, C2, m, com, tipogiunto, theta, dq
function i=funzione_coriolis(JJ,N,C2,m,com,tipogiunto,theta,dq,PI)
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
V=N; %valore dei dati della matrice B utilizzabili
%% Calcolo Derivate utili di C per la costruzione della funzione di coriolis
    % calcolo le derivate di C. i vari debug son utili per capire lo stato
    % della funzione. tutti i valori son già stati semplificati
    % precedentemente
    c3m=sym(zeros(N,V));
    c3mcx=sym(zeros(N,V));
    c3mcy=sym(zeros(N,V));
    c3mcz=sym(zeros(N,V));
    c3M=sym(zeros(N,V));
    c3Ixx=sym(zeros(N,V));
    c3Ixy=sym(zeros(N,V));
    c3Ixz=sym(zeros(N,V));
    c3Iyy=sym(zeros(N,V));
    c3Iyz=sym(zeros(N,V));;
    c3Izz=sym(zeros(N,V));
    n_c=11; %numero di derivate da calcolare, serve per il debug
    tot_int=N*n_c;
    for i=1:N;
                cunt=0;
        c3m(i,:)=diff(C2,m(i)); %derivata su mi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        c3mcx(i,:)=diff(c3m(i,1:V),com(1,i)); %derivata su mi e cxi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        c3mcy(i,:)=diff(c3m(i,1:V),com(2,i)); %derivata su mi e cyi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        c3mcz(i,:)=diff(c3m(i,1:V),com(3,i)); %derivata su mi e czi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        c3M(i,:)=c3m(i,1:V)-c3mcx(i,1:V)*com(1,i)-c3mcy(i,1:V)*com(2,i)-c3mcz(i,1:V)*com(3,i); %moltiplicatore di mi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        c3Ixx(i,:)=diff(C2,Jxx(i)); %moltiplicatore per Ixxi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        c3Iyy(i,:)=diff(C2,Jyy(i)); %moltiplicatore per Iyyi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        c3Izz(i,:)=diff(C2,Jzz(i)); %moltiplicatore di Izzi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        c3Ixy(i,:)=diff(C2,Jxy(i)); %moltiplicatore di Ixyi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        c3Ixz(i,:)=diff(C2,Jxz(i)); %moltiplicatore di Ixzi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        c3Iyz(i,:)=diff(C2,Jyz(i)); %moltiplicatore di Iyzi
                avanzamento_creazione_coriolis=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

    end

%% 
[JJxx, JJxy, JJxz, JJyz, JJzz, mcx, mcy]=risultati_dai_test(N,PI); 
C3=sym(zeros(V,1));
%% Calcolo di PI e Y
for i=1:N
    C3(:)=C3(:)+mcx(i)*c3mcx(i,:)'+mcy(i)*c3mcy(i,:)';
        
    C3(:)=C3(:)+c3Ixx(i,:)'*JJxx(i)+c3Ixy(i,:)'*JJxy(i)+c3Ixz(i,:)'*JJxz(i)+c3Izz(i,:)'*JJzz(i)+c3Iyz(i,:)'*JJyz(i);
end
i=1;
matlabFunction(C3, 'File', 'Vettore_coriolis.m', 'Vars', {dq, theta});

end
