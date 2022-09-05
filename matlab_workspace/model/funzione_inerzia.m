% la funzione crea la funzione inerzia,
% richiede N, tipogiunto, theta, JJ, B3, com, m
function i=funzione_inerzia(N,tipogiunto,theta,JJ, B2, com, m,PI)

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
%% creazione matrice s=(B2;G2)
%% Calcolo Derivate utili di S per la costruzione del regressore di S
    % calcolo le derivate di S. i vari debug son utili per capire lo stato
    % della funzione. tutti i valori son già stati semplificati
    % precedentemente
    b4m=sym(zeros(N,V));
    b4mcx=sym(zeros(N,V));
    b4mcy=sym(zeros(N,V));
    b4mcz=sym(zeros(N,V));
    B4M=sym(zeros(N,V));
    b4Ixx=sym(zeros(N,V));
    b4Ixy=sym(zeros(N,V));
    b4Ixz=sym(zeros(N,V));
    b4Iyy=sym(zeros(N,V));
    b4Iyz=sym(zeros(N,V));;
    b4Izz=sym(zeros(N,V));
    n_c=11; %numero di derivate da calcolare, serve per il debug
    tot_int=N*n_c;
    for i=1:N;
                cunt=0;
        b4m(i,:)=diff(B2,m(i)); %derivata su mi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        b4mcx(i,:)=diff(b4m(i,1:V),com(1,i)); %derivata su mi e cxi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        b4mcy(i,:)=diff(b4m(i,1:V),com(2,i)); %derivata su mi e cyi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        b4mcz(i,:)=diff(b4m(i,1:V),com(3,i)); %derivata su mi e czi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        B4M(i,:)=b4m(i,1:V)-b4mcx(i,1:V)*com(1,i)-b4mcy(i,1:V)*com(2,i)-b4mcz(i,1:V)*com(3,i); %moltiplicatore di mi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        b4Ixx(i,:)=diff(B2,Jxx(i)); %moltiplicatore per Ixxi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        b4Iyy(i,:)=diff(B2,Jyy(i)); %moltiplicatore per Iyyi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        b4Izz(i,:)=diff(B2,Jzz(i)); %moltiplicatore di Izzi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        b4Ixy(i,:)=diff(B2,Jxy(i)); %moltiplicatore di Ixyi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        b4Ixz(i,:)=diff(B2,Jxz(i)); %moltiplicatore di Ixzi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        b4Iyz(i,:)=diff(B2,Jyz(i)); %moltiplicatore di Iyzi
                avanzamento_diff_regressore_reverse=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

    end

%% Accorpazione dei parametri secondo le regole di??
[JJxx, JJxy, JJxz, JJyz, JJzz, mcx, mcy]=risultati_dai_test(N,PI);
B4=sym(zeros(V,1));
%% Calcolo di PI e Y
for i=1:N
    B4(:)=B4(:)+mcx(i)*b4mcx(i,:)'+mcy(i)*b4mcy(i,:)';
        
    B4(:)=B4(:)+b4Ixx(i,:)'*JJxx(i)+b4Ixy(i,:)'*JJxy(i)+b4Ixz(i,:)'*JJxz(i)+b4Izz(i,:)'*JJzz(i)+b4Iyz(i,:)'*JJyz(i);
end
      

B5=lineartomatrixI(B4,N);
matlabFunction(B5, 'File', 'Matrice_inerzia.m', 'Vars', {theta});
end
