% funzione che semplifica B. B2 è il vettore dei valori della matrice superiore, B3 è la matrice B semplificata
% è del tipo [B2,B3]=simplify_inertia(B,JJ,N,com,I,a,alphad,theta,d,m,tipogiunto)
function [B2,B3]=simplify_inertia(B,JJ,N,com,I,a,alphad,theta,d,m,tipogiunto)

% ricreazione dei parametri J dei centri di inerzia 
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

%% Calcolo Derivate e Semplificazione di B
    B1=triu(B); % essendo B simmetrica uso solo i valori della matrice superiore per i nostri scopi
    B1=reshape(B1,1,N*N);
    B1=nonzeros(B1); %ho solo i valori della matrice superiore di B in riga
    % calcolo le derivate di B. i vari debug son utili per capire lo stato
    % della funzione. I valori I son già semplificati sul momento.
    bm=sym(zeros(N,V));
    bmx=sym(zeros(N,V));
    bmy=sym(zeros(N,V));
    bmz=sym(zeros(N,V));
    bmxx=sym(zeros(N,V));
    bmyy=sym(zeros(N,V));
    bmzz=sym(zeros(N,V));
    bmxy=sym(zeros(N,V));
    bmyz=sym(zeros(N,V));
    bmxz=sym(zeros(N,V));
    bmcx=sym(zeros(N,V));
    bmcy=sym(zeros(N,V));
    bmcz=sym(zeros(N,V));
    BM=sym(zeros(N,V));
    Ixx=sym(zeros(N,V));
    Ixy=sym(zeros(N,V));
    Ixz=sym(zeros(N,V));
    Iyy=sym(zeros(N,V));
    Iyz=sym(zeros(N,V));;
    Izz=sym(zeros(N,V));
    n_c=20;
    tot_int=V*n_c;
    for i=1:N;
                cunt=0;
        bm(i,1:V)=diff(B1,m(i)); %derivata su mi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        bmx(i,1:V)=diff(bm(i,1:V),com(1,i)); %derivata su mi e cxi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmy(i,1:V)=diff(bm(i,1:V),com(2,i)); %derivata su mi e cyi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmz(i,1:V)=diff(bm(i,1:V),com(3,i)); %derivata su mi e czi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmxx(i,1:V)=diff(bmx(i,1:V),com(1,i)); %moltiplicatore di mi*cxi^2
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmxy(i,1:V)=diff(bmx(i,1:V),com(2,i)); %moltiplicatore di mi*cxi*cyi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmxz(i,1:V)=diff(bmx(i,1:V),com(3,i)); %moltiplicatore di mi*cxi*czi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmyz(i,1:V)=diff(bmy(i,1:V),com(3,i)); %moltiplicatore di mi*cyi*czi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmyy(i,1:V)=diff(bmy(i,1:V),com(2,i)); %moltiplicatore di mi*cyi^2
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmzz(i,1:V)=diff(bmz(i,1:V),com(3,i)); %moltiplicatore di mi*czi^2
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmcx(i,1:V)=bmx(i,1:V)-bmxx(i,1:V)*com(1,i)-bmxy(i,1:V)*com(2,i)-bmxz(i,1:V)*com(3,i); %moltiplicatore di mi*cxi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmcy(i,1:V)=bmy(i,1:V)-bmyy(i,1:V)*com(2,i)-bmxy(i,1:V)*com(1,i)-bmyz(i,1:V)*com(3,i); %moltiplicatore di mi*cy1
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        bmcz(i,1:V)=bmz(i,1:V)-bmzz(i,1:V)*com(3,i)-bmxz(i,1:V)*com(1,i)-bmyz(i,1:V)*com(2,i); %moltiplicatore di mi*czi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        BM(i,1:V)=bm(i,1:V)-bmcx(i,1:V)*com(1,i)-bmcy(i,1:V)*com(2,i)-bmcz(i,1:V)*com(3,i)-bmxx(i,1:V)/2*com(1,i)^2- bmxy(i,1:V)*com(1,i)*com(2,i)- bmxz(i,1:V)*com(1,i)*com(3,i)- bmyy(i,1:V)/2*com(2,i)^2- bmzz(i,1:V)*com(3,i)^2/2-bmyz(i,1:V)*com(2,i)*com(3,i); %moltiplicatore di mi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        Ixx(i,1:V)=simplify(expand(diff(B1,I(1,3*i-2)))); %moltiplicatore per Ixxi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        Iyy(i,1:V)=simplify(expand(diff(B1,I(2,3*i-1)))); %moltiplicatore per Iyyi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        Izz(i,1:V)=simplify(expand(diff(B1,I(3,3*i)))); %moltiplicatore di Izzi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        Ixy(i,1:V)=simplify(expand(diff(B1,I(2,3*i-2)))); %moltiplicatore di Ixyi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        Ixz(i,1:V)=simplify(expand(diff(B1,I(3,3*i-2)))); %moltiplicatore di Ixzi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        Iyz(i,1:V)=simplify(expand(diff(B1,I(3,3*i-1)))); %moltiplicatore di Iyzi
                avanzamento_diff=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

    end
%semplifico anche i Valori rimanenti per poter ottenere una matrice B
%semplificata e anche le derivate semplificate. è fatto punto per punto per
%permettere una maggior controllo nel debug

    n_c=10;
    tot_int=N*V*n_c;
    
    for i=1:N
        for j=1:V
            cunt=0;
            BM(i,j)=simplify(expand(BM(i,j)));
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
                
            bmcx(i,j)=simplify(expand(bmcx(i,j)));
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
            
            bmcy(i,j)=simplify(expand(bmcy(i,j)));
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
        
            bmcz(i,j)=simplify(expand(bmcz(i,j)));
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
       
            bmxx(i,j)=simplify(expand(bmxx(i,j)));
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
        
            bmxy(i,j)=simplify(expand(bmxy(i,j)));
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;

            bmxz(i,j)=simplify(expand(bmxz(i,j)));
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
        
             bmyy(i,j)=simplify(expand(bmyy(i,j)));
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
        
             bmyz(i,j)=simplify(expand(bmyz(i,j)));
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
                
             bmzz(i,j)=simplify(expand(bmzz(i,j))); 
                avanzamento_simp=[100*(V*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
        end
    end
%% creazione vettore matrice superiore di inerzia semplificata    
B2=sym(zeros(V,1)); %B2 sarà la Matrice B semplificata.
        
for i=1:N
    B2(:)=B2(:)+BM(i,:)'*m(i)+bmcx(i,:)'*m(i)*com(1,i)+bmcy(i,:)'*m(i)*com(2,i)+bmcz(i,:)'*m(i)*com(3,i);
        
    B2(:)=B2(:)+Ixx(i,:)'*Jxx(i)+Ixy(i,:)'*Jxy(i)+Ixz(i,:)'*Jxz(i)+Iyy(i,:)'*Jyy(i)+Izz(i,:)'*Jzz(i)+Iyz(i,:)'*Jyz(i);
       
end

% il concetto è semplice: ogni valore di p=(p1, p2, p3) moltiplica dei seni
% e coseni degli angoli non semplificabili per un altro valore di p.
% la semplificazione viene fatta quindi per ogni valore di p (utilizzando le derivazioni)
% e poi vengono moltiplicate le funzioni semplificate per ogni valore di p
% corrispondente e sommati insieme. questo creerà il vettore B2
% semplificato con un calcolo computazionale minore
%% ricreazione matrice quadrata di inerzia.
B3=lineartomatrixI(B2,N);
