% la funzione calcola la matrice di coriolis C (non semplificata) e il vettore Cq=C*dq semplificato
% è della forma  [C,C2]=calcolo_C(B,N,JJ,com,a,alphad,theta,d,m,tipogiunto,dq)
function [C,C2]=calcolo_C(B,N,JJ,com,a,alphad,theta,d,m,tipogiunto,dq)

%% inizializzazione
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
%% calcolo C
C=sym(zeros(N,N));
 for i=1:N
    for j=1:N
 for k=1:N
    C(i,j)= C(i,j)+(0.5*(diff(B(i,j),theta(k))+diff(B(i,k),theta(j))-diff(B(j,k),theta(i))))*dq(k);
    yaaiiii=[i j k]
 end
    end
 end


%% calcolo CQ
Cq=C*dq';

%% Calcolo Derivate e Semplificazione di Cq
    cm=sym(zeros(N,N));
    cmx=sym(zeros(N,N));
    cmy=sym(zeros(N,N));
    cmz=sym(zeros(N,N));
    cM=sym(zeros(N,N));
    cIxx=sym(zeros(N,N));
    cIxy=sym(zeros(N,N));
    cIxz=sym(zeros(N,N));
    cIyy=sym(zeros(N,N));
    cIyz=sym(zeros(N,N));;
    cIzz=sym(zeros(N,N));
    n_c=11;
    tot_int=N*n_c;
    for i=1:N;
        cunt=0;
        cm(i,1:N)=diff(Cq,m(i)); %derivata su mi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        cmcx(i,1:N)=diff(cm(i,1:N),com(1,i)); %derivata su mi e cxi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        cmcy(i,1:N)=diff(cm(i,1:N),com(2,i)); %derivata su mi e cyi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        cmcz(i,1:N)=diff(cm(i,1:N),com(3,i)); %derivata su mi e czi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        CM(i,1:N)=cm(i,1:N)-cmcx(i,1:N)*com(1,i)-cmcy(i,1:N)*com(2,i)-cmcz(i,1:N)*com(3,i); %moltiplicatore di mi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        cIxx(i,1:N)=diff(Cq,Jxx(i)); %moltiplicatore per Jxxi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        cIyy(i,1:N)=diff(Cq,Jyy(i)); %moltiplicatore per Jyyi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        cIzz(i,1:N)=diff(Cq,Jzz(i)); %moltiplicatore di Jzzi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        cIxy(i,1:N)=diff(Cq,Jxy(i)); %moltiplicatore di Jxyi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        cIxz(i,1:N)=diff(Cq,Jxz(i)); %moltiplicatore di Jxzi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        cIyz(i,1:N)=diff(Cq,Jyz(i)); %moltiplicatore di Iyzi
                avanzamento_diff_C=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
    end
%semplifico anche i valori rimanenti per poter ottenere una matrice CQ
%semplificata e anche le derivate semplificate. è fatto punto per punto per
%permettere una maggior controllo nel debug
    n_c=10;
    tot_int=N*N*n_c;
    for i=1:N
        for j=1:N
            cunt=0;
            CM(i,j)=simplify_this_C(CM(i,j),dq,N);
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
        
            cmcx(i,j)=simplify_this_C(cmcx(i,j),dq,N);
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
        
            cmcy(i,j)=simplify_this_C(cmcy(i,j),dq,N);
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
        
            cmcz(i,j)=simplify_this_C(cmcz(i,j),dq,N);
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
       
            cIxx(i,j)=simplify_this_C(cIxx(i,j),dq,N); %moltiplicatore per Ixxi
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
             cIyy(i,j)=simplify_this_C(cIyy(i,j),dq,N); %moltiplicatore per Iyyi
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
             cIzz(i,j)=simplify_this_C(cIzz(i,j),dq,N); %moltiplicatore di Izzi
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
             cIxy(i,j)=simplify_this_C(cIxy(i,j),dq,N); %moltiplicatore di Ixyi
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
             cIxz(i,j)=simplify_this_C(cIxz(i,j),dq,N); %moltiplicatore di Ixzi
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
             cIyz(i,j)=simplify_this_C(cIyz(i,j),dq,N); %moltiplicatore di Iyzi
                avanzamento_simp_C=[100*(N*n_c*(i-1)+n_c*(j-1)+cunt)/tot_int]
                cunt=cunt+1;
        end
    end
%% costruzione matrice Cq semplificata    
C2=sym(zeros(N,1)); %C2 sarà la Matrice Cq=C*dq semplificata.
for i=1:N
    C2(:)=C2(:)+CM(i,:)'*m(i)+cmcx(i,:)'*m(i)*com(1,i)+cmcy(i,:)'*m(i)*com(2,i)+cmcz(i,:)'*m(i)*com(3,i);
    C2(:)=C2(:)+cIxx(i,:)'*Jxx(i)+cIxy(i,:)'*Jxy(i)+cIxz(i,:)'*Jxz(i)+cIyy(i,:)'*Jyy(i)+cIzz(i,:)'*Jzz(i)+cIyz(i,:)'*Jyz(i);
end
end
