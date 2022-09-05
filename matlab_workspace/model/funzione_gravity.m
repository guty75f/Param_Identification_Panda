%funzione che crea la funzione gravità da G2
%Richiede in input theta, tipogiunto, G2, N, com e m;
function i=funzione_gravity(theta, tipogiunto, G2, N, com, m,g,PI)
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

gm=sym(zeros(N,N));
gmx=sym(zeros(N,N));
gmy=sym(zeros(N,N));
gmz=sym(zeros(N,N));
M=sym(zeros(N,N));
n_c=6;
tot_int=n_c*N;
    for  i=1:N
            cunt=0;
        gm(i,:)=(diff(G2,m(i))); %derivata per mi
                avanzamento_gravity_function=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        gmx(i,:)=diff(gm(i,:),com(1,i)); % coefficenti di mi*cxi
              avanzamento_gravity_function=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        gmy(i,:)=diff(gm(i,:),com(2,i)); % coefficenti di mi*cyi
                avanzamento_gravity_function=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        gmz(i,:)=diff(gm(i,:),com(3,i)); % coefficenti di mi*czi
                avanzamento_gravity_function=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        M(i,:)=gm(i,:)-gmx(i,:)*com(1,i)-gmy(i,:)*com(2,i)-gmz(i,:)*com(3,i); % coefficenti di mi
                avanzamento_gravity_function=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
end
[JJxx, JJxy, JJxz, JJyz, JJzz, mcx, mcy]=risultati_dai_test(N,PI);
G3=sym(zeros(7,1));
%% Calcolo di PI e Y
for i=1:N
    G3(:)=G3(:)+mcx(i)*gmx(i,:)'+mcy(i)*gmy(i,:)';       
end
matlabFunction(G3, 'File', 'Vettore_Gravita.m', 'Vars', {theta,g});
end
