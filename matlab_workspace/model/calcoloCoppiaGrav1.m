% calcola il vettore di gravità G e quello semplificato G2
%è del tipo: [G,G2]=calcoloCoppiaGrav1(g0,J,m,N,com)
function [G,G2]=calcoloCoppiaGrav1(g0,J,m,N,com)

G=sym(zeros(1,N));


for i=1:N
    g(i,:)=m(i)*g0'*J(1:3,N*i-N+1:N*i);
    G=G+g(i,:);
end

%% semplificazione di G
G2=sym(zeros(1,N))
gm=sym(zeros(N,N));
gmx=sym(zeros(N,N));
gmy=sym(zeros(N,N));
gmz=sym(zeros(N,N));
M=sym(zeros(N,N));
n_c=6;
tot_int=n_c*N;
    for  i=1:N
            cunt=0;
        gm(i,:)=(diff(G,m(i))); %derivata per mi
                avanzamento_semplificazione_G=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        gmx(i,:)=simplify(expand(diff(gm(i,:),com(1,i)))); % coefficenti di mi*cxi
              avanzamento_semplificazione_G=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        gmy(i,:)=simplify(expand(diff(gm(i,:),com(2,i)))); % coefficenti di mi*cyi
                avanzamento_semplificazione_G=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        gmz(i,:)=simplify(expand(diff(gm(i,:),com(3,i)))); % coefficenti di mi*czi
                avanzamento_semplificazione_G=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        M(i,:)=simplify(expand(gm(i,:)-gmx(i,:)*com(1,i)-gmy(i,:)*com(2,i)-gmz(i,:)*com(3,i))); % coefficenti di mi
                avanzamento_semplificazione_G=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        G2(1,:)=G2(1,:)+(gmx(i,:)*com(1,i)+gmy(i,:)*com(2,i)+gmz(i,:)*com(3,i)+M(i,:))*m(i); %G semplificato
                avanzamento_semplificazione_G=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
    end
G2=G2';
end
