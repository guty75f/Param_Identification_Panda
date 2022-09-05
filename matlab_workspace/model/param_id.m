% la funzione riceve le incognite dal Main e ritorna i valori richiesti:
% A(:,4*i-3:4*i)        trasformazione omogenea tra giunti consecutivi: tra il giunto i-1 e i
% T(:,4*i-3:4*i)        trasormazione omogena tra il giunto 0 e il giunto i: 
% pl(:,i)               centro di massa del link i pl
% J(:,N*i-(N-1):N*i)    Jacobiano del link i
% gli input richiesti sono:
% a vettore dei parametri DH di a
% d vettore dei parametri DH di d
% alphad vettore dei parametri DH di alpha in deg 
% theta vettore dei parametri DH di theta 
% com posizione dei centri di massa rispetto al proprio link--> com(:,i)
%       riferito al giunto i
% DH tipo di convenzione. DH=0 se classica, DH=1 se modificata
% N numero di giunti
% Tool eventuale Tool del cobot
% tipogiunto vettore che indica il tipo di giunto. 0 se prismatico, 1 se
%       rotazionale
function [A,T,pl,J]=param_id(a,d,alphad,theta,com,DH,N,Tool,tipogiunto)
%%  descrizioni vecchie 

% B                     Matrice di Inerzia non semplificata
% G                     Matrice di gravità non semplificata

% a vettore dei parametri DH di a
% d vettore dei parametri DH di d
% alphad vettore dei parametri DH di alpha in deg 
% theta vettore dei parametri DH di theta 
% m vettore delle masse
% g0 vettore di gravità
% com posizione dei centri di massa rispetto al proprio link--> com(:,i)
%       riferito al giunto i
% DH tipo di convenzione. DH=0 se classica, DH=1 se modificata
% N numero di giunti
% Tool eventuale Tool del cobot
% tipogiunto vettore che indica il tipo di giunto. 0 se prismatico, 1 se
%       rotazionale



%% definizione coseni e seni

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

%% Trasformazioni A 
    %Classic DH
	if(DH==0)
     for i=1:N
    
            A(:,4*i-3:4*i)=[
                c(i),-s(i)*cosd(alphad(i)), s(i)*sind(alphad(i)),a(i)*c(i);
        
                s(i),c(i)*cosd(alphad(i)), -c(i)*sind(alphad(i)),a(i)*s(i);
        
                0,   sind(alphad(i)),       cosd(alphad(i)),     d(i);
        
                0,     0,                     0,              1;];
            avanzamento_calcolo_A=i %debug
            
     end
    end

    %Modified DH
    if (DH==1)

        for i=1:N
        
            A(:,4*i-3:4*i)=            [
                c(i),  -s(i),   0, a(i);
             
                s(i)*cosd(alphad(i)),c(i)*cosd(alphad(i)), -sind(alphad(i)),-d(i)*sind(alphad(i));
            
                sind(alphad(i))*s(i),   sind(alphad(i))*c(i),       cosd(alphad(i)),   cosd(alphad(i))*d(i);
           
              0,     0,                     0,              1;];
             avanzamento_calcolo_A=i %debug


        end
    end
%% Trasformazione omogenea rispetto a terna 0

T=calcolo_Transformation_matrix(A,N); %ritorna la matrice 4x(4*N) contenente tutti le trasformazioni omogenee
    
%% posizioni e z
 p0=[0; 0;  0;]; %origine terna 0
 z0=[0;0;1]; % rotazione asse z

 [p,z]=calcolo_pos_z(T,N);

 %%  pos centri di massa

for i=i:N
    h=[com(:,i); 1;];
    pl(:,i)=T(:,4*i-3:4*i)*h; %sfrutto il calcolo per trovare la posizione di un punto rispetto alla terna 0 pl(i)_0=T0_i*[pl(i)_i;i]
      avanzamento_calcolo_Pl=i %debug

end
%% Jacobiano
for i=1:N
   if (DH==0) %DH standard
    J(:,N*i-(N-1):N*i)=calculus_J_old(i,T,A,com,N,Tool,DH,tipogiunto);end %ritorna il jacobiano per ogni link
    if (DH==1) %MDH
    J(:,N*i-(N-1):N*i)=calculus_J(i,T,com,N,tipogiunto);end %ritorna il jacobiano per ogni link
end




end



