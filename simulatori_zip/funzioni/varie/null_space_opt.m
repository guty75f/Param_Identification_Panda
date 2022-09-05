function qp=null_space_opt(robot,u)
    %optimiziation of null_space motion with joint limits gradient 
    %% initialization
    n=length(u);
    qmin=[-2.8973 	-1.7628 	-2.8973 	-3.0718 	-2.8973 	-0.0175 	-2.8973 -pi/4];% max q raggiungibili
    qmax=[ 2.8973 	1.7628 	2.8973 	-0.0698 	2.8973 	3.7525 	2.8973 -pi/4]; %min q raggiungibili
    if n<8 || ( u(1)==0 && u(2)==0 && u(3)==0 && u(4)==0 && u(5)==0 && u(6)==0 ) % per inizializzare
        u= [1, -pi/4, 0, -3 * pi/4, 0, pi/2, pi/4, -pi/4]; % initial position
    end
    qh=(qmin+qmax)/2; %posizione media dei giunti
    %% computation
    J=robot.jacob0(u);  %calcolo il jacobiano con RTB
    Jinv=pseudo_inversa_jacob(robot,u); %calcolo il pseudo inverso del jacobiano con la nostra funzione
    I=eye(n);
    P=I-Jinv*J; %matrice di proiezione
    for i=1:7
        q(i)=-(u(i)-qh(i))/n/(qmax(i)-qmin(i))^2; %calcolo della funzione obbiettivo con l'idea di stare distanti dai limiti
     end
    q(8)=0; %forzato perch+ non si può muovere
    qp=P*q'; %funzione particolare
end

