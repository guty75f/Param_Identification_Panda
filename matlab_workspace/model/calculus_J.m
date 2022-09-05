% 
% del tipo J=calculus_J(n,T,com,N,tipogiunto)
function J=calculus_J(n,T,com,N,tipogiunto)
% n giunto di cui calcolare il jacobiano
% T trasformazioni omogenee
% com(:,i) centro di massa del link i
% N numero di giunti
% tipogiunto: 1 se rotazionale o 0 se prismatico

    J = sym(zeros(6, N)); %inizializzo J come matrice 6xN
    pippo=T(:,n*4-3:n*4)*[com(:,n);1]; %calcolo per trovare il centro di massa rispetto alla terna 0
    Pl=pippo(1:3); %elimino il 4* numero in quanto è sempre 1 e non fa parte del centro di massa rispetto alla terna 0
    clear pippo;

    for j=n:-1:1
        z=T(1:3,j*4-1); %rotazione su z del giunto j
        if (tipogiunto(j)==1) %se il giunto j è rotazionale 
       
            p=T(1:3,j*4); %posizione del link j
        
            d=simplify(expand(cross(z,Pl-p))); 
            % d che sono i primi 3 valori di J del giunto n rispetto al
            % giunto j è calcolato come come cross product tra la posizione del centro di
            % massa del giunto n meno la posizione del giunto j e la
            % rotazione su z nel caso rotazionale
        
            J(:,j) = [d; z]; % d già definito, gli altri tre valori sono la rotazione nel giunto j su z rispetto alla terna 0
        end
        if (tipogiunto(j)==0) % se il giunto j è prismatico
            J(:,j)=[z; 0; 0; 0]; % i primi tre valori sono la rotazione del giunto j su z rispetto alla terna 0, gli altri valori sono 0
        end
        avanzamento_calcolo_J=[n,j] % debug
    end
end
