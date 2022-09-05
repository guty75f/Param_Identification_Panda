% funzione che calcola la matrice di inerzia non approssimata
% è del tipo:B=calcolo_Inerzia1(m,J,I,T,N)
function B=calcolo_Inerzia1(m,J,I,T,N)
B=sym(zeros(7,7));
for i=1:N
    Jp=J(1:3,(N*i-6):N*i);
    Jo=J(4:6,(N*i-6):N*i);
    R = T(1:3,4*i-3:4*i-1);
    Ji=R*I(:,(3*i-2):3*i)*R';
    %Ji=I(:,(3*i-2):3*i);
    %
    B=B+(m(i)*Jp'*Jp+Jo'*Ji*Jo);
   

end

end
