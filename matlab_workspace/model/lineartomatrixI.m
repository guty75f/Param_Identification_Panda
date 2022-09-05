% funzione che trasforma il vettore dei valori della matrice superiore di inerzia nella matrice di inerzia simmetrica
% è del tipo: B=lineartomatrixI(A,N) con A il vettore
function B=lineartomatrixI(A,N)
k=0;
for i=1:N
    for j=1:N
        if(i>=j)
            k=k+1;
            B(j,i)=A(k);
            B(i,j)=A(k);
        end
    end
end
end
