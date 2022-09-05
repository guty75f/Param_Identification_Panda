% ritorna la trasformazione omogenea tra il giunto i e la terna 0;
% del tipo T=calcolo_Transformation_matrix(A,N)
function T=calcolo_Transformation_matrix(A,N)
% ritorna la trasformazione omogenea tra il giunto i e la terna 0;

T_temp(:,1:4)=A(:,1:4); % sul giunto 1 la matrice T corrisponde alla matrice A
  avanzamento_calcolo_T=1 %debug

for i=2:N
  T_temp(:,i*4-3:4*i)=T_temp(:,((i-1)*4-3): 4*(i-1))*A(:,4*i-3:4*i); % T0_i=T0_i-1*Ai-1_i
  avanzamento_calcolo_T=i %debug
end
T=T_temp;
end

