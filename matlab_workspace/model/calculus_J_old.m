% calcola il jacobiano 
% del tipo J=calculus_J_old(n,T,A,com,N,tool,DH,tipogiunto)
function J=calculus_J_old(n,T,A,com,N,tool,DH,tipogiunto)
% in breve la funzione calcola il jacobiano rispetto all'end effector e poi
% calcola da quest'ultimo il jacobiano rispetto alla terna 0. vecchio
% programma, usato solo se DH è classica, computazionalmente pesante
    J = sym(zeros(6, N));
    if(n==N)
        U = tool;
    else        
        U=A(:,4*n+1:4*n+4);
    end
    
        U=sym(U);
        UT=U;
         U(1:3,4)=com(:,n);
         
    for j=n:-1:1
        if (DH == 0)
		U = A(:,4*j-3:4*j) * U;
        end
        UT = U;
        delta = UT(3,1:3)';  % nz oz az
	if (tipogiunto(j)==1)
          	 d = [	-UT(1,1)*UT(2,4)+UT(2,1)*UT(1,4)
                	-UT(1,2)*UT(2,4)+UT(2,2)*UT(1,4)
                	-UT(1,3)*UT(2,4)+UT(2,3)*UT(1,4)];
	         J(:,j) = [d; delta];
        end
	if (tipogiunto(j)==0)
		J(:,j)=[delta; 0;0;0];
            

            % modified DH convention
            if (DH==1)
	
	            U = A(:,4*j-3:4*j) * U;
            end
            
    end
    if(n==N)
        T_temp=T(:,N*4-3:N*4)*tool;
    else
        T_temp=T(:,4*n+1:4*n+4);
    end
    R=T_temp(1:3,1:3);
    J0= [R zeros(3,3); zeros(3,3) R] *J;
    J=J0;
  % J= simplify(expand(J0));
end
