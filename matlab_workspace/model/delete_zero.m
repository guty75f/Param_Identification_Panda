% cancella gli zeri dei regressori e dei suoi parametri PI ( tau=Y*PI)
% del tipo [Ys,PI]=delete_zero(Y,PI,N)
function [Yrs,PIs]=delete_zero(Yr,PI,N)
     M=length(Yr(1,:));
     V=length(Yr(:,1));
     j=1;
     Yrs_temp=sym(zeros(V,M));
     PIs_temp=sym(zeros(1,M));
    for i=1:M
        if not(isequal(Yr(:,i), sym(zeros(V,1))))
            Yrs_temp(:,j)=Yr(:,i);
            PIs_temp(1,j)=PI(1,i);
            j=j+1
        end
        i
    end
    Yrs=sym(zeros(V,j-1));
    PIs=sym(zeros(1,j-1));
 Yrs=Yrs_temp(:,1:j-1);
 PIs=PIs_temp(:,1:j-1);
end
            
