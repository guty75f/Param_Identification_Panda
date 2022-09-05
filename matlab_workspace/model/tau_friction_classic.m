%calcola il modello di attrito classico (coulomb+viscoso)
% del tipo tf=tau_friction_classic(dq,f0,fc,fv,N)
function tf=tau_friction_classic(dq,f0,fc,fv,N)
    for i=1:N
        tf(i)=sign(dq(i))*fc(i)+f0(i)+fv(i)*dq(i);
    end
end
