%calcolo dell'attrito secondo i parametri dati e il modello:
% phi1*(1/(1+exp(-phi2*(phi3+dqx(i))))-1/(1+exp(-phi2*phi3)))
function tauf=friction_model_single(phi1, phi2, phi3)
    
global TESTNUM dq
    dqx=dq(:,TESTNUM);
for i=1:length(dqx)
    tauf(i)=phi1*(1/(1+exp(-phi2*(phi3+dqx(i))))-1/(1+exp(-phi2*phi3)));
end
