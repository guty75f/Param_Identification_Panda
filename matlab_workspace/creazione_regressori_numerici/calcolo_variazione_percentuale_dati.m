%funzione che calcola la varianzia percentuale
function [perc,sigma2]=calcolo_variazione_percentuale_dati(PI,Tau_list, Ytot)
 SN=length(Tau_list); %n link* n samples
 nb=length(PI); % n base parameter
 sigma2=norm(Tau_list-(Ytot*PI))^2/(SN-nb);
 sigma=sqrt(diag(sigma2*inv(Ytot'*Ytot)));
 for i=1:length(sigma) 
     perc(i)=(sigma(i)/abs(PI(i)))*100; %calcolo std percentuale
 end 
end