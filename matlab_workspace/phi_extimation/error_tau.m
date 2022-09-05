%restituisce l'errore per il giunto preso in considearazione
function error_tau=error_tau()
    global EE TESTNUM
    error_tau=EE(:,TESTNUM);
end
        
