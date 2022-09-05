%funzione per calcolare il loss con un dato set di parametri x;
%in poche parole la funzione obbiettivo da minimizzare in fminsearch
function [loss] = try_this_things(x)
phi1 = x(1);
phi2 = x(2);
phi3 = x(3);
Tau_f = friction_model_single(phi1, phi2, phi3);

etau=error_tau();
error = etau' - Tau_f;
loss = (error*error');
end