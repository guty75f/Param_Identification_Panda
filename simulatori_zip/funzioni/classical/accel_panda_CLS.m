function qdd = accel_panda_CLS(Q)
%compute the acceleration of panda robot from torque imposed
%% inizialization
    n=7; %DOF 
    if length(Q)<(3*n) %% se non è ancora partita l'intera simulazione tiene ferme le acc
        qdd=zero(n,1);
        return
    end
    Q = Q(:)';   % make it a row vector
 	q = Q(1:n); %position
	dq = Q(n+1:2*n); %vel
 	torque = Q(2*n+1:3*n); %torque imposed by control
%% computation
    M = Funzione_inerzia_classical(q);  	% compute current manipulator inertia
    % compute gravity, coriolis and friction torque 
    tau=Vettore_coriolis_classical(dq,q)+Gravita_classical(q,9.81)+function_friction_CLS(dq)';
    % torque-tau is the torque computed for the control without compensation
   qdd=(M\(torque-tau')'); %acc computed with rigid dynamic model
end
