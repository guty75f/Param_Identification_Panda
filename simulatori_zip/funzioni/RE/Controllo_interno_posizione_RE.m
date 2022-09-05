function Q=controllo_interno_RE(K,ddqmax, U)
   %% inizialization 
        n=7;
        U = U(:)';   
 	    q = U(1:n) %real_position
        q_des=U(n+1:2*n)         %position desired
		dq = U(2*n+1:3*n); %dq_des
 		dq_des = U(3*n+1:4*n); %real_velocities
        ddq_des=U(4*n+1:5*n)         %acc desired

        de=dq_des-dq; % velocity error
        e=q_des-q %position error
        D=2*sqrt(K); %Damping vector from stiffness
        Kp=diag(K); % Stiffness matrix
        Kd=diag(D); % Damping matrix
        y=ddq_des'+Kd*de'+Kp*e'; % desired torque 
        for i=1:n;
        y(i)=min(max(-ddqmax(i),y(i)),ddqmax(i));
        end
        comp=Vettore_coriolis(dq,q)+gravita(q)'+function_friction_RE(dq); %compensazione degle effetti fisici
        Q=Funzione_inertia(q)*y+comp; %torque imposed on robot
end