function Q=controllo_interno_RNE_OSI(robot,K,ddqmax, U)
   %% inizialization 
        n=7;
         if length(U)<(5*n) %% se non è ancora partita l'intera simulazione tiene ferme le acc
         U=zeros(5*n,1);
         end
        U = U(:)';   
 	    q = U(1:n); %real_position
        q_des=U(n+1:2*n);         %position desired
		dq = U(2*n+1:3*n); %dq_des
 		dq_des = U(3*n+1:4*n); %real_velocities
        ddq_des=U(4*n+1:5*n);         %acc desired

    
        de=dq_des-dq; % velocity error
        e=q_des-q; %position error
        D=2*sqrt(K); %Damping vector from stiffness
        Kp=diag(K); % Stiffness matrix
        Kd=diag(D); % Damping matrix
        Y=ddq_des'+Kd*de'+Kp*e'; % desired torque 
        y=zeros(7,1);
        for i=1:n
        y(i)=min(max(-ddqmax(i),Y(i)),ddqmax(i));
        end
        comp=rne(robot,q,dq,zeros(1,n))+function_friction_OSI(dq); %compensazione degle effetti fisici
         Q=robot.rne(ones(n,1)*q,zeros(n,n),eye(n),'gravity',[0,0,0])*y+comp'; %torque imposed on robot
end