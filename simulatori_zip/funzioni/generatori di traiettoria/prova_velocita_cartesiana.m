function dx=prova_velocita_cartesiana(run_time,vel_max,angle,u)
% generatore di velocità cartesiane sinusoidali.
    run_time=run_time/2; %  per rendere la frequenza pari a 1/runtime
    cycle=floor((-1)^((u-mod(u,run_time))/run_time)); % per rendere il moto cw o ccw
    v=cycle*vel_max/2*(1-cos(2*pi/run_time*u)); %ampiezza della sinusoide nell'istante di tempo
    v_x=cos(angle)*v; 
    v_z=-sin(angle)*v;
    %dx=[v_x 0 v_z 0 0 0];
    dx=[0 0 0 0 v_x v_z]; %velocità cartesiana finale 
end
