function dq=prova_velocita_sinusoidale(run_time,omega_max,u)
    run_time=run_time/2;
    cycle=floor((-1)^((u-mod(u,run_time))/run_time));
    omega=cycle*omega_max/2*(1-cos(2*pi/run_time*u));
    dq_temp=[0 0 0 omega omega omega omega 0];
    dq=dq_temp;
end