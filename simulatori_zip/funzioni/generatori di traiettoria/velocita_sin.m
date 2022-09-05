function dq=velocita_sin(A,T, dqmax,t)
n=7;
    for i=1:n
        dq(i)=max(-dqmax(i),(min(A(i)*sin((2*pi*t)/T(i)),dqmax(i))));
    end
    
end