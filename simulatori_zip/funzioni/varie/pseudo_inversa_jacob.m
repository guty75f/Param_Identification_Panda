function Jinv=pseudo_inversa_jacob(robot, u)
%computation of pseudo inverse jacobian
   %% inizialization
    n=length(u);
    if n<8 || ( u(1)==0 && u(2)==0 && u(3)==0 && u(4)==0 && u(5)==0 && u(6)==0 ) % per inizializzare
        u= [1, -pi/4, 0, -3 * pi/4, 0, pi/2, pi/4, -pi/4]; %initial position
    end
   %% Weight Matrix
    W=eye(n); 
    W(5,5)=5;
   % W(4,4)=3;
    W(3,3)=1000;
    W(1,1)=0.9;
    W(8,8)=100000;
    %% Computation
    J=robot.jacob0(u); %robot geometry jacobian by RTB
    Jinv=inv(W)*J'*inv(J*inv(W)*J'); %pseudo inverse weighted jacobian
end