% calcola le posizioni e gli assi z da DH
%è del tipo: [p,z]= calcolo_pos_z(T,N)
function [p,z]= calcolo_pos_z(T,N)
% N numero di giunti
% T trasformazione omogena rispetto a 0
% p(:,i) ritorna il punto di origine della terna i rispetto alla terna 0
% z(:,i) ritorna la rotazione su z della terna i rispetto alla terna 0

for i=1:N
    p(:,i)=T(1:3,4*i);
    z(:,i)=T(1:3,(4*i-1));
    avanzamento_calcolo_z=i;
end

end
