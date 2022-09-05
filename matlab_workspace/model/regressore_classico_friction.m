% calcola il regressore classico e i valori di PI. quelli s sono senza i valori nulli tau=Y*PI.
% del tipo:[Yc,PIc,Ycs,PIcs]=regressore_classico_friction(JJ,N,B,ddq,C2,G2,m,com,tipogiunto,theta,alphad,a,d, fc, fv ,f0,dq);
function [Yc,PIc,Ycs,PIcs]=regressore_classico_friction(JJ,N,B,ddq,C2,G2,m,com,tipogiunto,theta,alphad,a,d, fc, fv ,f0,dq)
tau_fc=tau_friction_classic(dq,f0,fc,fv,N);
tau_class=B*ddq'+C2+G2+tau_fc';
Jxx=JJ(:,1)';
Jyy=JJ(:,2)';
Jzz=JJ(:,3)';
Jxy=JJ(:,4)';
Jxz=JJ(:,5)';
Jyz=JJ(:,6)';
for i=1:N
    if (tipogiunto(i)==1) % se rotazionale creo i coseni e i seni come incognite
   
        c(i)=cos(theta(i));
    
        s(i)=sin(theta(i));
    end
    if (tipogiunto(i)==0) % se prismatico i coseni e i seni son calcolati in deg per evitare che sin(0)!=0 e sin(pi/2)!=1
        thetad=rad2deg(theta(i));
        c(i)=cosd(thetad);
        s(i)=sind(thetad);
    end

end
    n_c=11;
    tot_int=N*n_c;
    for i=1:N
                cunt=0;
        tm(i,:)=diff(tau_class,m(i)); %derivata su mi
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        tmcx(i,:)=diff(tm(i,:),com(1,i)); %derivata su mi e cxi
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        tmcy(i,:)=diff(tm(i,:),com(2,i)); %derivata su mi e cyi
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;

        tIxx(i,:)=diff(tau_class,Jxx(i)); %moltiplicatore per Jxxi
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        tIzz(i,:)=diff(tau_class,Jzz(i)); %moltiplicatore di Jzzi
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        tIxy(i,:)=diff(tau_class,Jxy(i)); %moltiplicatore di Jxyi
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        tIxz(i,:)=diff(tau_class,Jxz(i)); %moltiplicatore di Jxzi
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
        tIyz(i,:)=diff(tau_class,Jyz(i)); %moltiplicatore di Iyzi
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
         tfc(i,:)=diff(tau_class,fc(i)); %moltiplicatore di fc
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
         tfv(i,:)=diff(tau_class,fv(i)); %moltiplicatore di fc
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
         tf0(i,:)=diff(tau_class,f0(i)); %moltiplicatore di fc
                avanzamento_diff_regr_class=[100*(n_c*(i-1)+cunt)/tot_int]
                cunt=cunt+1;
    end
mz=sym(zeros(1,N));
my=sym(zeros(1,N));
mx=sym(zeros(1,N));
mm=m;
f=alphad;
for i=1:N
 mz(i)=mm(i)*com(3,i);
    my(i)=mm(i)*com(2,i);
    mx(i)=mm(i)*com(1,i);
end
for i=N:-1:2
   
    Jxx(i)=Jxx(i)-Jyy(i);
    Jxx(i-1)=Jxx(i-1)+Jyy(i)+2*d(i)*mz(i)+d(i)^2*mm(i);
    Jxy(i-1)=Jxy(i-1)+a(i)*sind(f(i))*mz(i)+a(i)*d(i)*sind(f(i))*mm(i);
    Jxz(i-1)=Jxz(i-1)-a(i)*cosd(f(i))*mz(i)-a(i)*d(i)*cosd(f(i))*mm(i);
    Jyy(i-1)=Jyy(i-1)+cosd(f(i))^2*Jyy(i)+2*d(i)*cosd(f(i))^2*mz(i)+mm(i)*(a(i)^2+d(i)^2*cosd(f(i))^2);
    Jyz(i-1)=Jyz(i-1)+cosd(f(i))*sind(f(i))*Jyy(i)+2*d(i)*cosd(f(i))*sind(f(i))*mz(i)+d(i)^2*cosd(f(i))*sind(f(i))*mm(i);
    Jzz(i-1)=Jzz(i-1)+sind(f(i))^2*Jyy(i)+2*d(i)*sind(f(i))^2*mz(i)+mm(i)*(a(i)^2+d(i)^2*sind(f(i))^2);
    mx(i-1)=mx(i-1)+a(i)*mm(i);
    my(i-1)=my(i-1)-sind(f(i))*mz(i)-d(i)*sind(f(i))*mm(i);
    mz(i-1)=mz(i-1)+cosd(f(i))*mz(i)+d(i)*cosd(f(i))*mm(i);
    mm(i-1)=mm(i-1)+mm(i);
end
V=N+3;
for i=1:N
        PIc(1,(V*i-(V-1)):V*i)=[Jxx(i), Jxy(i), Jxz(i), Jyz(i), Jzz(i), mx(i), my(i),fc(i), fv(i), f0(i)];
         Yc(:,(V*i-(V-1)):V*i)=[tIxx(i,:)',tIxy(i,:)', tIxz(i,:)', tIyz(i,:)', tIzz(i,:)',tmcx(i,:)',tmcy(i,:)', tfc(i,:)', tfv(i,:)', tf0(i,:)']; 
end
[Ycs,PIcs]=delete_zero(Yc,PIc,N);
end
