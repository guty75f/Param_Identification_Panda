%plotta 4 figure e restituisce il MSE ( errore quadratico medio)
%le 4 figure sono:
%l'errore misurato sottratto dell'errore stimato
%l'errore misurato e stimato in funzione della velocità
%l'erorre rimanente
% le coppie simulate e reali
function [MSE]=plot_risultati(phi1,phi2,phi3,tau_stack,taux)

global dq TESTNUM
dqx=-2.1:0.01:2.1;

figure("name", "error and estimated tau")
for T=1:7
    TESTNUM=T;
    etau(:,TESTNUM)=error_tau;
    Tau_f(TESTNUM,:) = friction_model_single(phi1(TESTNUM), phi2(TESTNUM), phi3(TESTNUM)); 
    subplot(4,2,TESTNUM)
        hold on
        grid on
        plot(etau(:,TESTNUM))
        plot(Tau_f(TESTNUM,:))
        xlabel('samples [#]');
        ylabel('torque [Nm]');
        legend('error','estimated friction');
end

figure("Name","tau from dq")
for T=1:7
    TESTNUM=T;
    for i=1:length(dqx)
        tauf(i)=phi1(TESTNUM)*(1/(1+exp(-phi2(TESTNUM)*(phi3(TESTNUM)+dqx(i))))-1/(1+exp(-phi2(TESTNUM)*phi3(TESTNUM))));
    end
    subplot(4,2,TESTNUM)
        grid on
        hold on
        plot(dq(:,TESTNUM),etau(:,TESTNUM),":","Color",'r')
        plot(dqx',tauf',"Color",'b')
        xlabel('dq [rad/s]');
        ylabel('torque [Nm]');
        legend('real point','estimated friction');
end

figure ("Name","errore rimasto")
for T=1:7
    TESTNUM=T;
    Err(:,TESTNUM)=etau(:,TESTNUM)-Tau_f(TESTNUM,:)'; %errore rimanente dopo il modello del attriro
    subplot(4,2,TESTNUM)
        hold on
        grid on
        plot(Err(:,TESTNUM))
        xlabel('samples [#]');
        ylabel('torque [Nm]');
        legend('error');
end
 
figure ("Name","confronto tau misurata e modellizzata")
for T=1:7
    TESTNUM=T;
    subplot(4,2,TESTNUM)
        hold on
        grid on
        plot(tau_stack(1:10000,TESTNUM))
        plot(taux(1:10000,TESTNUM)+Tau_f(TESTNUM,1:10000)')
        xlabel('samples [#]');
        ylabel('torque [Nm]');
        legend('real tau','estimated tau');
end
for T=1:7
    MSE(T)=immse(etau(:,T),Tau_f(T,:)'); %calcolo errore quadratico medio
end
end
