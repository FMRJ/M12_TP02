%% Trabajo práctico 2 - Mecánica aplicada
clear all 
clc
% Datos
T1=310;         %[K]
T0=500;         %[K]
T2=310;         %[K]
rho=1100;       %[kg/m^3]
Cp=280;         %[J/kg*K]
k=101;          %[W/K*m]
L=2.1;          %[m]
alfa=k/(rho*Cp);
DeltaT=(L*2)/4*alfa;

%% Solución analítica
Ti=12000;
dt=2;
dL=L/40;
Dt=0:dt:Ti;
Dx=0:dL:L;
Mt=max(size(Dt));
Mx=max(size(Dx));
[T] = analitica(Ti,dt,dL);

%% Solución MDF
Q=zeros(Mx,Mt);
Th=zeros(Mx,Mt);
Th(1,:)=310;
Th(Mx,:)=310;
for x=2:Mx-1
    Th(x,1)=500;
end
for t=1:Mt-1
    for x=2:Mx-1
        Th(x,t+1)=Th(x,t)+(dt/(rho*Cp))*(k*((Th(x+1,t)-2*Th(x,t)+Th(x-1,t))/(dL^2))+Q(x,t));
    end
end

% Con termino fuente
Q(:,:)=-27400; %[J/s*m^3]
Thq=zeros(Mx,Mt);
Thq(1,:)=310;
Thq(Mx,:)=310;
for x=2:Mx-1
    Thq(x,1)=500;
end
for t=1:Mt-1
    for x=2:Mx-1
        Thq(x,t+1)=Thq(x,t)+(dt/(rho*Cp))*(k*((Thq(x+1,t)-2*Thq(x,t)+Thq(x-1,t))/(dL^2))+Q(x,t));
    end
end
%% Comparativa de resultados
[~,X, Y]=textread('Solidworks1.txt', '%f%f%f','delimiter',',');

% figure(1)
% subplot(2,1,1)
% plot(Dt,T(11,:),'r-')
% hold on
% grid on
% plot(X,Y,'-g')
% title('Solución analítica vs Simulación en  Solidworks');
% ylabel('Temperatura [K]');
% xlabel('Tiempo [s]');
% legend('Analítica', 'Solidworks');
% 
% subplot(2,1,2)
% plot(Dt,T(11,:),'r-')
% hold on
% grid on
% plot(Dt,Th(11,:),'-b')
% title('Solución analítica vs Aproximación MDF');
% ylabel('Temperatura [K]');
% xlabel('Tiempo [s]');
% legend('Analítica', 'Diferencias finitas');

%% Errores
% MDF
for j=1:Mt
    ErrM(j)=100*abs(T(11,j)-Th(11,j))/abs(T(11,j));
end

% SW
dt2=50;
dL2=0.06804354;
Dt2=0:dt:Ti;
Dx2=0:dL:L;
[T_2] = analitica(Ti,dt2,dL2);
for j=1:240
    ErrS(j)=100*abs(T_2(16,j)-Y(j,1))/abs(T_2(1,j));
end
fprintf(' Máximo error relativo porcentual en MDF: %4.3f %% \n',max(ErrM))
fprintf(' Máximo error relativo porcentual en SolidWorks: %4.3f %% \n',max(ErrS))

%% Variaciones temporales

[~,T50,X50,~,~]=textread('T50.txt', '%f%f%f%f%f','delimiter',',');
[~,T2000,X2000,~,~]=textread('T2000.txt', '%f%f%f%f%f','delimiter',',');
[~,T6000,X6000,~,~]=textread('T6000.txt', '%f%f%f%f%f','delimiter',',');

X50=X50+1050;
X50=flip(X50);
T50=flip(T50);
T2000=flip(T2000);
T6000=flip(T6000);

% % T=50
% figure(2)
% movegui southwest
% plot(X50*10^-3,T50,'g-')
% hold on 
% grid on
% plot(Dx,T(:,25),'r-')
% plot(Dx,Th(:,25),'b-')
% title('Soluciones para t=50');
% ylabel('Temperatura [K]');
% xlabel('Longitud [m]');
% legend('SolidWorks','Analítica','Diferencias finitas');
% 
% % T=2000
% figure(3)
% movegui center
% plot(X50*10^-3,T2000,'g-')
% hold on 
% grid on
% plot(Dx,T(:,1000),'r-')
% plot(Dx,Th(:,1000),'b-')
% title('Soluciones para t=2000');
% ylabel('Temperatura [K]');
% xlabel('Longitud [m]');
% legend('SolidWorks','Analítica','Diferencias finitas');
% 
% % T=6000
% figure(4)
% movegui southeast
% plot(X50*10^-3,T6000,'g-')
% hold on 
% grid on
% plot(Dx,T(:,3000),'r-')
% plot(Dx,Th(:,3000),'b-')
% title('Soluciones para t=6000');
% ylabel('Temperatura [K]');
% xlabel('Longitud [m]');
% legend('SolidWorks','Analítica','Diferencias finitas');

%% Caso con Término fuente

[~,Xq, Yq]=textread('Qt.txt', '%f%f%f','delimiter',',');
[~,Tq,XXq,~,~]=textread('QL.txt', '%f%f%f%f%f','delimiter',',');

figure(5)
subplot(2,1,1)
plot(Xq,Yq,'r-')
hold on
grid on
plot(Dt,Thq(21,:),'-b')
title('Simulación SolidWorks vs Aproximación MDF');
ylabel('Temperatura [K]');
xlabel('Tiempo [s]');
legend('SolidWorks', 'Diferencias finitas');

