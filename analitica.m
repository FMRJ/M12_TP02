function [T,Dt,Dx] = analitica(Ti,dt,dL)
% Calcula la solución analítica de la temperatura
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
Dt=0:dt:Ti;
Dx=0:dL:L;
Mt=max(size(Dt));
Mx=max(size(Dx));
T=zeros(Mx,Mt);
T(1,:)=310;
T(Mx,:)=310;
for x=2:Mx-1
    T(x,1)=500;
end
for t=2:Mt
    for x=2:Mx-1
        est=T1+((T2-T1)*Dx(x)/L);
        for n=1:100
            if mod(n,2)==1
                a_n(n)=(2/(n*pi))*(2*T0-T1-T2);
            elseif mod(n,2)==0
                a_n(n)=(2/(n*pi))*(T2-T1);
            end
            ex(n)=alfa*(n^2)*(pi^2)*Dt(t)/(L^2);
            S(n)=a_n(n)*exp(-1*ex(n))*sin(n*pi*Dx(x)/L);
        end
        Sn=sum(S);
        T(x,t)=est+Sn;
    end
end

end