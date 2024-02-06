clc
clear
close all

x0=0;
x1=1;
x3=3;
x4=4;
x6=6;
x7=7;
x75=7.5;
x85=8.5;
xtot=10;
step=0.001;

%prima alzata da x0 a x1
[a1,v1,h1,t1,dt1,step1]= ProgettoPMF_function_ravaglia(5*(x1-x0),10*(x1-x0),15*(x1-x0),20*(x1-x0),25*(x1-x0),30*(x1-x0),36*(x1-x0),20,0,step);

%seconda alzata da x3 a x4
[a3,v3,h3,t3,dt2,step3]= ProgettoPMF_function_ravaglia(5*(x4-x3),10*(x4-x3),15*(x4-x3),20*(x4-x3),25*(x4-x3),30*(x4-x3),36*(x4-x3),30,0,step);

%prima discesa da x6 a x7
[a5,v5,h5,t5,dt5,step5]= ProgettoPMF_function_ravaglia(5*(x7-x6),10*(x7-x6),15*(x7-x6),20*(x7-x6),25*(x7-x6),30*(x7-x6),36*(x7-x6),-10,0,step);

%seconda discesa da x75 a x85
[a7,v7,h7,t7,dt7,step7]= ProgettoPMF_function_ravaglia(5*(x85-x75),10*(x85-x75),15*(x85-x75),20*(x85-x75),25*(x85-x75),30*(x85-x75),36*(x85-x75),-40,0,step);

%accelerazioni
Ae=[a1 linspace(0,0,(360*(x3-x1)/xtot)/step1) a3(2:end) linspace(0,0,(360*(x6-x4)/xtot)/step3) a5(2:end) linspace(0,0,(360*(x75-x7)/xtot)/step5) a7(2:end) linspace(0,0,(360*(xtot-x85)/xtot)/step7)];
A=Ae(1:end-1);
%velocità
Ve=[v1 linspace(0,0,(360*(x3-x1)/xtot)/step1) v3(2:end) linspace(0,0,(360*(x6-x4)/xtot)/step3) v5(2:end) linspace(0,0,(360*(x75-x7)/xtot)/step5) v7(2:end) linspace(0,0,(360*(xtot-x85)/xtot)/step7)];
V=Ve(1:end-1);
%alzata
He=[h1 linspace(h1(end),h1(end),(360*(x3-x1)/xtot)/step1) h3(2:end)+h1(end) linspace(h3(end)+h1(end),h3(end)+h1(end),(360*(x6-x4)/xtot)/step3) h3(end)+h1(end)+h5(2:end) linspace(40,40,(360*(x75-x7)/xtot)/step5) h3(end)+h1(end)+h5(end)+h7(2:end) linspace(0,0,(360*(xtot-x85)/xtot)/step7)];
H=He(1:end-1);
%scala dei tempi
T=linspace(0,360-step1,length(A));
T2=linspace(0,10-step1,length(A));

%Diagramma accelerazione
figure(1)
plot(T,A,'r','linewidth',1.5)
xlabel('Angolo [deg]')
ylabel('Accelerazione [mm/s^2]')
title('Diagramma Accelerazione','interpreter','latex','fontsize',16)
xlim([0 T(end)+step1])
grid on

%Diagramma accelerazione in gradi
figure(2);
plot(T2,A,'r','linewidth',1.5)
xlabel('Tempo [s]')
ylabel('Accelerazione [mm/s^2]')
title('Diagramma Accelerazione','interpreter','latex','fontsize',16)
xlim([0 T2(end)+step1])
grid on

%Diagramma Velocità
figure(3);
plot(T,V,'g','linewidth',1.5)
xlabel('Angolo [deg]')
ylabel('Velocità [mm/s]')
title('Diagramma Velocita','interpreter','latex','fontsize',16)
xlim([0 T(end)+step1])
grid on

figure(4);
plot(T2,V,'g','linewidth',1.5)
xlabel('Tempo [s]')
ylabel('Velocità [mm/s]')
title('Diagramma Velocita','interpreter','latex','fontsize',16)
xlim([0 T2(end)+step1])
grid on

%Diagramma Alzate
figure(5)
plot(T,H,'b','linewidth',1.5)
xlabel('Angolo [deg]')
ylabel('Alzata [mm]')
title('Diagramma Delle Alzate','interpreter','latex','fontsize',16)
xlim([0 T(end)+step1])
ylim([0 70])
grid on

figure(6);
plot(T2,H,'b','linewidth',1.5)
xlabel('Tempo [s]')
ylabel('Alzata [mm]')
title('Diagramma Delle Alzate','interpreter','latex','fontsize',16)
xlim([0 T2(end)+step1])
ylim([0 70])
grid on

%definisco il raggio di base della camma e il raggio della rotella -->
%valori a caso da cambiare
%[entrambi in mm]
Rb=77;
Rr=20;

r=Rb+H; %raggio della camma

xr=r.*cos(T*pi/180); %scrivo in coord. circolari l'andamento di x e y rotella
yr=r.*sin(T*pi/180);

%definisco l'angolo Beta da sottrarre il raggio della rotella alla camma
beta=linspace(1,1,length(A)); %faccio un vettore di partenza con tutti 1
beta(1)=pi/2+atan2((yr(1)-yr(length(A))),(xr(1)-xr(length(A)))); %determino il primo valore di beta per fare iterazione

for i=2:length(A)
    beta(i)=pi/2+atan2((yr(i)-yr(i-1)),(xr(i)-xr(i-1)));
end

%scrivo quindi i valori di X e Y corretti
X=xr+Rr*cos(beta);
Y=yr+Rr*sin(beta);
X1=medfilt1(X);
X1=X1';
Y1=medfilt1(Y);
Y1=Y1';
ap=zeros(1,length(A));

%determino l'angolo di pressione e verifico sia sotto i 35°
ap(1)=180/pi*atan2((H(1)-H(end)),sqrt(((xr(1)-xr(end))^2)+((yr(1)-yr(end))^2))); 
for i=1:length(H)-1
    ap(i+1)=180/pi*atan2((H(i+1)-H(i)),sqrt(((xr(i+1)-xr(i))^2)+((yr(i+1)-yr(i))^2)));
end
ap1=medfilt1(ap);

%massimo ap
[apmax,index]=max(ap1);
apmin=min(ap1);

if apmax > 35
    disp('L''angolo di pressione massimo apmax è maggiore di 35°,riguardare le dimensioni')
else
    fprintf('L''angolo apmax è %g°\n',apmax)
end


%grafico AP
figure(8);
plot(T,ap1,'c','linewidth',1.5)
hold on
plot(T(index),apmax,'*r')
ylim([-40 40])
xlim([0 360])
xlabel('Angolo [deg]')
ylabel('Angolo di pressione [deg]')
title('Diagramma Angolo Di Pressione','interpreter','latex','fontsize',16)
grid on
angolomax=linspace(35,35,length(A)); %linea limite angolo di pressione 35°
plot(T,angolomax,'m','linewidth',1.5)
cx=linspace(-Rr/2,Rr/2,1000);
cy=linspace(-Rr/2,Rr/2,1000);
const=linspace(0,0,1000);

%grafico Camma
figure(9);
plot((Rb)*sin(T),(Rb)*cos(T),'b');
hold on
axis square
plot(xr,yr,'k')
hold on
plot(X1,Y1,'g','linewidth',2)
plot(cx,const,'k','linewidth',1.5)
plot(const,cy,'k','linewidth',1.5)
xlim([-150 150])
ylim([-150 150])
title('Profilo Camma','interpreter','latex','fontsize',16)
legend('Raggio di Base','Profilo Centro Rotella','Profilo Camma','location','best')
grid on


%% Curvatura

% Inizializzazione vettori per calcolo curvatura
ds1=zeros(1,length(T));
ds2=zeros(1,length(T));
dang=zeros(1,length(T));
segno=zeros(1,length(T));
rho=zeros(1,length(T));
rho_abs=zeros(1,length(T));

for n=2:(length(T)-2)

    ds1(n)=sqrt((X1(n)-X1(n+1))^2+((Y1(n)-Y1(n+1))^2));
    ds2(n)=sqrt((X1(n+1)-X1(n+2))^2+((Y1(n+1)-Y1(n+2))^2));
    dang(n)=acos((((X1(n)-X1(n+1))*(X1(n+1)-X1(n+2)))+((Y1(n)-Y1(n+1))*(Y1(n+1)-Y1(n+2))))/(ds1(n)*ds2(n)));
    segno(n)=sign((X1(n+1)-X1(n))*(Y1(n+1)-Y1(n))-(X1(n)-r(n-1))*(Y1(n)-Y1(n-1)));
    rho(n)=(segno(n)*(ds1(n)/dang(n)));
    rho_abs(n)=abs((segno(n)*(ds1(n)/dang(n))));

end

% Fattore di sottocampionamento
sottocampionamento = 69;

% Creazione del nuovo vettore alfa con sottocampionamento
alpha_sottocampionato = T(1:sottocampionamento:end);
rho_sottocampionato = rho(1:sottocampionamento:end);
rho_abs_sottocampionato = rho_abs(1:sottocampionamento:end);

% Creazione del grafico
figure
subplot(2,1,1)
plot(alpha_sottocampionato, rho_sottocampionato, 'r', 'LineWidth', 1.25);
title('Raggio di Curvatura','interpreter','latex','fontsize',16);
axis([0 360 -300 200]);
xlabel('Alpha [°]');
ylabel('Rho [mm]');
grid on;
subplot(2,1,2)
plot(alpha_sottocampionato, rho_abs_sottocampionato, 'LineWidth', 1.25)
title('Modulo Raggio di Curvatura','interpreter','latex','fontsize',16);
axis([0 360 0 200]);
xlabel('Alpha [°]');
ylabel('abs(Rho) [mm]');
grid on;

%% Function PMF
% function [a,v,h,T,dt,step] = ProgettoPMF_function_ravaglia(t1,t2,t3,t4,t5,t6,t7,H,hs,step)
% dt=t7;
% perc=1;
% ts=0; %time shift
% ANG=(pi/2)*perc; %questa e' la percentuale di angolo coinvolta 
%                  %col seno:90 per la legge trapezoidale modificata o 
%                  %sinusoidale,1 per la trapezia
% t=ts:step:t7; %definisco il vettore dei tempi
% 
% %Prima parte:determinare il ratio
% a=0;
% v=0;
% T=ts;
% 
% index_t1=find(t<=t1,1,'last');
% index_t2=find(t<=t2,1,'last');
% index_t3=find(t<=t3,1,'last');
% index_t4=find(t<=t4,1,'last');
% index_t5=find(t<=t5,1,'last');
% index_t6=find(t<=t6,1,'last');
% index_t7=find(t<=t7,1,'last');
% 
% for i=2:index_t1
%     T(i)=T(end)+step;
%     a(i)=sin((ANG)*(T(i)/t1));
%     v(i)=v(end)+a(i)*step;
% end
% 
% for i=index_t1:index_t2
%     T(i)=T(i-1)+step;
%     a(i)=a(end);
%     v(i)=v(end)+a(i)*step;
% end
% 
% for i=index_t2:index_t3
%     T(i)=T(i-1)+step;
%     a(i)=sin((ANG)*(1-(T(i)-t2)/(t3-t2)));
%     v(i)=v(end)+a(i)*step;
% end
% 
% % % Parte grafica
% % figure(1)
% % plot(T,a,'r');
% % title('Primo Tratto fino a t3')
% % grid
% % hold on
% % plot(T,v,'g');
% % hold on
% 
% for i=index_t3:index_t4
%     T(i)=T(i-1)+step;
% end
% 
% v1=v;
% V1=v1(end);
% 
% %reinizializzo variabili
% a=0;
% v=0;
% 
% for i=index_t4:index_t5
%     T(i)=T(i-1)+step;
%     a(i)=sin((ANG)*((T(i)-t4)/(t5-t4)));
%     v(i)=v(end)+a(i)*step;
% end
% 
% for i=index_t5:index_t6
%     T(i)=T(i-1)+step;
%     a(i)=a(end);
%     v(i)=v(end)+a(i)*step;
% end
% 
% for i=index_t6:index_t7
%     T(i)=T(i-1)+step;
%     a(i)=sin((ANG)*(1-(T(i)-t6)/(t7-t6)));
%     v(i)=v(end)+a(i)*step;
% end
% 
% v2=v;
% V2=v2(end);
% if V2==0
%     disp('non si può eseguire calcolo ratio, lo prendo uguale a 1 per vedere come prosegue');
%     ratio=1;
% else
%     ratio=V1/V2;
% end
% 
% % % Parte grafica
% % figure(2)
% % plot(T,a,'r');
% % title('Secondo Tratto')
% % grid
% % hold on
% % plot(T,v,'g');
% % hold on
% 
% %seconda parte:scale
% a=0;
% v=0;
% h=hs;
% T=ts;
% 
% for i=2:index_t1
%     T(i)=T(end)+step;
%     a(i)=sin((ANG)*(T(i)/t1));
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t1:index_t2
%     T(i)=T(i-1)+step;
%     a(i)=a(end);
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t2:index_t3
%     T(i)=T(i-1)+step;
%     a(i)=sin((ANG)*(1-(T(i)-t2)/(t3-t2)));
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t3:index_t4
%     T(i)=T(i-1)+step;
%     a(i)=a(end);
%     v(i)=v(end);
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t4:index_t5
%     T(i)=T(i-1)+step;
%     a(i)=-ratio*sin((ANG)*((T(i)-t4)/(t5-t4)));
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t5:index_t6
%     T(i)=T(i-1)+step;
%     a(i)=a(end);
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t6:index_t7
%     T(i)=T(i-1)+step;
%     a(i)=-ratio*sin((ANG)*(1-(T(i)-t6)/(t7-t6)));
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% hf=h(end);
% scale=H/hf;
% % % Parte grafica
% % figure(3)
% % plot(T,a,'r');
% % title('Legge di moto a 7 tratti senza scale')
% % grid
% % hold on
% % plot(T,v,'g');
% % hold on
% % plot(T,h,'b');
% % hold on
% 
% 
% %parte finale:dati corretti
% a=0;
% v=0;
% h=hs;
% T=ts;
% 
% for i=2:index_t1
%     T(i)=T(end)+step;
%     a(i)=scale*sin((ANG)*(T(i)/t1));
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t1:index_t2
%     T(i)=T(i-1)+step;
%     a(i)=a(end);
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t2:index_t3
%     T(i)=T(i-1)+step;
%     a(i)=scale*sin((ANG)*(1-(T(i)-t2)/(t3-t2)));
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t3:index_t4
%     T(i)=T(i-1)+step;
%     a(i)=a(end);
%     v(i)=v(end);
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t4:index_t5
%     T(i)=T(i-1)+step;
%     a(i)=-ratio*scale*sin((ANG)*((T(i)-t4)/(t5-t4)));
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t5:index_t6
%     T(i)=T(i-1)+step;
%     a(i)=a(end);
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% for i=index_t6:index_t7
%     T(i)=T(i-1)+step;
%     a(i)=-ratio*scale*sin((ANG)*(1-(T(i)-t6)/(t7-t6)));
%     v(i)=v(end)+a(i)*step;
%     h(i)=h(end)+v(i)*step;
% end
% 
% end

Z1=linspace(0,0,length(X1));
Z1=Z1';