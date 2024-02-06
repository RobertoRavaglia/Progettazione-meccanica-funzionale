function [a,v,h,T,dt,step] = ProgettoPMF_function_ravaglia(t1,t2,t3,t4,t5,t6,t7,H,hs,step)
dt=t7;
perc=1;
ts=0; %time shift
ANG=(pi/2)*perc; %questa e' la percentuale di angolo coinvolta 
                 %col seno:90 per la legge trapezoidale modificata o 
                 %sinusoidale,1 per la trapezia
t=ts:step:t7; %definisco il vettore dei tempi

%Prima parte:determinare il ratio
a=0;
v=0;
T=ts;

index_t1=find(t<=t1,1,'last');
index_t2=find(t<=t2,1,'last');
index_t3=find(t<=t3,1,'last');
index_t4=find(t<=t4,1,'last');
index_t5=find(t<=t5,1,'last');
index_t6=find(t<=t6,1,'last');
index_t7=find(t<=t7,1,'last');

for i=2:index_t1
    T(i)=T(end)+step;
    a(i)=sin((ANG)*(T(i)/t1));
    v(i)=v(end)+a(i)*step;
end

for i=index_t1:index_t2
    T(i)=T(i-1)+step;
    a(i)=a(end);
    v(i)=v(end)+a(i)*step;
end
    
for i=index_t2:index_t3
    T(i)=T(i-1)+step;
    a(i)=sin((ANG)*(1-(T(i)-t2)/(t3-t2)));
    v(i)=v(end)+a(i)*step;
end

% % Parte grafica
% figure(1)
% plot(T,a,'r');
% title('Primo Tratto fino a t3')
% grid
% hold on
% plot(T,v,'g');
% hold on

for i=index_t3:index_t4
    T(i)=T(i-1)+step;
end

v1=v;
V1=v1(end);

%reinizializzo variabili
a=0;
v=0;

for i=index_t4:index_t5
    T(i)=T(i-1)+step;
    a(i)=sin((ANG)*((T(i)-t4)/(t5-t4)));
    v(i)=v(end)+a(i)*step;
end

for i=index_t5:index_t6
    T(i)=T(i-1)+step;
    a(i)=a(end);
    v(i)=v(end)+a(i)*step;
end

for i=index_t6:index_t7
    T(i)=T(i-1)+step;
    a(i)=sin((ANG)*(1-(T(i)-t6)/(t7-t6)));
    v(i)=v(end)+a(i)*step;
end

v2=v;
V2=v2(end);
if V2==0
    disp('non si può eseguire calcolo ratio, lo prendo uguale a 1 per vedere come prosegue');
    ratio=1;
else
    ratio=V1/V2;
end

% % Parte grafica
% figure(2)
% plot(T,a,'r');
% title('Secondo Tratto')
% grid
% hold on
% plot(T,v,'g');
% hold on

%seconda parte:scale
a=0;
v=0;
h=hs;
T=ts;

for i=2:index_t1
    T(i)=T(end)+step;
    a(i)=sin((ANG)*(T(i)/t1));
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end

for i=index_t1:index_t2
    T(i)=T(i-1)+step;
    a(i)=a(end);
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end
    
for i=index_t2:index_t3
    T(i)=T(i-1)+step;
    a(i)=sin((ANG)*(1-(T(i)-t2)/(t3-t2)));
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end

for i=index_t3:index_t4
    T(i)=T(i-1)+step;
    a(i)=a(end);
    v(i)=v(end);
    h(i)=h(end)+v(i)*step;
end

for i=index_t4:index_t5
    T(i)=T(i-1)+step;
    a(i)=-ratio*sin((ANG)*((T(i)-t4)/(t5-t4)));
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end

for i=index_t5:index_t6
    T(i)=T(i-1)+step;
    a(i)=a(end);
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end

for i=index_t6:index_t7
    T(i)=T(i-1)+step;
    a(i)=-ratio*sin((ANG)*(1-(T(i)-t6)/(t7-t6)));
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end
    
hf=h(end);
scale=H/hf;
% % Parte grafica
% figure(3)
% plot(T,a,'r');
% title('Legge di moto a 7 tratti senza scale')
% grid
% hold on
% plot(T,v,'g');
% hold on
% plot(T,h,'b');
% hold on


%parte finale:dati corretti
a=0;
v=0;
h=hs;
T=ts;

for i=2:index_t1
    T(i)=T(end)+step;
    a(i)=scale*sin((ANG)*(T(i)/t1));
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end

for i=index_t1:index_t2
    T(i)=T(i-1)+step;
    a(i)=a(end);
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end
    
for i=index_t2:index_t3
    T(i)=T(i-1)+step;
    a(i)=scale*sin((ANG)*(1-(T(i)-t2)/(t3-t2)));
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end

for i=index_t3:index_t4
    T(i)=T(i-1)+step;
    a(i)=a(end);
    v(i)=v(end);
    h(i)=h(end)+v(i)*step;
end

for i=index_t4:index_t5
    T(i)=T(i-1)+step;
    a(i)=-ratio*scale*sin((ANG)*((T(i)-t4)/(t5-t4)));
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end

for i=index_t5:index_t6
    T(i)=T(i-1)+step;
    a(i)=a(end);
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end

for i=index_t6:index_t7
    T(i)=T(i-1)+step;
    a(i)=-ratio*scale*sin((ANG)*(1-(T(i)-t6)/(t7-t6)));
    v(i)=v(end)+a(i)*step;
    h(i)=h(end)+v(i)*step;
end

end
