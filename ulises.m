%% PR08 
%
% Agregar la expresión análitica de la serie en ambos formatos (hay dos series 
% en esta tarea).
%
% _ *1.-* _ Encuentra la expresion de la serie de Fourier, para la señal f(t) = t, 
% en el intervalo [−1, 1].Observa que es una señal impar.
%%

figure(1)
hFig = figure(1);
set(hFig, 'Position', [0 0 900 900])
f=@(u) u;
u=-1:0.01:1;
plot(u,f(u))
%% 
% Con:
% 
% $$T_0 =2,\;\;\;\;\;\;\;W_0 =\pi$$
% 
% Siendo la serier de fourier compacta:
% 
% $$S_f \left(t\right)=C_0 \;\mathrm{cos}\left(\phi_0 \right)+\sum_{n=1}^{\infty 
% } \left(C_n \;\mathrm{cos}\left(n\varpi_0 t\;+\phi_n \right)\right)$$
% 
% Asi obteniendo los valores necesarios para la serie de fourier trigonometrica 
% compacta:
%%

syms t n 
assume n positive integer

T0=2;
w0=2*pi/T0;

a0=(1/T0)*int(t,t,-1,1);
an=simplify((2/T0)*int(t*cos(n*w0*t),t,-1,1));
bn=simplify((4/T0)*int(t*sin(n*w0*t),t,0,1));

c0=abs(a0)
Cn=sqrt(simplify(an^2+bn^2))
cn=@(n) 2./(pi.*n);
if(a0<0)
    
fi0=pi; 

else
    fi0=0;
end

%%
%como fin es algo entre 0, se reescribira para adecuarla

fin=@(n) pi/2.*((-1).^n);

%% 
% a.- 
% 
% Realizamos la grafica
% 
% $$\begin{array}{l}S_f \left(t\right)=C_0 \;\mathrm{cos}\left(\phi_0 \right)+\sum_{n=1}^{\infty 
% } \left(C_n \;\mathrm{cos}\left(n\varpi_0 t\;+\phi_n \right)\right)\\=\sum_{n=1}^{\infty 
% } \left(\frac{2}{n\pi }\;\mathrm{cos}\left(n\pi t+\phi_n \right)\right)\end{array}$$
%%
t0=-1;
tf=1;
armo=15;
a=-7;
b=7;

sf=c0*cos(fi0);
t=a:0.01:b;

for m=1:armo
    theta=w0*m*t + fin(m);
    sf=sf + cn(m)*cos(theta);
end
figure(2)
hFig = figure(2);
set(hFig, 'Position', [0 0 900 900])
subplot(3,2,1)
plot(t,sf,'LineWidth',2)
grid on
legend('Serie de Fourier','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)

sf=c0*cos(fi0);
t1=t0:0.01:tf;

for m=1:armo
    theta=w0*m*t1 + fin(m);
    sf=sf + cn(m)*cos(theta);
end

subplot(3,2,2)
plot(t1,f(t1),'b','LineWidth',2)
grid on
hold on
plot(t1,sf,'LineWidth',2)
legend('Función original','Serie de Fourier ','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
nn=-0:armo;
axis auto

subplot(3,2,4)
e=f(t1)-sf;
plot(t1,e,'LineWidth',2)
title('Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

subplot(3,2,6)
e=f(t1)-sf;
area(t1,e.^2)
legend('Energia del error','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

cnM=zeros(1,length(nn));
cont=1;
for i =-armo:armo
    if i==0
        cnM(cont)=c0;
    end
    cnM(cont)=cn(i);
    cont=cont+1;
end

subplot(3,2,3)
stem(w0*nn,cn(nn),'LineWidth',2)
title('Espectro de magnitud Cn','FontWeight','bold','FontSize',16)
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on

subplot(3,2,5)
stem(w0*nn,fin(nn),'LineWidth',2)
title('Espectro de fase, \phi_n ','FontWeight','bold','FontSize',16) % % 
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on

%% 
% b. Serie de Fourier Exponencial Compleja:
%%


syms t n 
assume n positive integer
T0=2;
W0=2*pi/T0;
f=t;
Dn=simplify((1/T0)*int(f*exp(-n*j*W0*t),-T0/2,T0/2))

%% 
% Las Serie de Fourier queda de la siguiente forma:
% 
% 
% 
% Realizamos la grafica
%%

m=15;
Dn=@(n) j*(-1)^(n)/(n*pi);
D0=0;
t=(5)*(-T0/2):0.01:(5)*(T0/2);
Sfc=D0;
for n=1:m
    Sfc=Sfc+Dn(-n)*exp(W0*-n*t*j)+Dn(n)*exp(W0*n*t*j);
end
figure(3)
hFig = figure(3);
set(hFig, 'Position', [0 0 900 900])

subplot(3,2,1)
plot(t,Sfc,'LineWidth',2)
grid on
legend('S. Fourier Exp. Compleja','Location','southeast')
xlabel('t','FontWeight','bold','FontSize',16)
axis([-6 6 -2 2])
Sfc=D0;
t=(-T0/2):0.01:(T0/2);
for n=1:m
    Sfc=Sfc+Dn(-n)*exp(W0*-n*t*j)+Dn(n)*exp(W0*n*t*j);
end
subplot(3,2,2)
f=@(t) t;
plot(t,f(t),'r','LineWidth',1)
grid on
hold on
plot(t,Sfc,'g','LineWidth',0.85)
legend('Función original','Serie de Fourier ','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
axis([-1.2 1.2 -1.2 1.2])
subplot(3,2,4)
Ec=f(t)-Sfc;
plot(t,Ec,'LineWidth',2)
title('Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on
subplot(3,2,6)
area(t,Ec.^2)
legend('Energia del error','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on
nn=-m:m;
absdn=zeros(1,length(nn));
cont=1;
for i =-m:m
    if i==0
        absdn(cont)=D0;
    end
    
    absdn(cont)=Dn(i);
    cont=cont+1;
end
subplot(3,2,3)
stem(W0*nn,abs(absdn),'LineWidth',2)
title('Espectro de magnitud D_n ','FontWeight','bold','FontSize',16)
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on
axis auto
subplot(3,2,5)
stem(W0*nn,angle(absdn),'LineWidth',2)
title('Espectro de fase, \angle de D_n ','FontWeight','bold','FontSize',16)
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on
axis auto

%% 
% _ * 2.- * _ Encuentra la expresion de la serie de Fourier, para la señal f(t) = t^2  
% en el intervalo [−2, 2]. Observa que es una señ̃al par.
%%


figure(4)
hFig = figure(4);
set(hFig, 'Position', [0 0 900 900])
f=@(u) u.^2;
u=-2:0.01:2;
plot(u,f(u))

%% 
% Realizando lo mismo que el anterior
% 
% a.-Serie de Fourier Trigonometrica Compacta
%%
syms n t 

assume(n,["positive" "integer"])

T0=4;
w0=2*pi/T0;

a0=(1/T0)*int(t^2,t,-2,2);
an=simplify((2/T0)*int(t^2*cos(n*w0*t),t,-2,2));
bn=simplify((4/T0)*int(t^2*sin(n*w0*t),t,-2,2));

c0=abs(a0)
Cn=sqrt(simplify(an^2+bn^2))
cn=@(n) 16./((n.*pi).^2);
fin=-atan(bn/an);

if(a0<0)
    
fi0=pi
else
    fi0=0;
end
%como fin es o entre algo, se reescribira para adecuarla
fin=@(n) pi.*((mod(n,1)==0)) -pi.*((mod(n,2)==0));

%% 
% Las Serie de Fourier queda de la siguiente forma:
% 
% $$\begin{array}{l}S_f \left(t\right)=C_0 \;\mathrm{cos}\left(\phi_0 \right)+\sum_{n=1}^{\infty 
% } \left(C_n \;\mathrm{cos}\left(n\varpi_0 t\;+\phi_n \right)\right)\\=\frac{4}{3}+\sum_{n=1}^{\infty 
% } \left(\frac{16{\left(-1\right)}^n }{n^2 \pi^2 }\;\mathrm{cos}\left(\frac{n\pi 
% t\;}{2}+\phi_n \right)\right)\end{array}$$
% 
% 
% 
% Realizamos la grafica
%%
t0=-2;
tf=2;
armo=15;
a=-10;
b=10;

sf=c0*cos(fi0);
t=a:0.01:b;

for m=1:armo
    theta=w0*m*t + fin(m);
    sf=sf + cn(m)*cos(theta);
end
figure(5)
hFig = figure(5);
set(hFig, 'Position', [0 0 900 900])
subplot(3,2,1)
plot(t,sf,'LineWidth',2)
grid on
legend('Serie de Fourier','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)

sf=c0*cos(fi0);
t1=t0:0.01:tf;
i=1;
for m=1:armo
    theta=w0*m*t1 + fin(m);
    sf=sf + cn(m)*cos(theta);
end

subplot(3,2,2)
plot(t1,f(t1),'b','LineWidth',2)
grid on
hold on
plot(t1,sf,'LineWidth',2)
legend('Función original','Serie de Fourier ','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
nn=0:armo;
axis auto

subplot(3,2,4)
e=f(t1)-sf;
plot(t1,e,'LineWidth',2)
title('Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

subplot(3,2,6)
e=f(t1)-sf;
area(t1,e.^2)
legend('Energia del error','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

cnM=zeros(1,length(nn));
cont=1;
for i =-armo:armo
    if i==0
        cnM(cont)=c0;
    end
    cnM(cont)=cn(i);
    cont=cont+1;
end

subplot(3,2,3)
stem(w0*nn,cn(nn),'LineWidth',2)
title('Espectro de magnitud Cn','FontWeight','bold','FontSize',16)
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on

subplot(3,2,5)
stem(w0*nn,fin(nn),'LineWidth',2)
title('Espectro de fase, \phi_n ','FontWeight','bold','FontSize',16)
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on


%% 
% b. Serie de Fourier Exponencial Compleja:
%%

syms t n 
assume(n,'integer')
T0=4;
W0=2*pi/T0;
f2=t*t;
Dn=simplify((1/T0)*int(f2*exp(-n*j*W0*t),-T0/2,T0/2))
D0=simplify((1/T0)*int(f2,-T0/2,T0/2))
%% 
% Las Serie de Fourier queda de la siguiente forma:
% 
% $$S_f \left(t\right)=\sum_{n=-\infty }^{\infty } D_n e^{\mathrm{nj}\omega_0 
% t} =\frac{4}{3}+\sum_{n=-\infty }^{-1} \frac{8{\left(-1\right)}^n }{n^2 \pi^2 
% }e^{\frac{\pi }{2}\mathrm{njt}} +\sum_1^{\infty } \frac{8{\left(-1\right)}^n 
% }{n^2 \pi^2 }e^{\frac{\pi }{2}\mathrm{njt}}$$
% 
% Realizamos la grafica
%%
m=15;
Dn=@(n) (8*(-1)^n)/((n*pi)^2);
D0=4/3;
t=(5)*(-T0/2):0.01:(5)*(T0/2);
Sfc=D0;
for n=1:m
    Sfc=Sfc+Dn(-n)*exp(W0*-n*t*j)+Dn(n)*exp(W0*n*t*j);
end
figure(6)
hFig = figure(6);
set(hFig, 'Position', [0 0 900 900])
subplot(3,2,1)
plot(t,Sfc,'LineWidth',2)
grid on
legend('S. Fourier','Location','southeast')
xlabel('t','FontWeight','bold','FontSize',16)
axis([-10 10 -0.5 4.5])

Sfc=D0;
t=(-T0/2):0.01:(T0/2);
for n=1:m
    Sfc=Sfc+Dn(-n)*exp(W0*-n*t*j)+Dn(n)*exp(W0*n*t*j);
end

subplot(3,2,2)
f2=@(t) t.*t;
plot(t,f2(t),'r','LineWidth',1)
grid on
hold on
plot(t,Sfc,'b','LineWidth',0.85)
legend('Función original','Serie de Fourier ','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
axis([-2 2 -0.5 4.5])
subplot(3,2,4)
Ec=f2(t)-Sfc;
plot(t,Ec,'LineWidth',2)
title('Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

subplot(3,2,6)
area(t,Ec.^2)
legend('Energia del error','Location','Best')
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on
nn=-m:m;
absdn=zeros(1,length(nn));
cont=1;
for i =-m:m
    if i==0
        absdn(cont)=D0;
    end
    
    absdn(cont)=Dn(i);
    cont=cont+1;
end

subplot(3,2,3)
stem(W0*nn,abs(absdn),'LineWidth',2)
title('Espectro de magnitud D_n ','FontWeight','bold','FontSize',16)
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on
axis auto
subplot(3,2,5)
stem(W0*nn,angle(absdn),'LineWidth',2)
title('Espectro de fase, \angle de D_n ','FontWeight','bold','FontSize',16)
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on
axis auto
