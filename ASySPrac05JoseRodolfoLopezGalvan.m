%% Práctica 5: Series de Fourier en tiempo continuo
% Integrantes:
%% 
% * Arellano Paz Angel Ulises
% * Cardoso Arias Javier
% * Gachuz Hernández Karla Denisse
% * García Arteaga Alejandro
% * López Galván José Rodolfo
%

%% Objetivos
%
% * Realizar gráficas de series de Fourier exponenciales y trigonométricas en 
% tiempo continuo
% * Manipulación de instrucciones en MATLAB
%

%% Ejercicio PR10
%
% Sea una señal $f\left(t\right)\;$de periodo $T$, su descripción en el intervalo 
% $\left(-\frac{T}{2},\frac{T}{2}\right)$ es:
% 
% $$f\left(t\right)={\textrm{ae}}^{-a\left|t\right|}$$
% 
% Su Serie de Fourier en dicho intervalo es:
% 
% $$S_f \left(t\right)=\sum_{n=-\infty }^{\infty } \frac{a^2 \left(1-e^{-a} 
% \cos \left(\pi n\right)\right)}{a^2 +n^2 \pi^2 }e^{\textrm{jn}\pi t}$$ 
% 
% *a)* Determine el valor de $T$
% 
% *b*) ¿Cúal es el valor promedio de $f\left(t\right)$? (realice dos procedimientos). 
% El valor promedio de una sẽnal periódica se define como:
% 
% $$V_p \left(f\right)=\frac{1}{T}\int_{<\textrm{To}>} f\left(t\right)\textrm{dt}$$
% 
% *c)* La componente de$f\left(t\right)$ en cierta frecuencia se puede ecpresar 
% como $\textrm{Acos}\left(3\pi t\right)$. Determine el valor de $A$.
% 
% *d)* Calcule la Serie de Fourier para la señal $f\left(t\right)$ con el dato 
% encontrado en *a)* y verifique que coincida con la proporcionada 
% 
% 
% 
% Solución:
% 
% Ya que $T$ es el perido de $f\left(t\right)$haremos la igualdad $T=\textrm{T0}$ 
% esto con la finalidad de poder respetar la nomenclatura de la Serie de Fourier
% 
% *a)* 
% 
% De la $S_f \left(t\right)\;$que nos fue proporcionada podemos notar que $\omega_0 
% =\pi$ ya que:
% 
% $$S_f \left(t\right)=D_0 +\sum_{n=-\infty }^{\infty } D_n e^{{\textrm{jn}\omega 
% }_0 t} =\sum_{n=-\infty }^{\infty } \frac{a^2 \left(1-e^{-a} \cos \left(\pi 
% n\right)\right)}{a^2 +n^2 \pi^2 }e^{\textrm{jn}\pi t} \;\;$$
% 
% Analizando las exponenciales
% 
% $${e^{\textrm{jn}\omega_0 t} \;=e}^{\textrm{jn}\pi t} \;\Longrightarrow \omega_0 
% =\pi \;$$

W0=pi;
T0=2*pi/W0;
T0
%% 
% *b)* Declaramos las variables simbolicas a utilizar y escribimos la integral

syms a t T0
f=a*exp(-a*abs(t))
Vp=simplify(((1/T0)*int(f,-T0/2,T0/2)))
%% 
% Intentamos simplificar pero el software no da para más
% 
% Si ocupamos $T=2\;$que fue obtenido en el inciso anterior, el valor promedio 
% de $f\left(t\right)$ es:

Vp=subs(Vp,T0,2)
%% 
% *c)* Notemos que $\textrm{Acos}\left(3\pi t\right)$ tiene la forma $\textrm{Acos}\left(n\pi 
% t\right)$ con $n=3$, dicha forma cumple con la forma Trigonometrica Compacta, 
% por ende $\textrm{Acos}\left(3\pi t\right)$ 
% 
% seria el tercer armonico y $A$ seria el tercer coeficiente, al cual llamaremos 
% $C_3$

syms a n 
syms a n
assume a n real
assumptions
Dn=(((a^2)*(1-exp(-a)*cos(pi*n)))/((a^2)+(pi*n)^2));
C3=2*abs(subs(Dn,n,3))
%% 
% Si buscamos obtener $\phi_n \;$ usamos $\phi_n =\angle \left(D_n \right)$

phi3=angle(subs(Dn,n,3))
%% 
% *d)* Usaremos $\textrm{T0}=2$.

syms n t a 
assume n integer
assume a real
T0=2;
W0=pi;
Dn=simplify((1/T0)*int(f*exp(-n*W0*j*t),-T0/2,T0/2))
%% 
% La expresión obtenida para los  $D_n$ es una expresión equivalente a la proporcionada
% 
% Para poder graficar haremos $a=2$

Dn=subs(Dn,a,2)
%% 
% Con lo que la Serie de Fourier exponenial compleja queda de la siguiente forma
% 
% $$S_f \left(t\right)=\sum_{n=-\infty }^{\infty } D_n e^{\textrm{jn}\pi t} 
% =\sum_{n=-\infty }^{\infty } \frac{4e^{-2} \left(e^2 -{\left(-1\right)}^n \right)}{4+\pi^2 
% n^2 }e^{\textrm{jn}\pi t}$$
% 
% Para obtener la Trigonometrica Compacta.

Cn=2*abs(simplify((Dn)));
Cn=simplify(Cn)
Phi_n=angle(Dn)
C0=abs(subs(Dn,n,0))
Phi_0=angle(subs(Dn,n,0))
%% 
% Utilizando los $C_n ,{\;\phi }_{n\;} ,\;\phi_0$ y $C_0$ calculados ya podemos 
% expresar la serie de Fourier en su forma Trigonometrica Compacta
% 
% $$S_f \left(t\right)=C_0 \cos \left(\phi_0 \right)+\sum_{n=1}^{\infty } C_n 
% \cos \left(n\omega_0 t+\phi_n \right)=e^{-2} \left(e^2 -1\right)+\sum_{n=1}^{\infty 
% } \frac{{8e}^{-2} \left(e^{-2} -{\left(-1\right)}^n \right)}{4+\pi^2 n^2 }\cos 
% \left(n\pi t\right)$$
% 
% Existe algo llamado Error que se define de la siguiente manera
% 
% $$E=f\left(t\right)-S_f \left(t\right)$$
% 
% El cual tambien sera graficado
% 
% _Exponencial Compleja_

T0=2;
W0=2*pi/T0;
a=2;
m=15;
f=@(t) a*exp(-a*abs(t));
Dn=@(n) ((4*exp(-2)*(exp(2)-(-1)^n))/(4+(n*pi)^2));
D0=exp(-2)*(exp(2)-1);
t= (2.5)*(-T0):0.001:(2.5)*(T0);
Sfc=D0;
for n=1:m
    Sfc=Sfc+Dn(-n)*exp(W0*-n*t*j)+Dn(n)*exp(W0*n*t*j);
end

figure (1);
hFig = figure(1);
set(hFig, 'Position', [0 0 900 900])
subplot(3,2,1)
plot(t,Sfc,'LineWidth',2)
grid on
legend('S. Fourier Exp. Compleja','Location','southeast')
xlabel('t','FontWeight','bold','FontSize',16)
axis([-5 5 -0.5 2.5])

Sfc=D0;
t=(-T0/2):0.0001:(T0/2);
for n=1:m
    Sfc=Sfc+Dn(-n)*exp(W0*-n*t*j)+Dn(n)*exp(W0*n*t*j);
end

subplot(3,2,2)
plot(t,f(t),'r','LineWidth',0.75)
grid on
hold on
plot(t,Sfc,'b','LineWidth',1.5)
legend('Función original','Serie de Fourier ','Location','northeast')
xlabel('t','FontWeight','bold','FontSize',16)
axis auto

subplot(3,2,4)
Ec=f(t)-Sfc;
plot(t,Ec,'LineWidth',2)
title('Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

subplot(3,2,6)
area(t,Ec.^2)
legend('Energia del error','Location','northeast')
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

subplot(3,2,5) % % 
stem(W0*nn,angle(absdn),'LineWidth',2) % % 
title('Espectro de fase, \angle de D_n ','FontWeight','bold','FontSize',16) % % 
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on
%% 
% _Trigonometrica compacta_

figure();
T0=2;
W0=2*pi/T0;
a=2;
m=15;
f=@(t) a*exp(-a*abs(t));
Cn=@(n) ((8*exp(-2)*(exp(2)-(-1)^n))/(4+(n*pi)^2));
C0=exp(-2)*(exp(2)-1);
t= (2.5)*(-T0):0.001:(2.5)*(T0);
phin=0;
phi0=0;
Sft=C0*cos(phi0);
for n=1:m
    Sft=Sft+Cn(n)*cos(n*pi*t+phin);
end
figure (2);
hFigg = figure(2);
set(hFigg, 'Position', [0 0 900 900])
subplot(3,2,1)
plot(t,Sft,'LineWidth',2)
grid on
legend('S. Fourier Trigo. Compacta','Location','southeast')
xlabel('t','FontWeight','bold','FontSize',16)
axis([-5 5 -0.5 2.5])

Sft=D0;
t=(-T0/2):0.0001:(T0/2);
for n=1:m
    Sft=Sft+Cn(n)*cos(n*pi*t+phin);
end

subplot(3,2,2)
plot(t,f(t),'r','LineWidth',0.75)
grid on
hold on
plot(t,Sft,'b','LineWidth',1.5)
legend('Función original','Serie de Fourier ','Location','northeast')
xlabel('t','FontWeight','bold','FontSize',16)
axis auto

subplot(3,2,4)
Et=f(t)-Sft;
plot(t,Et,'LineWidth',2)
title('Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

subplot(3,2,6)
area(t,Et.^2)
legend('Energia del error','Location','northeast')
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

nn=0:m;
abscn=zeros(1,length(nn));
cont=1;
for i =0:m
    if i==0
        abscn(cont)=C0;
    end
    
    abscn(cont)=Cn(i);
    cont=cont+1;
end
nn=0:m;
absdnn=zeros(1,length(nn));
cont=1;
for i =0:m
    if i==0
        absdnn(cont)=D0;
    end
    
    absdnn(cont)=Dn(i);
    cont=cont+1;
end


subplot(3,2,3)
stem(W0*nn,abs(abscn),'LineWidth',2)
title('Espectro de magnitud C_n ','FontWeight','bold','FontSize',16)
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on

subplot(3,2,5) % % 
stem(W0*nn,angle(absdnn),'LineWidth',2) % % 
title('Espectro de fase \phi_n ','FontWeight','bold','FontSize',16) % % 
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on


%% Ejercicio 3
%Ejercicio de Karla

clear;
figure();

T0=2*pi;
W0=2*pi/T0;
m=15;
f=@(t) (heaviside(t+pi/2)-heaviside(t-pi/2));
Dn=@(n) (sin((n*pi)/2)/(pi*n));
D0=1/2;
t= (2.5)*(-T0):0.001:(2.5)*(T0);
Sfc=D0;
for n=1:m
    Sfc=Sfc+Dn(-n)*exp(W0*-n*t*j)+Dn(n)*exp(W0*n*t*j);
end

figure (3);
hFiggg = figure(3);
set(hFiggg, 'Position', [0 0 900 900])
subplot(3,2,1)
plot(t,Sfc,'LineWidth',2)
grid on
legend('S. Fourier Exp. Compleja','Location','southeast')
xlabel('t','FontWeight','bold','FontSize',16)
axis([-2.5*T0 2.5*T0 -0.5 1.5])

Sfc=D0;
t=(-T0/2):0.0001:(T0/2);
for n=1:m
    Sfc=Sfc+Dn(-n)*exp(W0*-n*t*j)+Dn(n)*exp(W0*n*t*j);
end

subplot(3,2,2)
plot(t,f(t),'r','LineWidth',0.75)
grid on
hold on
plot(t,Sfc,'b','LineWidth',1.5)
legend('Función original','Serie de Fourier ','Location','northeast')
xlabel('t','FontWeight','bold','FontSize',16)
axis auto

subplot(3,2,4)
Ec=f(t)-Sfc;
plot(t,Ec,'LineWidth',2)
title('Error','FontWeight','bold','FontSize',16)
xlabel('t','FontWeight','bold','FontSize',16)
axis auto
grid on

subplot(3,2,6)
area(t,Ec.^2)
legend('Energia del error','Location','northeast')
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

subplot(3,2,5) % % 
stem(W0*nn,angle(absdn),'LineWidth',2) % % 
title('Espectro de fase, \angle de D_n ','FontWeight','bold','FontSize',16) % % 
xlabel('\omega','FontWeight','bold','FontSize',16)
grid on


%%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%% PR08 
%
% Agregar la expresión análitica de la serie en ambos formatos (hay dos series 
% en esta tarea).
%
% _ *1.-* _ Encuentra la expresion de la serie de Fourier, para la señal f(t) = t, 
% en el intervalo [−1, 1].Observa que es una señal impar.
%%

clearvars
figure(4)
hFigggg = figure(4);
set(hFigggg, 'Position', [0 0 900 900])
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
figure(5)
hFiggggg = figure(5);
set(hFiggggg, 'Position', [0 0 900 900])
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
figure(6)
hFigggggg = figure(6);
set(hFigggggg, 'Position', [0 0 900 900])

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


figure(7)
hFiggggggg = figure(7);
set(hFiggggggg, 'Position', [0 0 900 900])
f=@(u) u.^2;
u=-2:0.01:2;
plot(u,f(u))

%% 
% Realizando lo mismo que el anterior
% 
% a.-Serie de Fourier Trigonometrica Compacta
%%
syms n t

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
figure(8)
hFigggggggg = figure(8);
set(hFigggggggg, 'Position', [0 0 900 900])
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
figure(9)
hFiggggggggg = figure(9);
set(hFiggggggggg, 'Position', [0 0 900 900])
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
