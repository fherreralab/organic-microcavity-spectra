clc;
clear all;
%% Constantes
c_speed = 3e17;% nm
hbar = 6.582119514e-16; %eV*s
J = 2.6544e-3;% factor admitancia

Cavidad = {'sustrato','plata','silica','TDBC','silica','plata','aire'};
Ncapas = length(Cavidad);
polarizacion = 2; %TM = 1, TE = 2

dEntrada = 0;
dAg = 30;
dTDBC = 20;
dSilicaL = 50;
dSilicaR = 50;
dS = 0;

n0 = 1;
nEntrada = 1.4587;
nS = 1;

%% Definición de parámetros
x = 200;
y = 200;
Emin = 1.5; % eV
Emax = 3.2; % eV
lambda = linspace(2*pi*c_speed*hbar/Emax,2*pi*c_speed*hbar/Emin,y);

theta = (pi/180)*linspace(0,90,x);
w1 = 2*pi*c_speed./lambda;
E = hbar*w1;
ko = 2*pi./lambda;

kx = zeros(y,x);
R = zeros(length(theta),length(lambda));
T = zeros(length(theta),length(lambda));
A = zeros(length(theta),length(lambda));

for i = 1:y
    for ii = 1:x
kx(i,ii) = n0*ko(i)*sin(theta(ii));
    end
end

for k = 1:length(theta)
for l = 1:length(lambda)
%% Variables
w = 2*pi*c_speed/lambda(l);

eps5a = lorentzAg(w); % silver
eps5a_real = real(eps5a);
eps5a_im = imag(eps5a);

eps5g = lorentzTDBC(w); % TDBC
eps5g_real = real(eps5g);
eps5g_im = imag(eps5g);

% Indice de refracción Ag y TBDC
nAg = (eps5a_real+1i*eps5a_im)^0.5;
nTDBC = (eps5g_real+1i*eps5g_im)^0.5;
n_silica = sellmeier(lambda(l));

theta_Entrada = theta(k);
theta_Ag1 = asin((nEntrada/nAg)*sin(theta_Entrada));
theta_silicaL = asin((nAg/n_silica)*sin(theta_Ag1));
theta_TDBC = asin((n_silica/nTDBC)*sin(theta_silicaL));
theta_silicaR = asin((nTDBC/n_silica)*sin(theta_TDBC));
theta_Ag2 = asin((n_silica/nAg)*sin(theta_silicaR));
theta_S = asin((nAg/nS)*sin(theta_Ag2));

deltaEntrada = 2*pi*nEntrada*dEntrada*cos(theta_Entrada)/lambda(l);
deltaAg1 = 2*pi*nAg*dAg*cos(theta_Ag1)/lambda(l);
delta_silicaL = 2*pi*n_silica*dSilicaL*cos(theta_silicaL)/lambda(l);
deltaTDBC = 2*pi*nTDBC*dTDBC*cos(theta_TDBC)/lambda(l);
delta_silicaR = 2*pi*n_silica*dSilicaR*cos(theta_silicaR)/lambda(l);
deltaAg2 = 2*pi*nAg*dAg*cos(theta_Ag2)/lambda(l);
deltaS = 2*pi*nS*dS*cos(theta_S)/lambda(l);

if polarizacion == 1;
etaEntrada = nEntrada/cos(theta_Entrada);
etaAg1 = nAg/cos(theta_Ag1);
eta_silicaL = n_silica/cos(theta_silicaL);
etaTDBC = nTDBC/cos(theta_TDBC);
eta_silicaR = n_silica/cos(theta_silicaR);
etaAg2 = nAg/cos(theta_Ag2);
etaS = nS/cos(theta_S);
eta_m = etaS;
eta_o = n0/cos(theta(k));
else
etaEntrada = nEntrada*cos(theta_Entrada);
etaAg1 = nAg*cos(theta_Ag1);
eta_silicaL = n_silica*cos(theta_silicaL);
etaTDBC = nTDBC*cos(theta_TDBC);
eta_silicaR = n_silica*cos(theta_silicaR);
etaAg2 = nAg*cos(theta_Ag2);
etaS = nS*cos(theta_S);
eta_m = etaS;
eta_o = n0*cos(theta(k));
end

delta = [deltaEntrada deltaAg1 delta_silicaL deltaTDBC delta_silicaR deltaAg2 deltaS];
eta = [etaEntrada etaAg1 eta_silicaL etaTDBC eta_silicaR etaAg2 etaS];

%% Matriz de matrices de cada capa
Matrices = [cos(deltaEntrada),-1i*sin(deltaEntrada)/etaEntrada; % primera matriz
           -1i*etaEntrada*sin(deltaEntrada),cos(deltaEntrada)];

for i=2:Ncapas
    Matrices(:,:,i) = [cos(delta(i)),-1i*sin(delta(i))/eta(i); % matrices de cada capa
                       -1i*eta(i)*sin(delta(i)),cos(delta(i))];
end

%% Matriz del conjunto
Mcon = Matrices(:,:,1);

for i=2:Ncapas
    Mcon = Mcon*Matrices(:,:,i); % matriz del conjunto
end

%% Cálculo de la reflectancia
BC = Mcon*[1;etaS];
Y = BC(2)/BC(1);
R(k,l) = ((eta_o-Y)/(eta_o+Y))*conj((eta_o-Y)/(eta_o+Y)); % reflectancia
A(k,l) = (4*eta_o*real(BC(1)*conj(BC(2))-eta_m))/((eta_o*BC(1)+BC(2))*conj(eta_o*BC(1)+BC(2))); % absorbancia
T(k,l) = (4*eta_o*real(eta_m))/((eta_o*BC(1)+BC(2))*conj(eta_o*BC(1)+BC(2))); % transmitancia

end
end

%% Curva de Dispersión
% figure (1)
% E = repmat(E,x,1);
% h = surf(kx',E,R);
% colormap (parula)
% ylim([Emin,Emax])
% xlim([0,max(max(kx))+0.05*max(max(kx))])
% set(h,'edgecolor','none');
% colorbar
% view(2)
% ylabel('Energía [eV]')
% xlabel('k_{||} [nm^{-1}]')

%% Reflectancia, Absortancia y Transmitancia
figure (2)
plot(lambda,T(100,:),'b',lambda,R(100,:),'r',lambda,A(100,:),'g')
legend('\fontsize{9}Transmitancia','\fontsize{9}Reflectancia','\fontsize{9}Absortancia')
xlabel('Longitud de onda [nm]')
ylabel('Intensidad')
