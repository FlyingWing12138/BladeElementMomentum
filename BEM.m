clear
clc
close all

Nb = 4;                 % Number of blades
Clalpha = 5.73;         % NACA 0012, [rad^-1]
HingeOffset = 0.0508;   % [m]
RootCutout = 0.2096;    % [m]
thetatw = -8*2*pi/360;  % Linear Twist [rad]
R = 0.8606;             % Radius [m]
c = 0.0660;             % Airfoil chord [m]
sigma = 0.0977;         % Rotor solidity, Nb*c/(pi*R)
rho = 1.293;            % Air density kg/m^3

Cd = 0.0002;            % Rotor drag coefficient
alphas = -3*2*pi/360;   % Angle between rotor disk and free-stream velcity [rad]
muinf = 0.149;          % rotor advance ratio, Uinf*cos(alphas)/Vtip
Ct = 0.0063;            % Rotor thrust coefficient, T/(rho*A*Vtip^2)
Cq = 0.00036;           % Rotor torque coefficient, Q/(rho*A*R*Vtip^2)
Vinf = 28.5002;         % [m/s]
Vtip = 190.2866;        % Rotor blade tip velocity [m/s]
Mtip = 0.5533;          %
% mu = Vinf*cos(alphas)/(Omega*R);
Omega = Vtip/R;           % Rotor speed [rpm]

psi = 0:2*pi/144:2*pi;
r = (0.2:0.01:1)';

% theta0 = 9.37;          % Collective pitch at r/R=0.75 [deg]
% theta1c = -1.11;        % Lateral pitch (cosine term) [deg]
% theat1s = 3.23;         % Longitudinal pitch (sine term) [deg]
% theta0 = 6.26;          % Collective pitch at r/R=0.75 [deg]
% theta1c = -2.08;        % Lateral pitch (cosine term) [deg]
% theat1s = 1.96;         % Longitudinal pitch (sine term) [deg]
theta0 = 5.49;              % Collective pitch at r/R=0.75 [deg]
theta1c = -2.23;            % Lateral pitch (cosine term) [deg]
theat1s = 2.09;             % Longitudinal pitch (sine term) [deg]
% theta0 = 6.873;          % Collective pitch at r/R=0.75 [deg]
% theta1c = -2.014;        % Lateral pitch (cosine term) [deg]
% theat1s = 1.073;  
theta75 = (theta0-theta1c*cos(psi)-theat1s*sin(psi))*2*pi/360;    % Pitch control [rad]
thetar=theta75+(r-0.75).*thetatw;

beta0 = 1.5;            % Coning angle [deg]
beta1c = 0;             
beta1s = 0;
beta = (beta0+beta1c*cos(psi)+beta1s*cos(psi))*2*pi/360;        % Flapping [rad]
betadot = 0;

modelCell = {'Coleman','Drees','Payne','WhiteBlake','PittPeters','Howlett'};
figureVisibility = 'on';

%% Linear Inflow Model Plotting

for i = 1:size(modelCell,2)

    model = modelCell{i};
    lambda = linearInflow(muinf, alphas, Ct, psi, r, model);

    Dynamics = bladeElement(r, R, Omega, muinf, psi, lambda, beta, betadot, thetar, rho, c, Clalpha, Cd, Nb);

    a=figure('Visible',figureVisibility);
    polarPcolor(r',psi*180/pi,lambda,'Ncircles',5,'Nspokes',13,'typeRose','default')
    % sgtitle('Linear Inflow')
    view([-270, 90])
    clim([-0.015 0.075])
    a.Children(1).Ticks=-0.01:0.01:0.07;

    inflowTitle=[model,'LinearInflow'];
%     saveas(a, inflowTitle, 'fig')
%     saveas(a, inflowTitle, 'png')

    b=figure('Visible',figureVisibility);
    polarPcolor(r',psi*180/pi,Dynamics.dCt,'Ncircles',5,'Nspokes',13,'typeRose','default')
    % sgtitle('Linear Inflow Thrust Coefficient')
    view([-270, 90])
    clim([-1e-5 21e-5])
    b.Children(1).Ticks=0:2e-5:2e-4;

    thrustTitle = [model,'InflowThrust'];
%     saveas(b, thrustTitle, 'fig')
%     saveas(b, thrustTitle, 'png')

    cc=figure('Visible',figureVisibility);
    polarPcolor(r',psi*180/pi,Dynamics.alpha/2/pi*360,'Ncircles',5,'Nspokes',13,'typeRose','default')
    % sgtitle('Linear Inflow Angle of Attack')
    view([-270, 90])
    clim([-4 7.5])
    cc.Children(1).Ticks=-4:1:7;

    alphaTitle = [model,'InflowAlpha'];
%     saveas(cc, alphaTitle, 'fig')
%     saveas(cc, alphaTitle, 'png')

    fileName = [model,'Dynamics'];
%     save(fileName, "Dynamics", "lambda", "r", "psi");

    fprintf('%s Model Complete. \n', model)

end

%% ThetaR Configuration Plotting

if 1 == 1

    d=figure('Visible',figureVisibility);
    polarPcolor(r',psi*180/pi,thetar/2/pi*360,'Ncircles',5,'Nspokes',9,'typeRose','default')
    % sgtitle('Theta Configuration')
    view([-270, 90])
    clim([1.5 13.5])
    d.Children(1).Ticks=2:1:13;

%     saveas(d,'ThetaConfiguration','fig')
%     saveas(d,'ThetaConfiguration','png')

    fprintf('Theta Configuration Complete. \n')

end

%% Uniform Inflow Model Plotting

if 1 == 1

    lambda = uniformInflow(muinf, alphas, Ct);
    lambda = repmat(lambda, size((r.*psi)));

    Dynamics = bladeElement(r, R, Omega, muinf, psi, lambda, beta, betadot, thetar, rho, c, Clalpha, Cd, Nb);

    e=figure('Visible',figureVisibility);
    polarPcolor(r',psi*180/pi,lambda,'Ncircles',5,'Nspokes',9,'typeRose','default')
    % sgtitle('Uniform Inflow')
    view([-270, 90])
    clim([-0.015 0.075])
    e.Children(1).Ticks=-0.01:0.01:0.07;

%     saveas(e, 'UniformInflow','fig')
%     saveas(e, 'UniformInflow','png')

    f=figure('Visible',figureVisibility);
    polarPcolor(r',psi*180/pi,Dynamics.dCt,'Ncircles',5,'Nspokes',9,'typeRose','default')
    % sgtitle('Uniform Inflow Thrust Coefficient')
    view([-270, 90])
    clim([-1e-5 21e-5])
    f.Children(1).Ticks=0:2e-5:2e-4;

%     saveas(f,'UniformInflowThrust','fig')
%     saveas(f,'UniformInflowThrust','png')

    g=figure('Visible',figureVisibility);
    polarPcolor(r',psi*180/pi,Dynamics.alpha/2/pi*360,'Ncircles',5,'Nspokes',9,'typeRose','default')
    % sgtitle('Uniform Inflow Thrust Coefficient')
    view([-270, 90])
    clim([-4 7.5])
    g.Children(1).Ticks=-4:1:7;

%     saveas(g,'UniformInflowAlpha','fig')
%     saveas(g,'UniformInflowAlpha','png')

%     save("UniformDynamics", "Dynamics", "lambda", "r", "psi");

    fprintf('Uniform Model Complete. \n')

end

%% Momentum Theory Inflow Model

function lambda = uniformInflow(mu, alpha, Ct)

lambda0 = sqrt(Ct/2);
epsilon0 = 0.0005;
epsilon = 1;
i = 0;
while epsilon>=epsilon0

    lambda = mu*tan(-alpha)+Ct/(2*sqrt(mu^2+lambda0^2));
    epsilon = norm((lambda-lambda0)/lambda);

    lambda0 = lambda;
    i = i+1;
end

end

function lambda = linearInflow(mu, alpha, Ct, psi, r, model)

lambda0 = uniformInflow(mu, alpha, Ct);

chi = atan(mu/lambda0);

switch model
    case 'Coleman'
        kx = tan(chi/2);
        ky = 0;

    case 'Drees'
        kx = (4/3)*(1-cos(chi)-1.8*mu^2)/sin(chi);
        ky = -2*mu;

    case 'Payne'
        kx = (4/3)*(mu/lambda0/(1.2+mu/lambda0));
        ky = 0;

    case 'WhiteBlake'
        kx = sqrt(2)*sin(chi);
        ky = 0;

    case 'PittPeters'
        kx = (15*pi/23)*tan(chi/2);
        ky = 0;

    case 'Howlett'
        kx = sin(chi)^2;
        ky = 0;

    otherwise
        fprintf('Input Correct Linear Inflow Model\n')

end

lambda = lambda0*(1+kx*r.*cos(psi)+ky*r.*sin(psi));

end

%% Blade Element Theory

function Dynamics = bladeElement(r, R, Omega, muinf, psi, lambda, beta, betadot, thetar, rho, c, Clalpha, Cd, Nb)

Dynamics = struct;

y = r*R;
% dr = r(2)-r(1);
dy = y(2)-y(1);
A = pi*R^2;
dpsi = psi(2)-psi(1);

Ut = Omega*y+muinf*Omega*R*sin(psi);
Up = lambda*Omega*R+y*betadot+muinf*Omega*R*beta.*cos(psi);
U = sqrt(Ut.^2+Up.^2);

phi = atan(Up./Ut);
alpha = thetar - phi;

dL = 0.5*rho*U.^2.*c.*Clalpha.*alpha.*dy;
dD = 0.5*rho*U.^2.*c.*Cd.*dy;
dFz = dL.*cos(phi)-dD.*sin(phi);
dFx = dL.*sin(phi)+dD.*cos(phi);
dT = Nb*dFz;
dQ = Nb*dFx.*y;
dP = Nb*dFx*Omega.*y;

F = prandtlTipLoss(Nb, r, lambda);

dCt=F.*dT/(rho*A*(Omega*R)^2); 
dCq=F.*dQ/(rho*A*(Omega*R)^2*R);
dCp=F.*dP/(rho*A*(Omega*R)^3);

Ct = sum(sum(dCt.*r*dpsi))/pi;
Cq = sum(sum(dCq.*r*dpsi))/pi;
Cp = sum(sum(dCp.*r*dpsi))/pi;

Dynamics.Ctright = sum(sum(dCt(:,1:floor(size(dCt,2)/2)).*r*psi(2)))/pi;
Dynamics.Ctleft = sum(sum(dCt(:,ceil(size(dCt,2)/2):end).*r*psi(2)))/pi;

Dynamics.Ct = Ct;
Dynamics.Cq = Cq;
Dynamics.Cp = Cp;

Dynamics.dCt = dCt;
Dynamics.dCq = dCq;
Dynamics.dCp = dCp;

Dynamics.alpha = alpha;

end

function F = prandtlTipLoss(Nb, r, lambda)

phi=lambda./r;
% phi(phi(:,:)<0)=0;
f=abs(Nb/2*(1-r)./(r.*phi));
F=(2/pi)*acos(exp(-f));

end