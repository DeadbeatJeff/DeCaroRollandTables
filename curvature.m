% Only works with all revolute joint robots, currently; will need to generalize

clear all; clc; close all;

syms g;
% g = 9.81

% dhConvention = [theta alpha r d]

dhParams = [0   0	1  	0;
            0	0   0   1];
            % 0.0203	-pi/2	0.15005	0];
            % 0   	pi/2	0.4318	0;
            % 0       -pi/2	0   	0;
            % 0       0       0       0];

T02 = vpa(DH_HTM(dhParams,'r'),3) ;

% dhConvention = [theta alpha r d m Ixx Ixy Ixz Iyy Iyz Izz] in MKS

[nJeff ~] = size(dhParams);

Theta=sym('Theta',[1,nJeff]);
syms(Theta);

thetaDot=sym('thetaDot',[1,nJeff]);
syms(thetaDot);

doubleThetaDot=sym('doubleThetaDot',[1,nJeff]);
syms(doubleThetaDot);

m=sym('m',[1,nJeff]);
syms(m);

Ixx=sym('Ixx',[1,nJeff]);
syms(Ixx);

Ixy=sym('Ixy',[1,nJeff]);
syms(Ixy);

Ixz=sym('Ixz',[1,nJeff]);
syms(Ixz);

Iyy=sym('Iyy',[1,nJeff]);
syms(Iyy);

Iyz=sym('Iyz',[1,nJeff]);
syms(Iyz);

Izz=sym('Izz',[1,nJeff]);
syms(Izz);

L=sym('L',[1,nJeff]);
syms(L);

for i = 1:nJeff
    m(i) = 0.181;
    Ixx(i) = 0.093805;
    Ixy(i) = 0;
    Txz(i) = 0;
    Iyy(i) = 0.001;
    Iyz(i) = 0;
    Izz(i) = 0.001;
    L(i) = 0.250;
end 

drParamsSym = [Theta1   0	L1  	0 m(1) Ixx(1) Ixy(1) Ixz(1) Iyy(1) Iyz(1) Izz(1);
               Theta2	0   0   L2 m(2) Ixx(2) Ixy(2) Ixz(2) Iyy(2) Iyz(2) Izz(2)];
               % 0.0203	-pi/2	0.15005	0];
               % 0   	pi/2	0.4318	0;
               % 0       -pi/2	0   	0;
               % 0       0       0       0];  

RHS=sym('RHS',[1,nJeff]);
syms(RHS);

x=sym('x',[1,nJeff]);
syms(x);

y=sym('y',[1,nJeff]);
syms(y);

Tau=sym('Tau',[1,nJeff]);
syms(Tau);

% L = ((1/6)*drParamsSym(1,5)*drParamsSym(1,3)^2 + (1/6)*drParamsSym(2,5)*drParamsSym(2,4)^2 + (1/2)*drParamsSym(2,5)*drParamsSym(1,3)^2 + (1/2)*drParamsSym(1,3)*drParamsSym(1,3)*drParamsSym(2,4)*cos(Theta2))*Theta1^2 + (1/6)(1/6)*drParamsSym(2,5)*drParamsSym(2,4)^2*Theta2^2

x1 = (1/2)*drParamsSym(1,3)*cos(drParamsSym(1,1));
y1 = (1/2)*drParamsSym(1,3)*sin(drParamsSym(1,1));

syms thetaSum;

thetaSum = drParamsSym(1,1);

x(1) = drParamsSym(1,3)*cos(drParamsSym(1,1));
y(1) = drParamsSym(1,3)*sin(drParamsSym(1,1));
for i = 1:nJeff
    xSum = x(1);
    ySum = y(1);
    thetaSum = drParamsSym(1,1);
    for j = 2:i-1
        thetaSum = thetaSum + drParamsSym(j,1);
        xSum = xSum + drParamsSym(j,4)*cos(thetaSum);
        ySum = ySum + drParamsSym(j,4)*sin(thetaSum);
    end
    for j = i:i
        thetaSum = thetaSum + drParamsSym(j,1);
        xSum = xSum + (1/2)*drParamsSym(j,4)*cos(thetaSum);
        ySum = ySum + (1/2)*drParamsSym(j,4)*sin(thetaSum);
    end
    x(i) = xSum;
    y(i) = ySum;
end

xdot=sym('xdot',[1,nJeff]);
syms(xdot);

ydot=sym('ydot',[1,nJeff]);
syms(ydot);

syms VarJeff;

xdot(1) = diff(x1,drParamsSym(1,1))*thetaDot1;
ydot(1) = diff(y1,drParamsSym(1,1))*thetaDot1;

for i = 2:nJeff
    xdotSum = VarJeff;
    ydotSum = VarJeff;
    for j = 1:i
        xdotSum = xdotSum + diff(x(i),drParamsSym(j,1))*thetaDot(j);
        ydotSum = ydotSum + diff(y(i),drParamsSym(j,1))*thetaDot(j);
    end
    xdot(i) = xdotSum - VarJeff;
    ydot(i) = ydotSum - VarJeff;
end

syms KE V VTerm termSum;

KETerm=sym('KETerm',[1,nJeff]);
syms(KETerm);

thetaDotSum = thetaDot(1);
v(1) = sqrt(xdot(1)^2 + ydot(1)^2);
KETerm(1) = (1/2)*drParamsSym(1,5)*v(1)^2 + (1/2)*drParamsSym(1,11)*thetaDot(1)^2;
KE = KETerm(1);
for i = 2:nJeff
    thetaDotSum = thetaDotSum + thetaDot(i);
    v(i) = sqrt(xdot(i)^2 + ydot(i)^2);
    KETerm(i) = (1/2)*drParamsSym(i,5)*v(i)^2 + (1/2)*drParamsSym(i,11)*thetaDotSum^2;
    KE = KE + KETerm(i);
end
% KE = KE - VarJeff;

VTerm=sym('VTerm',[1,nJeff]);
syms(VTerm);

V = VarJeff;
VTerm(1) = drParamsSym(1,5)*g*(1/2)*drParamsSym(1,3)*sin(drParamsSym(1,1));
V = VTerm(1);
for i = 2:nJeff
    thetaSum = drParamsSym(1,1);
    VTerm(i) = VarJeff;
    VTerm(i) = VTerm(i) + drParamsSym(i,5)*g*drParamsSym(1,3)*sin(thetaSum);
    for j = 2:i-1
        thetaSum = thetaSum + drParamsSym(j,1);
        VTerm(i) = VTerm(i) + drParamsSym(i,5)*g*drParamsSym(j,4)*sin(thetaSum);
    end
    thetaSum = thetaSum + drParamsSym(i,1);
    VTerm(i) = VTerm(i) + (1/2)*drParamsSym(i,5)*g*drParamsSym(i,4)*sin(thetaSum);
    VTerm(i) = VTerm(i) - VarJeff;
    V = V + VTerm(i);
end
V = V - VarJeff;

Lagrangian = KE - V;

Lagrangian1=sym('Lagrangian1',[1,nJeff]);
syms(Lagrangian1)

for i = 1:nJeff
    Lagrangian1(i) = diff(Lagrangian,thetaDot(i));
end

syms varJeff

for i = 1:nJeff
    RHS(i) = varJeff;
    for j = 1:nJeff
        RHS(i) = RHS(i) + diff(Lagrangian1(i),Theta(j))*thetaDot(j) + diff(Lagrangian1(i),thetaDot(j))*doubleThetaDot(j);
    end;
    RHS(i) = RHS(i) - diff(Lagrangian,Theta(i)) - varJeff;
end

MJeffVec=sym('MJeff',[1,nJeff^2]);
syms(MJeffVec);
MJeff = reshape(MJeffVec,[nJeff nJeff]).';

for i = 1:nJeff
    for j = 1:nJeff
PlaceholderJeff = fliplr(coeffs(RHS(i),doubleThetaDot(j)));
MJeff(i,j) = PlaceholderJeff(1);
    end
end

MJeff; % I have the Riemannian metric

MJeffinv = inv(MJeff);

matlabFunction(MJeffinv, 'File', 'Minvmatrix');

MmJeff = MJeff*doubleThetaDot.';

Gamma1st = sym(zeros(nJeff,nJeff,nJeff));

for i=1:nJeff
    for j = 1:nJeff
        for k = 1:nJeff
            Gamma1st(i,j,k) = (1/2)*(diff(MJeff(i,k),Theta(j))+diff(MJeff(j,k),Theta(i))-diff(MJeff(i,j),Theta(k)));
        end
    end
end

Gamma1st = simplify(Gamma1st);

Gamma2nd = sym(zeros(nJeff,nJeff,nJeff));

for i=1:nJeff
    for j = 1:nJeff
        for l = 1:nJeff
            Gamma2nd(i,j,l) = VarJeff;
            for k = 1:nJeff
                Gamma2ndTerm = (1/2)*MJeffinv(l,k)*(diff(MJeff(i,k),Theta(j))+diff(MJeff(j,k),Theta(i))-diff(MJeff(i,j),Theta(k)));
                Gamma2nd(i,j,l) = Gamma2nd(i,j,l) + Gamma2ndTerm;
            end
            Gamma2nd(i,j,l) = Gamma2nd(i,j,l) - VarJeff;
        end
    end
end

Gamma1st = simplify(Gamma1st);

matlabFunction(Gamma1st, 'File', 'Gammamatrix');

VmJeff=sym('VmJeff',[nJeff,1]);
syms(VmJeff);
for i = 1:nJeff
    VmJeff(i,1) = simplify(thetaDot*Gamma1st(:,:,i)*thetaDot.');
end

G1= fliplr(coeffs(RHS(1), g));
G2 =fliplr(coeffs(RHS(2), g));

GmJeff = -[G1(1); G2(1)]*g;

matlabFunction(GmJeff, 'File', 'Gvector');

% EoM: Tau = MmJeff + VmJeff + GmJeff
% EoM: doubleThetaDot = MJeffInv*(Tau - VmJeff - GmJeff)

Riemann = sym(zeros(nJeff,nJeff,nJeff,nJeff));
syms GammaSum;

for i=1:nJeff
    for j = 1:nJeff
        for k = 1:nJeff
            % RiemannSum = VarJeff;
            for l = 1:nJeff
                Riemann(i,j,k,l) = diff(Gamma1st(i,k,l),Theta(j)) - diff(Gamma1st(j,k,l),Theta(i));
                GammaSum = VarJeff;
                for p = 1:nJeff
                    GammaSum = GammaSum + Gamma1st(i,k,p)*Gamma1st(j,p,l) - Gamma1st(j,k,p)*Gamma1st(i,p,l);
                end
                GammaSum = GammaSum - VarJeff;
                Riemann(i,j,k,l) = Riemann(i,j,k,l) + GammaSum;
            end
        end
    end
end % I have the Riemann tensor

read = eye(nJeff);

K = sym(zeros(nJeff,nJeff));

Vec1 = sym(zeros(nJeff,nJeff,nJeff));

for i = 1:nJeff
    for j = 1:nJeff
        Vec1(:,i,j) = VarJeff*ones(nJeff,1);
        for p = 1:nJeff
            Vec1(:,i,j) = Vec1(:,i,j) + Riemann(p,i,j,j)*read(:,p:p);
        end
        Vec1(:,i,j) = Vec1(:,i,j)-VarJeff*ones(nJeff,1);
        Vec1(:,i,j);
        K(i,j) = read(:,i:i).'*Vec1(:,i,j);
    end
end % I have the sectional curvature

KCSpace = K(1,2)



% sum 1 2 2
