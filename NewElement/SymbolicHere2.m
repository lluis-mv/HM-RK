%clc; 
clear all; close all;

% define some constants
syms h positive
syms M positive
syms k positive
syms dt positive
syms QBiot positive
syms AlphaStab real
% syms DeltaPW positive

% Define shape functions (to then integrate)
syms x positive
Np = [(h-x)/h; x/h];
Nu = inv(h^2/2)*[ (h/2-x)*(h-x);
    (-x)*(h/2-x)
    -2*(-x)*(h-x);
    ];
DN_DX = diff(Nu, x);

DN_DXp = diff(Np, x);
DN_DXp = DN_DX;
Np = Nu

% Define the elemental system


% Internal forces. Effective stress forces
K = +int( M *  DN_DX * (DN_DX'), x, 0,  h);
% Internal forces. Water pressure forces
Qt = -int( DN_DX * Np', x, 0,  h);


% Mass conservation
% Volume change
Q = int( Np*DN_DX', x, 0, h);
% Darcy law
H = int( DN_DXp*k*DN_DXp', x, 0, h);



Mstab = AlphaStab*h/12*[1, -1; -1 1];



Ce = [K, Qt; Q, zeros(3,3)];
Ke = [0*K, 0*Qt; 0*Q, H];
