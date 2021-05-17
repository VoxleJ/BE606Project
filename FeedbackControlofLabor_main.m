%% Normal Spontaneous Labor over 18 hours

clear all
close all

% Parameters from Gefferie et al 2018
P1 = 0.740; %mU/min Baseline oxytocin release rate in hypothalamus
P2 = 50; %mU/(min cm) Gain from cervical dilation to oxytocin release rate in hypothalamus
P3 = 0.0693; %1/min Pharmacokinetic parameter: Elimination rate constant
P4 = 18700; %mL Pharmacokinetic parameter: Volume of distribution
P5 = 0.5; %1/min Maximal concentration frequency
P6 = 7.9; %mU/mL Pharmacodynamic parameter: Oxytocin concentration resulting in 50% effect
P7 = 1.11; %dimensionless Pharmacodynamic parameter: Slope of sigmoidal curve
P8 = 40; %mm Hg Baseline contraction amplitude
P9 = 40; %mm Hg Maximal contraction amplitude
P10 = 0.001; %cm/min Dilation increase due to pressure exerted by the fetus on the cervix
P11 = 0.019; %cm/mm Hg Dilation increase due to contraction frequency and amplitude

%x(1) = d
%x(2) = m
%x(3) = b
%x(4) = g
%dm/dt = -P3 * m + (P1 + P2* d)
%dd/dt = P10 + P11* (P5 * (m/P4)^P7/(P6^P7 + (m/P4)^P7)) * (P8 +
%P9*(m/P4)^P7/(P6^P7 + (m/P4)^P7))


% define ODEs
g = @(t,x)[ P10 + P11 * (P5 * ((x(2) / P4)^P7 / (P6^P7+(x(2) / P4)^P7))) * (P8 + P9 * ((x(2) / P4)^P7 / (P6^P7+(x(2) / P4)^P7))); 
    -P3 * x(2) + (P1 + P2 * x(1));
    ((1/1+exp(-(Qv(x(1))+54.23)/9.88)) - x(3) )/(0.45 + 3.9/(1+((Qv(x(1)) + 66)/26)^2));
    ((0.02 + (0.98/(1+exp(-(Qv(x(1))+72.98)/4.64))))- x(4)) / (150 - 150/(1+exp((Qv(x(1))-417.43)/203.18)) * (1+exp(-(Qv(x(1))+61.11)/8.07)))];

%2cm dilation, 275mU Oxytocin
[t, xa] = ode45(@(t,x) g(t,x), [0 1200], [2,275,0,0]);

n = length(xa);

figure;
plot(linspace(0,18,n), xa(:,1))
title('Dilation')
xlabel('Time (hours)')
ylabel('Dilation (cm)')

figure;
plot(linspace(0,18,n), xa(:,2))
title('Oxytocin Mass')
xlabel('Time (hours)')
ylabel('Oxytocin (mU)')

% figure;
% plot(t, xa(:,3)) %sanity checks for b and g

% figure;
% plot(t, xa(:,4))

b2 = xa(:,3).^2;
V42 = Qv(xa(:,2)) - 42;
gg = xa(:,4);


Icat = 0.058 .* b2 .* gg .* V42;

figure;

plot(linspace(0,18,n), Icat)
title('T type Calcium Current')
xlabel('Time (hours)')
ylabel('Current (mA)')

% Ica = 0.058 .* xa(:,3).^2 .* xa(:,4).(Qv(xa(:,2)) - 42);

%% Stability Analysis on Normal 

% Phase Portrait Code is adapted from: 
% http://matlab.cheme.cmu.edu/2011/08/09/phase-portraits-of-a-system-of-odes/
d = linspace(0,12,100);
m = linspace(0,7700, 100);
b = linspace(0, 1.5, 100);
gg = linspace(-4, 8, 100); 


[D, M] = meshgrid(d,m);
[B, GG] = meshgrid(b,gg);
u = zeros(size(D));
v = zeros(size(D));
w = zeros(size(D));
z = zeros(size(D));

t=0;
for i = 1:numel(D)
    Y = g(t,[D(i); M(i); B(i); GG(i)]);
    u(i) = Y(1);
    v(i) = Y(2);
    w(i) = Y(3);
    z(i) = Y(4);
end


figure;

quiver(D,M,u,v,'r')
hold on

for Dilations = [2 3 4 5 6 7 8]
    [tphase,yphase] = ode45(@(t,x) g(t,x), [0 100], [Dilations,275,0,0]);
    plot(yphase(:,1),yphase(:,2))
    plot(yphase(1,1),yphase(1,2),'bo') % starting point
    plot(yphase(end,1),yphase(end,2),'ks') % ending point
end

title('Phase Portrait and Selected Solutions for Normal Labor')
xlabel('Dilation (cm)')
ylabel('Oxytocin (mU)')
hold off


%linearize system about equilibria
%solve for jacobian and equilibria using symbolic representations

syms m d b g V

dxdt = [ P10 + P11 * (P5 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))) * (P8 + P9 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))); 
        -P3 * m + (P1 + P2 * d);
        ((1/1+exp(-(V+54.23)/9.88)) - b )/(0.45 + 3.9/(1+((V + 66)/26)^2));
        ((0.02 + (0.98/(1+exp(-(V+72.98)/4.64))))- g) / (150 - 150/(1+exp((V-417.43)/203.18)) * (1+exp(-(V+61.11)/8.07)))];
    
    

J = jacobian(dxdt, [m d b g V]);

dxdt0 = [ P10 + P11 * (P5 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))) * (P8 + P9 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))) == 0; 
        -P3 * m + (P1 + P2 * d) == 0;
        ((1/1+exp(-(V+54.23)/9.88)) - b )/(0.45 + 3.9/(1+((V + 66)/26)^2)) == 0;
        ((0.02 + (0.98/(1+exp(-(V+72.98)/4.64))))- g) / (150 - 150/(1+exp((V-417.43)/203.18)) * (1+exp(-(V+61.11)/8.07))) == 0];
    
equilibria = solve(dxdt, [m d b g V]); 

Jvals = subs(J, equilibria);

eignvals = double(eig(Jvals));
disp('Normal Eigenvalues')
disp(eignvals)




disp('Normal Complete')
%% Induced Labour Simulation

clear g
%induced labour changes oxytocin ODE, introduction of synthetic
%oxytocin
%process begins with application of prostaglandin medications to the cervix
%Assumption: Cervix is ready; the general ODE of oxytocin should remain the
%same general parameters

%starting oxytocin, Penfield 2017, alternate high dose
%4mU/min -->dose increase: 4mU/min every 15 minutes
%certain protocols cap oxytocin adminstration at a max cumulative dose at
%10 U
%simulate this by adding (t/15*4) to dm/dt and starting at 279

g = @(t,x)[ P10 + P11 * (P5 * ((x(2) / P4)^P7 / (P6^P7+(x(2) / P4)^P7))) * (P8 + P9 * ((x(2) / P4)^P7 / (P6^P7+(x(2) / P4)^P7))); 
    -P3 * x(2)  + (P1 + P2 * x(1)) + t/15 * 4;
    ((1/1+exp(-(Qv(x(1))+54.23)/9.88)) - x(3) )/(0.45 + 3.9/(1+((Qv(x(1)) + 66)/26)^2));
    ((0.02 + (0.98/(1+exp(-(Qv(x(1))+72.98)/4.64))))- x(4)) / (150 - 150/(1+exp((Qv(x(1))-417.43)/203.18)) * (1+exp(-(Qv(x(1))+61.11)/8.07)))];

%2cm dilation, 279mU Oxytocin for Induced Labor
[t, xb] = ode45(@(t,x) g(t,x), [0 840], [2,279,0,0]);

n = length(xb);

figure;
plot(linspace(0,14,n), xb(:,1))
title('Dilation w/ Induced Labor')
xlabel('Time (hours)')
ylabel('Dilation (cm)')

figure;
plot(linspace(0,14,n), xb(:,2))
title('Oxytocin Mass w/ Induced Labor')
xlabel('Time (hours)')
ylabel('Oxytocin (mU)')

b2 = xb(:,3).^2;
V42 = Qv(xb(:,2)) - 42;
gg = xb(:,4);


Icat = 0.058 .* b2 .* gg .* V42;

figure;
plot(linspace(0,14,n), Icat)
title('T type Calcium Current w/ Induced Labor')
xlabel('Time (hours)')
ylabel('Current (mA)')


%% Stability Analysis of Induced Labor
%linearize system about equilibria
%solve for jacobian and equilibria using symbolic representations
%because there is only the addition of a time dependent constant, stability
%should be the same
clear m d b g V tt

syms m d b g V tt

dxdt1 = [ P10 + P11 * (P5 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))) * (P8 + P9 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))); 
        -P3 * m + (P1 + P2 * d) + tt/15 * 4;
        ((1/1+exp(-(V+54.23)/9.88)) - b )/(0.45 + 3.9/(1+((V + 66)/26)^2));
        ((0.02 + (0.98/(1+exp(-(V+72.98)/4.64))))- g) / (150 - 150/(1+exp((V-417.43)/203.18)) * (1+exp(-(V+61.11)/8.07)))];
    
    

J1 = jacobian(dxdt1, [m d b g V]);

dxdt2 = [ P10 + P11 * (P5 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))) * (P8 + P9 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))) == 0; 
        -P3 * m + (P1 + P2 * d) + tt/15 * 4 == 0;
        ((1/1+exp(-(V+54.23)/9.88)) - b )/(0.45 + 3.9/(1+((V + 66)/26)^2)) == 0;
        ((0.02 + (0.98/(1+exp(-(V+72.98)/4.64))))- g) / (150 - 150/(1+exp((V-417.43)/203.18)) * (1+exp(-(V+61.11)/8.07))) == 0];
    
equilibria1 = solve(dxdt2, [m d b g V]); 

Jvals1 = subs(J1, equilibria1);

eignvals1 = double(eig(Jvals1));
disp('Induced Labor Eigenvalues')
disp(eignvals1)



disp('Induced Complete')

%% Oxytocin Resistance Simulation
%change parameters such that there is a lower oxytocin receptor density
%Parameters to change--Gefferie mentioned P2 and P11 for patient specifics
clear g 
%we can assume a lower P2, in this case a decrease of half to exaggerate
P2 = 40;

g = @(t,x)[ P10 + P11 * (P5 * ((x(2) / P4)^P7 / (P6^P7+(x(2) / P4)^P7))) * (P8 + P9 * ((x(2) / P4)^P7 / (P6^P7+(x(2) / P4)^P7))); 
    -P3 * x(2) + (P1 + P2 * x(1));
    ((1/1+exp(-(Qv(x(1))+54.23)/9.88)) - x(3) )/(0.45 + 3.9/(1+((Qv(x(1)) + 66)/26)^2));
    ((0.02 + (0.98/(1+exp(-(Qv(x(1))+72.98)/4.64))))- x(4)) / (150 - 150/(1+exp((Qv(x(1))-417.43)/203.18)) * (1+exp(-(Qv(x(1))+61.11)/8.07)))];

%2cm dilation, 275mU Oxytocin
[t, xa] = ode45(@(t,x) g(t,x), [0 1200], [2,275,0,0]);

n = length(xa);

figure;
plot(linspace(0,18,n), xa(:,1))
title('Dilation w/ Decreased Oxytocin Receptors')
xlabel('Time (hours)')
ylabel('Dilation (cm)')

figure;
plot(linspace(0,18,n), xa(:,2))
title('Concentration w/ Decreased Oxytocin Receptors')
xlabel('Time (hours)')
ylabel('Oxytocin (mU)')

b2 = xa(:,3).^2;
V42 = Qv(xa(:,2)) - 42;
gg = xa(:,4);


Icat = 0.058 .* b2 .* gg .* V42;

figure;
plot(linspace(0,18,n), Icat)

title('T type Calcium Current w/ Decreased Oxytocin Receptors')
xlabel('Time (hours)')
ylabel('Current (mA)')

%% Stability Analysis on Oxytocin Resistant Labor
% The Stability of the System should not change, as only a single Parameter was changed
% Nevertheless for completeness:
syms m d b g V

dxdt4 = [ P10 + P11 * (P5 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))) * (P8 + P9 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))); 
        -P3 * m + (P1 + P2 * d);
        ((1/1+exp(-(V+54.23)/9.88)) - b )/(0.45 + 3.9/(1+((V + 66)/26)^2));
        ((0.02 + (0.98/(1+exp(-(V+72.98)/4.64))))- g) / (150 - 150/(1+exp((V-417.43)/203.18)) * (1+exp(-(V+61.11)/8.07)))];
    
    

J4 = jacobian(dxdt4, [m d b g V]);

dxdt5 = [ P10 + P11 * (P5 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))) * (P8 + P9 * ((m / P4)^P7 / (P6^P7+(m / P4)^P7))) == 0; 
        -P3 * m + (P1 + P2 * d) == 0;
        ((1/1+exp(-(V+54.23)/9.88)) - b )/(0.45 + 3.9/(1+((V + 66)/26)^2)) == 0;
        ((0.02 + (0.98/(1+exp(-(V+72.98)/4.64))))- g) / (150 - 150/(1+exp((V-417.43)/203.18)) * (1+exp(-(V+61.11)/8.07))) == 0];
    
equilibria3 = solve(dxdt5, [m d b g V]); 

Jvals5 = subs(J4, equilibria3);

eignvals5 = double(eig(Jvals5));
disp('Oxytocin Resistance Eigenvalues')
disp(eignvals5)

%Phase portraits
disp('Oxytocin Resistance Complete')



%% functions

%This function takes an input an through cosine, relates it to the 
%increasing voltage of the action potential
function [qt] = Qv(t) 
    if t <=1
    qt = -55* cos(t*25/2*pi);
    else    
    qt = 55* cos(t*25/2*pi);
    end
end

