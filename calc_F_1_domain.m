function[Fcum] = calc_F_1_domain(data, Ea, lnDoa2, a)

% This function calculates the gas loss for a sphereical single domain with
% following inputs:

% data = two column array with step heat schedule 
%           col 1 = step number
%           col 2 = Temp (deg C)
%           col 3 = t (hr)
% Ea = activation energy in kcal/mol
% lnDoa2 = the y-intercept from your calculation of activation energy on a
%          arrhenius plot
% a = the relative size of the domain being modeled relative to ref_a
% below, which is tied to the Ea and lnDoa2 you input.

% This calcualtion of F uses the analytical solution for a sphere published
% by Farely et al., 2010 in G-cubed

time = data(:,3);
Temp = data(:,2);
Ea = Ea * 1000; % convert to cal from kcal
t = time * 3600; % time in seconds
TK = Temp + 273.15; % Temp in Kelvin
ref_a = 1; % radius of reference domain --> handle for given kinetics
% a = radius of domain being modeled --> relative to input kinetics
R = 1.987; % gas constant for Ea in kcal/mol

%calc 10000/K

for i = 1:length(TK)
    inv_TK(i) = 10000/TK(i);
end

% calc Do

lnDo = lnDoa2*ref_a^2; % note that this uses the referance domain size
Do = exp(lnDo);

% calc Dt

D = [];
Dt11 = [];
Dt11_cum = [];

Doa2 = Do/a^2;  % note that this uses the input size for this specific domain 
EaR = Ea/R;

for i = 1:length(t)
    D(i) = Doa2 * exp(-EaR * inv_TK(i)/10000);
    Dt11(i) = D(i)*t(i);
end

% calc cummulative Dt
Dt11_cum(1) = Dt11(1);

for i = 2:length(Dt11)
    Dt11_cum(i) = Dt11_cum(i-1)+Dt11(i);
end

% calculate Radial F for each step -- 

% NOTE that this code has the hooks for accomodating different 
% diffusivities for radial and axial pathways...
% this is not currently implemented, diffusion is isotropic here.
Dt33_cum = Dt11_cum;

F11 = [];
for i = 1:length(Dt11_cum)
    Dt = Dt11_cum(i);
    
    if Dt <= 0.2
        aa = 4*((Dt/pi)^0.5);
        bb = -1*(Dt^1.5/(3*(pi^0.5)));
        cc = -1*(0.2244122552*Dt^2);
        F11(i) = aa-Dt+bb+cc;
    else
        zz = exp(-5.783185963*Dt);
        F11(i) = 1-(0.691660276*zz);
    end
end

% calculate Axial F for each step

F33 = [];
for i = 1:length(Dt33_cum)
    Dt = Dt33_cum(i);
    
    if Dt <= 0.053
        F33(i) = 4*(Dt/pi)^0.5;
    else
        kk = exp(-1*(pi^2)*Dt);
        F33(i) = 1-(8*kk/(pi^2));
    end
end

% Now combine radial and axial components of F

Fcum = [];

for i = 1:length(F11)
    Fcum(i) = F11(i) + F33(i) - F11(i)*F33(i);
end

    
    

  



