function[arrDATA] = arrhenius_calc(data_in)

%  This function calculates the arrhenius data for an initially 
%  square profile using the equations from F+K 19XX

% The inputs for this function are a three column matrix of step heating 
% data with the following:
%   col 1 = step number
%   col 2 = step Temp (deg C)
%   col 3 = step duration (hr)
%   col 4 = step F (0-1) not %

% This function returns a two column array with 10000/K and LnDa^2 for each
% step.

time_s = data_in(:,3)*3600;
TK = data_in(:,2) + 273.15;
Fcum = data_in(:,4);
arrDATA =[];
lnDa2 =[];

% trap out Fcum = 1

for i = 1:length(Fcum)
    if Fcum(i)==1
        Fcum(i) = 0.999999;
    end
end


% Square Profile (1st step)

F = Fcum(1);
t = time_s(1);

if F < 0.1
    lnDa2(1) = log((F^2*pi)/(36*t)); % NOTE this is natural log!!  

elseif F > 0.9
    a = 1/((pi^2)*t);
    b = log((pi^2/6)*(1-F));
    lnDa2(1) = log(-a*b);
    
else
    a = 1/((pi^2)*t);
    b = (2*pi)-((pi^2/3)*F);
    c = (2*pi)*(sqrt(1-(pi/3))*F);
    lnDa2(1) = log(a*(b-c));
end

% Rounded Profile -- rest of steps

for i = 2:length(Fcum)
    
    F = Fcum(i);
    Fbefore = Fcum(i-1);
    t = time_s(i);
    
    if F < 0.1
        lnDa2(i) = log((((F^2)-(Fbefore^2))*pi)/(36*t));

    elseif F > 0.9
        lnDa2(i) = log((1/((pi^2)*t))*log((1-Fbefore)/(1-F)));
    
    else
        a = 1/((pi^2)*t);
        b = (-1*(pi^2/3))*(F-Fbefore);
        c = (2*pi)*((sqrt(1-((pi/3)*F)))-(sqrt(1-((pi/3)*Fbefore))));
        lnDa2(i) = log(a*(b-c));
    end
end
    

% calc 10000/T
for i = 1:length(Fcum)
    arrDATA(i,1) = 10000/TK(i);
end

arrDATA(:,2) = lnDa2;



