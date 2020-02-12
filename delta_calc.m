function[delta] = delta_calc(data, Ea, lnDoa2)

% calculates ln(r/ro) which is referred to here as Delta
% It does so using an input Ea and (in kcal/mol!) and ln(Do/a2) 


lnDa2 = data(:,2);
TK = data(:,1);
R = 1.987;  % gas constant for kcal/mol
slope = -Ea/(R*10);
Do = lnDoa2;

for i = 1:length(lnDa2)
    b = slope*TK(i);
    c = b+Do;
    d = c - lnDa2(i);
    
    delta(i) = d*0.5;
end



