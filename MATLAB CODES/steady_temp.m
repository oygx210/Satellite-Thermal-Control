function steady_temp(abs_sat1, emi_sat)


% defining constants
% all are in SI units
n = 3;
abs_p = 0.4;
emi_p = 0.6;
ref_p = 0.3;
mass_sat = 50;
tc_sat = 1000;
r_sat = 0.5;
H = 300000;
q_int = 100;
eff = 0;

if eff~=0
    pf = 0.3;
end


% Input data from excel sheet
X = xlsread('planets_data.xlsx');
distance = (X(:,1))';
mass = (X(:,2))';
temperature = (X(:,4))';
radius = (X(:,5))';


% Assigning corresponding data from excel sheet to new variables
temp_p = temperature(n);
mass_p = mass(n);
r_sp = distance(n);    %distance between sun and planet
r_p = radius(n);


%Defining constants
%All are in SI units
temp_s = 5778;
r_s = 6.9551e8;
sigma = 5.67e-8;  
abs_sat = 0;     

if eff~=0
      abs_sat = (abs_sat1 - eff * pf);
else
      abs_sat = abs_sat1;
end


%View factor calculations
vf_sat_s = 0.5 * (1-sqrt(1-((r_s^2)/(r_sp-r_p-H)^2)));
vf_s_sat = (vf_sat_s * (r_sat^2))/(r_s^2);
vf_sat_p = 0.5 * (1-sqrt(1-((r_p^2)/(H+r_p)^2)));
vf_p_sat = (vf_sat_p * (r_sat^2))/(r_p^2);
vf_sat_surr = 1;
vf_p_s = 0.5 * (1-sqrt(1-((r_s^2)/(r_sp^2))));
vf_s_p = (vf_p_s * (r_p^2))/(r_s^2);


%Planetary emission, solar and albedo radiations
q_s = (sigma * power(temp_s,4) * 4 * pi * power(r_s,2) * power(r_sat,2) * abs_sat) / (4 * power(r_sp,2)); %solar radiation
q_a = (sigma * power(temp_s,4) * 4 * pi * power(r_s,2) * vf_s_p * ref_p * vf_p_sat * abs_sat);            %albedo radiation
q_p = (sigma * emi_p * power(temp_p,4) * 4 * pi * power(r_p,2) * vf_p_sat * abs_sat );                    %planetary emission
q_sat0 = (sigma * emi_sat * 4 * pi * power(r_sat,2) * vf_sat_surr );                                      %(IR radiation emitted by satellite)/(T^4)


% Orbital time period of satellite calculation
t_o = 2 * pi * sqrt(power(r_p+H,3)/((6.67e-11)*mass_p));


% Albedo and eclipse functions in terms of satellite angular position
sympref('FloatingPointOutput',true); %Function used to convert symbolic display to decimal points
syms f_a f_e temp_mean_sat temp_amp_sat si phi 
f_e = (1 + cos(phi));   %cosine modulation over the average may be a suitable first approximation for albedo and eclipse functions
f_a = f_e;


% Temperature of satellite is assumed as approximate sinusoidal function for single node satellite
temp_sat = temp_mean_sat + temp_amp_sat * expand(cos(phi-si));


% Solving thermal balance equation for the satellite
eq1 = -expand((mass_sat * tc_sat * 2 * pi * temp_amp_sat * expand(sin(phi-si)))/t_o) ;       %Rate of energy stored in the body due to the thermal capacity of the satellite
eq2 = expand(q_s*f_e + q_a*f_a + q_p - expand(q_sat0*power(temp_mean_sat,3)*(temp_mean_sat+4*temp_amp_sat*expand(cos(phi-si))))) + q_int ; %Net heat interactions on the satellite


% Co-efficients of sin(?) and cos(?) are stored in co-efficient matrix along with the independent terms
c1 = coeffs(eq1, [sin(phi) cos(phi)]);        
c2 = coeffs(eq2, [sin(phi) cos(phi)]);
c0 = [0 c1(1) c1(2)];


% Equating co-efficient matrix gives phase lag(?) and Amplitude and mean temperatures of the satellite
[temp_mean_sat, temp_amp_sat, si] = solve(c0==c2, [temp_mean_sat, temp_amp_sat, si]);


%Printing output in a file
fileID = fopen('steadytemp.txt','a+');
fprintf(fileID,'Temperature of the satellite as a function of orbital position for absorptivity = %0.2f and emissivity = %0.2f is given by %.2f + %.2fcos(Ф-%.2f)\n', abs_sat1, emi_sat, abs(temp_mean_sat(1)), abs(temp_amp_sat(1)), abs(si(1)));
fclose(fileID);

end