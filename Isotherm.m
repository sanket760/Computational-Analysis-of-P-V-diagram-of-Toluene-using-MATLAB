Temp = 423;    % 150°C/423K.
Tc = 593;      % 320°C/593K.
Pc = 41;       % in bar...
R = 0.08314;   % for pressure in bar, volume in L, temperature in Kelvin....
T = 0;
for i = 1:20
    T = Temp + 10*(i-1) ;
    Tr = T/Tc;
    w = 0.263;                        % ascentric factor for toluene..
    c = 0.8792;                       % (0.48 + 1.564*w - 0.176*w^(2));
    b = 0.1042;                       % 0.08664*R*Tc/Pc;
    alpha = (1 + c*(1 - sqrt(Tr)))^2; % α_srk(Tr);
    a = 25.3645*alpha;                 % (0.42748*alpha*(R*Tc)^2)/Pc;
    
    hold on ;
    plot = fplot(@(v) (R*T/(v-b)) - a/(v*(v+b)), [0 10] ) ;
    xlabel 'Vm (in L/mol)';
    ylabel 'Pressure (in bar)';
    title 'P vs. Vm isotherms[T => 423-613K]';
    axis ([0 5 0 100]);
end