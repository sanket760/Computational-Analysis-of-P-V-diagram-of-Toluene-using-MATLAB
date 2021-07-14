Temp = 423;   % 150°C/423K.
tol = 0.001;  % tolerance limit for difference in chemical potentials.
Tc = 593;     % critical temperature from NIST database(in Kelvin).
Pc = 41;      % critical pressure from NIST database(in Bar).
A = 4.54436;  % Antoine equation parameters from NIST database.
B = 1738.123; % all these parameters are valid over T = [420,580].
C = 0.394;
R = 0.08314;                      % Gas constant in bar.L/mol.K.
t = 423:10:593;                   % Temperature values
p = 7:2:41;
w = 0.263;                        % ascentric factor for toluene..
c = 0.8792;                       % (0.48 + 1.564*w - 0.176*w^(2));
b = 0.1042;                       % 0.08664*R*Tc/Pc;
Ps = zeros(17,1);
v_l = zeros(17,1);
v_g = zeros(17,1);

for i = 1:17
    i
    t(i)
    Tr = t(i)/Tc;            % reduced temp.
    Pg = 10^(A - B/(t(i) + C))
    dif_meu = 10;            % random initialisation 
    a = 0.42748*((1 + c*(1 - sqrt(Tr)))*(R*Tc))^(2)/Pc
    c0 = -a*b;
    c1 = -(b^2 + R*t(i)*b - a);
    c2 = -R*t(i);
    c3 = Pg;
    x = roots([c3 c2 c1 c0]);
    y = x(imag(x) == 0);
    l = length(y);
    
    v_min = min(y);
    d_min = 1/v_min; % density corresponding to lowest volume root....
    v_max = max(y);
    d_max = 1/v_max; % density correaponding to highest volume root...
        
    z_l = (1/(1-d_min*b))-((a*d_min)/(R*t(i)*(1+d_min*b)));     
    z_g = (1/(1-d_max*b))-((a*d_max)/(R*t(i)*(1+d_max*b)));
    meu_l = -log(1-d_min*b)-((a)/(b*R*t(i)))*log(1+d_min*b)+z_l-1-log(z_l);  
    meu_g = -log(1-d_max*b)-((a)/(b*R*t(i)))*log(1+d_max*b)+z_g-1-log(z_g);
    dif_meu = (meu_g - meu_l);
    val = dif_meu;
    flag = 1;
    count = 0;
    p = 1;
    
    
    while(1)
        %count = count + 1;
        
        Pg = Pg - 0.001;
        c0 = -a*b;
        c1 = -(b^2 + R*t(i)*b - a);
        c2 = -R*t(i);
        c3 = Pg;
        x = roots([c3 c2 c1 c0]);
        y = x(imag(x) == 0);
        l = length(y);
        
        v_min = min(y);
        d_min = 1/v_min; % density corresponding to lowest volume root....
        v_max = max(y);
        d_max = 1/v_max; % density correaponding to highest volume root...
        
        z_l = (1/(1-d_min*b))-((a*d_min)/(R*t(i)*(1+d_min*b)));     
        z_g = (1/(1-d_max*b))-((a*d_max)/(R*t(i)*(1+d_max*b)));
        meu_l = -log(1-d_min*b)-((a)/(b*R*t(i)))*log(1+d_min*b)+z_l-1-log(z_l);  
        meu_g = -log(1-d_max*b)-((a)/(b*R*t(i)))*log(1+d_max*b)+z_g-1-log(z_g);
        dif_meu = (meu_g - meu_l);
        
        if val > dif_meu
            flag = 1;  %if dif_meu converges on decreasing Pg, then we will keep decreasing
        else
            flag = -1;
        end
        
        val = dif_meu;
        
        if (dif_meu < tol)
            break;
        end
       
    end
    
    Ps(i) = Pg % since the difference in chemical potentials was calculated before pressure update, the value for which it converged will be Pg-0.2. 
    v_l(i) = v_min
    v_g(i) = v_max
    
      
    %plotting.....
    hold on
    
    plot(v_l(i),Ps(i),'*');
    plot(v_g(i),Ps(i),'*');
    fplot (@(x) Ps(i),[v_l(i) v_g(i)], 'color', 'blue');
    
    hold on
    
end

p = plot(v_l, Ps, 'color','red');
p = plot(v_g, Ps, 'color','red');
title 'Saturation curve';
xlabel 'Vm, molar volume (in L/mol)';
ylabel 'P, pressure (in bar)';
grid 'on';


