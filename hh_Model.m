%% declaring variables :

global gK gNa gL vK vNa vL phi C Q ;


%% assigning values and equations :

gK = 36 ;               %% unit: mS/cm^2
vK = -72 ;              %% unit: mV
gNa = 120 ;             %% unit: mS/cm^2
vNa = 55 ;              %% unit: mV
gL = 0.3 ;              %% unit: mS/cm^2
vL = -49.4 ;            %% unit: mV
C = 1 ;                 %% unit: muF/cm^2
phi = 1 ; Q = 3 ;       %% Temperature coefficient-unitless ; phi = Q^((T-6.3)/10) ; here T is taken as 6.3 . You can vary temperature as well, if you want.
tspan = [0 100] ;


%% Defining opening rate(alpha) and closing rate(beta) of K channel(n), Na channel(m,h) :

an = @(V) (-0.01*(V+50))./(exp(-(V+50)/10)-1) ;
bn = @(V) 0.125*exp(-(V+60)/80) ;

am = @(V) (-0.1*(V+35))./(exp(-(V+35)/10)-1) ;
bm = @(V) 4*exp(-(V+60)/18) ;

ah = @(V) 0.07*exp(-(V+60)/20) ;
bh = @(V) 1./(exp(-(V+30)/10)+1) ;


%% Plots of gating Probabilities and membrane potential with various current injection value :

for I = 0:1:10 ;
    
    % defining hh model equation :
    
    func = @(t,p) [ (1/C)*(I-(gK*p(2).^4.*(p(1) - vK)) - (gNa*p(3).^3.*p(4).*(p(1)-vNa))-(gL.*(p(1)-vL))); ( (an(p(1)).*(1-p(2))) - (bn(p(1)).*p(2)) ) ; ( (am(p(1)).*(1-p(3))) - (bm(p(1)).*p(3)) ) ; ( (ah(p(1)).*(1-p(4))) - (bh(p(1)).*p(4)) ) ] ;
    
    [T, Parameters] = ode15s(func, tspan, [-60 0.317 0.0529 0.596]);
    
    figure
    plot(T,Parameters(:,1));
    xlabel('time');
    ylabel('membrane potential (mV)') ;
    title(['Membrane potential vs time for I = ' num2str(I) '\muA/cm^2']);

    figure
    plot(T,Parameters(:,2), T,Parameters(:,3), T, Parameters(:,4));
    legend('n(t)', 'm(t)', 'h(t)') ;
    xlabel('time') ;
    title(['n, m & h vs time for I = ' num2str(I) '\muA/cm^2']);

end


%% Effect of varying Leakage voltage on membrane potential :

for vL = [-49.2 -49.3 -49.4 -49.45 -49.5]
    
    I = 0 ; % taking external current zero for the time being
    
    % defining hodgkin-huxley model equation :
    func = @(t,p) [ (1/C)*(I-(gK*p(2).^4.*(p(1) - vK)) - (gNa*p(3).^3.*p(4).*(p(1)-vNa))-(gL.*(p(1)-vL))); ( (an(p(1)).*(1-p(2))) - (bn(p(1)).*p(2)) ) ; ( (am(p(1)).*(1-p(3))) - (bm(p(1)).*p(3)) ) ; ( (ah(p(1)).*(1-p(4))) - (bh(p(1)).*p(4)) ) ] ;
    
    [T, Parameters] = ode15s(func, tspan, [-60 0.317 0.0529 0.596]);
    
    if vL == -49.2 
        plot(T,Parameters(:,1) , 'b');
    end
    
    if vL == -49.3 
        plot(T,Parameters(:,1) , 'g');
    end
    
    if vL == -49.4 
        plot(T,Parameters(:,1) , 'y');
    end
    
    if vL == -49.45 
        plot(T,Parameters(:,1) , 'c');
    end
    
    if vL == -49.5 
        plot(T,Parameters(:,1) , 'm');
    end
    
    hold on ;
    
    xlabel('time');
    ylabel('membrane potential (mV)') ;
    legend('Eleak =-49.2 V' , 'Eleak =-49.3 V' , 'Eleak =-49.4 V' , 'Eleak =-49.45 V' , 'Eleak =-49.5 V');

end


%% Evaluating equilibrium points of HH equation : 

for I = 8:12 ;
    
    % assigning sym variables :
    syms V1 n1 m1 h1
    
    %% equilibrium points :
    [ equipt.V1, equipt.n1, equipt.m1, equipt.h1 ] = solve( ((I-gK*n1.^4.*(V1 - vK) - gNa*m1.^3.*h1.*(V1-vNa)-gL.*(V1-vL))/C)==0 , an(V1).*(1-n1)-bn(V1).*n1 ==0 , am(V1).*(1-m1)-bm(V1).*m1 ==0 , ah(V1).*(1-h1)-bh(V1).*h1 ==0 ) ;
    equi_v = equipt.V1 
    equi_n = equipt.n1 
    equi_m = equipt.m1 
    equi_h = equipt.h1 
    
    func1 = @(t,p) [ (1/C)*(I - (gK*p(2).^4.*(p(1) - vK)) - (gNa*p(3).^3.*p(4).*(p(1)-vNa)) - (gL.*(p(1)-vL))) ; ( (an(p(1)).*(1-p(2))) - (bn(p(1)).*p(2)) ) ; ( (am(p(1)).*(1-p(3))) - (bm(p(1)).*p(3)) ); ( (ah(p(1)).*(1-p(4))) - (bh(p(1)).*p(4)) ) ] ;
    [T, Parameters] = ode15s(func, tspan, [double(equipt.V1) double(equipt.m1) double(equipt.n1) double(equipt.h1)] );
    figure
    plot(T,Parameters(:,1));
    title(['for I=' num2str(I) '\muA/cm^2']);
end    
