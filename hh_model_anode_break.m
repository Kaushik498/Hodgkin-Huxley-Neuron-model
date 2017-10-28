%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       In the membrane the phenomenon of anode break excitation produces an action potential at the end of a
%       hyperpolarizing stimulus pulse. This model demonstrates anode break excitation (with a hyperpolarizing
%       current of -3 μA/cm2 for 20 ms).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%% Anode break :

for i=1:3 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   this loop is divided into three parts. first part when 0 external current is applied for 100 mSec , then negative 
    %   current is passed for 20 mSec , in the third part it is removed and again no current is applied for another 100 mSec. 
    %   we will observe the anode break at the begining of the third part.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % No external current :
    if i == 1 
        I = 0 ;
        func = @(t,p) [ (1/C)*(I-(gK*p(2).^4.*(p(1) - vK)) - (gNa*p(3).^3.*p(4).*(p(1)-vNa))-(gL.*(p(1)-vL))); ( (an(p(1)).*(1-p(2))) - (bn(p(1)).*p(2)) ) ; ( (am(p(1)).*(1-p(3))) - (bm(p(1)).*p(3)) ) ; ( (ah(p(1)).*(1-p(4))) - (bh(p(1)).*p(4)) ) ] ;
        [T1, P1] = ode15s(func, tspan, [-60 0.317 0.0529 0.596]); 
    
    % hyperpolarizing current :
    elseif i == 2
        I = -3 ;
        
        func = @(t,p) [ (1/C)*(I-(gK*p(2).^4.*(p(1) - vK)) - (gNa*p(3).^3.*p(4).*(p(1)-vNa))-(gL.*(p(1)-vL))); ( (an(p(1)).*(1-p(2))) - (bn(p(1)).*p(2)) ) ; ( (am(p(1)).*(1-p(3))) - (bm(p(1)).*p(3)) ) ; ( (ah(p(1)).*(1-p(4))) - (bh(p(1)).*p(4)) ) ] ;
        [T2, P2] = ode15s(func, [T1(length(T1)) 20+T1(length(T1))], [P1(length(P1) , 1) P1(length(P1) , 2) P1(length(P1) , 3) P1(length(P1) , 4)]);
    
    % hyperpolarizing current removed :
    else i == 3
        I = 0 ;
        
        func = @(t,p) [ (1/C)*(I-(gK*p(2).^4.*(p(1) - vK)) - (gNa*p(3).^3.*p(4).*(p(1)-vNa))-(gL.*(p(1)-vL))); ( (an(p(1)).*(1-p(2))) - (bn(p(1)).*p(2)) ) ; ( (am(p(1)).*(1-p(3))) - (bm(p(1)).*p(3)) ) ; ( (ah(p(1)).*(1-p(4))) - (bh(p(1)).*p(4)) ) ] ;
        [T3, P3] = ode15s(func, [T2(length(T2)) 100+T2(length(T2))], [P2(length(P2) , 1) P2(length(P2) , 2) P2(length(P2) , 3) P2(length(P2) , 4)]);
        
        
        figure
        plot ( [T1 ; T2 ; T3] , [P1(:,1) ; P2(:,1) ; P3(:,1)] ) ;
        title('Anode Break')
    end
end

