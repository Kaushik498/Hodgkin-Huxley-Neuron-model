%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       It is helpful to use a reduced model with only two state variables to study the effects of noninactivation. 
%       Because in classical HH model, we have 4 state variables (V, n, m, h). Simulating this model in any real life
%       problem specifically while observing a network, is not very efficient. So this model is often reduced to a two
%       variable system with n and V as the only state variables by appropriately setting m and h. It should show that
%       the reduced model has the same general behavior as the full model, with minor differences. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
bh = @(V) 1./(exp(-(V+30)/10)+1) 


%% reduced model : 

%% 
for I = [0 1]
    
    % classical hh model :
    func = @(t,p) [ (1/C)*(I-(gK*p(2).^4.*(p(1) - vK)) - (gNa*p(3).^3.*p(4).*(p(1)-vNa))-(gL.*(p(1)-vL))); ( (an(p(1)).*(1-p(2))) - (bn(p(1)).*p(2)) ) ; ( (am(p(1)).*(1-p(3))) - (bm(p(1)).*p(3)) ) ; ( (ah(p(1)).*(1-p(4))) - (bh(p(1)).*p(4)) ) ] ;
    [T, Parameters] = ode15s(func, tspan, [-60 0.317 0.0529 0.596]);

    % reduced hh model :
    func3 = @(t,y) [ (1/C)*(I-gK*y(2).^4.*(y(1) - vK) - gNa*(am(y(1))./(am(y(1))+bm(y(1)))).^3.*(ah(y(1))/(ah(y(1))+bh(y(1)))).*(y(1)-vNa)-gL.*(y(1)-vL) ); an(y(1)).*(1-y(2))-bn(y(1)).*y(2)];
    [Tn, parm] = ode15s(func3 , tspan, [-60 0.317] ) ;
    
    % plot :
    figure
    plot(Tn, parm(:,1) , T, Parameters(:,1))
    legend ('Reduced HH model' , 'Classical HH model') ;
    title (['membrane potential vs time for Iext = ' num2str(I) '\muA/cm^2'])

end
