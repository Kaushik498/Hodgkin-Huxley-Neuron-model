%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       One cause of a paralytic muscle disease called myotonia is loss of Na+ channel inactivation due to a genetic defect. 
%       The action potential will be shown to be very sensitive to small changes in Na+ channel inactivation. The above situation 
%       is modeled by assuming that a fraction fni of Na+ channels do not inactivate, so that the sodium current in this HH model 
%       is given by:
%       I_Na = gNa(1-fni)*(m^3)*h*(V - vNa) + gNa*fni*(m^3)*(V - vNa) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%% Including Sodium channel deactivation factor fni in hh model :

for fni = [0 0.1 0.17 0.2]
    
    func2 = @(t,p) [ (1/C)*(I - (gK*p(2).^4.*(p(1) - vK)) - (gNa*(1 - fni)*p(3).^3.*p(4).*(p(1)-vNa)) - (gNa*fni*p(3).^3.*(p(1)-vNa)) - (gL.*(p(1)-vL))) ; ( (an(p(1)).*(1-p(2))) - (bn(p(1)).*p(2)) ) ; ( (am(p(1)).*(1-p(3))) - (bm(p(1)).*p(3)) ); ( (ah(p(1)).*(1-p(4))) - (bh(p(1)).*p(4)) ) ] ;
    
    [T1, Parm1] = ode15s(func2, [0 10], [-53 0.317 0.0529 0.596] );
    
    [T2, Parm2] = ode15s(func2, [0 10], [-52 0.317 0.0529 0.596] );
    
    [T3, Parm3] = ode15s(func2, [0 10], [-51 0.317 0.0529 0.596] );
    
    [T4, Parm4] = ode15s(func2, [0 10], [-50 0.317 0.0529 0.596] );
    
    % membrane potential vs time :
    figure
    plot(T1,Parm1(:,1) , T2,Parm2(:,1) , T3,Parm3(:,1) , T4,Parm4(:,1) );
    title(['membrane potential for fni = ' num2str(fni)] ) ; 
    
    % Phase plane plot of the model :
    figure
    [T, P] = ode15s(func4, [0 100], [-53 0.317] );
    plot(P(:,2) , P(:,1) );
    title( 'phase plane plot for HH model including fni' ) ; 
    hold on ;

end
