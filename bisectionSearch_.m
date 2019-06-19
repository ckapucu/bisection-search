clear; clc;
NumPoints = 15000;

% Gi and Tc are for trying different operational conditions.
% If Gi=1000 and Tc=25 it will give same outputs with STC
Gi = 1000;  
Tc = 25;   

%constant data 
k=1.3806503e-23;                            % Boltzman constant
q=1.60217646e-19;                           % electron charge
Tn = 25+273.15;                             % temperature of the cell at STC: 298.15 K (K)
Gn = 1000;                                  % surface irradiance of the cell at STC:1000W/m2 (W/m2)

%parameters from manufacturer’s datasheet (KC200GT) Polycrystalline
Vmn = 26.3;                               % maximum power voltage at STC (V)
Imn = 7.61;                               % maximum power current at STC (A)
Vocn = 32.9;                              % open circuit voltage at STC (V)
Iscn = 8.21;                              % short circuit current at STC (A)
Kv = -0.123;                              % open circuit voltage temperature coefficient (V/*C)
Ki = 0.00318;                             % short circuit current temperature coefficient (A/*C)
Ns = 54;                                  % number of cells in series in one module
Np = 1;                                   % number of series in parallel in one module
Pmn = Vmn*Imn;                            % maximum power at STC (W)

%parameters from manufacturer’s datasheet (SP70) Monocrystalline
% Vmn = 16.5;                               % maximum power voltage at STC (V)
% Imn = 4.25;                               % maximum power current at STC (A)
% Vocn = 21.4;                              % open circuit voltage at STC (V)
% Iscn = 4.7;                               % short circuit current at STC (A)
% Kv = -0.076;                              % open circuit voltage temperature coefficient (V/*C)
% Ki = 0.002;                               % short circuit current temperature coefficient (A/*C)
% Ns = 36;                                  % number of cells in series in one module
% Np = 1;                                   % number of series in parallel in one module
% Pmn = Vmn*Imn;                            % maximum power at STC (W)

%parameters from manufacturer’s datasheet (ST40) Thin-film
% Vmn = 16.6;                               % maximum power voltage at STC (V)
% Imn = 2.41;                               % maximum power current at STC (A)
% Vocn = 23.3;                              % open circuit voltage at STC (V)
% Iscn = 2.68;                              % short circuit current at STC (A)
% Kv = -0.1;                                % open circuit voltage temperature coefficient (V/*C)
% Ki = 0.00035;                             % short circuit current temperature coefficient (A/*C)
% Ns = 36;                                  % number of cells in series in one module
% Np = 1;                                   % number of series in parallel in one module
% Pmn = Vmn*Imn;                            % maximum power at STC (W)

%parameters from manufacturer’s datasheet (MX60) Polycrystalline
% Vmn = 17.1;                               % maximum power voltage at STC (V)
% Imn = 3.5;                                % maximum power current at STC (A)
% Vocn = 21.1;                              % open circuit voltage at STC (V)
% Iscn = 3.8;                               % short circuit current at STC (A)
% Kv = -0.08;                               % open circuit voltage temperature coefficient (V/*C)
% Ki = 0.003;                               % short circuit current temperature coefficient (A/*C)
% Ns = 36;                                  % number of cells in series in one module
% Np = 1;                                   % number of series in parallel in one module
% Pmn = Vmn*Imn;                            % maximum power at STC (W)

%diode parameters
Vtn=k*Ns*Tn/q;                              % thermal voltage of the diode at STC
an=1.3;                                     % 1.3 has been suggested for monocrystalline and polycrystalline PV modules [11].
Vtn_an=Vtn*an;

%parameters under operation conditions
To = Tc + 273.15;                           % temperature of the cell (K)
Go = Gi;                                    % surface irradiance of the cell (W/m2)
%running parameters
nimax = 1000; %maximum iterations
tol = 0.0001; %tolerance 
dP_dV =zeros(1, 50);
% Start the timer
tic

%reverse saturation current Eq. (4)
Ion=Iscn/(exp(Vocn/Vtn_an)-1);
fprintf('Ion= %g ',Ion);
%Rsnr is the maximum value of Rs.n Eq. (9)

Rsnr=((Vtn_an*log((Iscn-Imn)/Ion+1)-Vmn)/Imn);
fprintf('Rsnr= %f ',Rsnr);

Rsn=0;
fprintf('Rsn = %f ',Rsn);
Rpn=(Vmn+Imn*Rsn-Iscn*Rsn)/(Iscn-Ion-Imn-Ion* exp((Vmn+Imn*Rsn)/Vtn_an));%Eq. (8) 
fprintf('Rpn = %f ',Rpn);

dPm_dVm_1=Imn+Vmn*(-Ion/Vtn_an*exp((Vmn+Imn*Rsn)/Vtn_an)-1/ Rpn)/(1+Ion*Rsn/Vtn_an*exp((Vmn+Imn*Rsn)/Vtn_an)+Rsn/Rpn);% Eq. (3)
fprintf('dPm_dVm_1 = %f\n',dPm_dVm_1);


ni = 1; %iteration counter
aa = 0;
bb = Rsnr;


while(ni<=nimax)
    fprintf('%2d. step  ',ni);
    Rsn = (aa + bb)/2; %new midpoint
    fprintf('Rsn = %f ',Rsn);
    Rpn=(Vmn+Imn*Rsn-Iscn*Rsn)/(Iscn-Ion-Imn-Ion*exp((Vmn+Imn*Rsn)/Vtn_an));%Eq. (8)
    fprintf('Rpn = %f ',Rpn);

    dP_dV(ni)=Imn+Vmn*(-Ion/Vtn_an*exp((Vmn+Imn*Rsn)/Vtn_an)- 1/Rpn)/(1+Ion*Rsn/Vtn_an*exp((Vmn+Imn*Rsn)/Vtn_an)+Rsn/Rpn); %Eq. (2)
    
    fprintf('dPm_dVm_2 = %f ',dP_dV(ni));
    fprintf('(bb-aa)/2 = %f ',(bb-aa)/2);
    fprintf('tol = %f\n',tol);

    if(abs(dP_dV(ni))<=tol && ((bb-aa)/2)<tol)
        fprintf('Rsn = %f ',Rsn);
        Rpn=(Vmn+Imn*Rsn-Iscn*Rsn)/(Iscn-Ion-Imn-Ion*exp((Vmn+Imn*Rsn)/Vtn_an)); %Eq. (8)
        fprintf('Rpn = %f ',Rpn);
        ILn=(Rsn+Rpn)/Rpn*Iscn+Ion*exp((Iscn*Rsn)/Vtn_an)-Ion; %Eq. (6)
        dP_dV(ni)=Imn+Vmn*(-Ion/Vtn_an*exp((Vmn+Imn*Rsn)/Vtn_an)- 1/Rpn)/(1+Ion*Rsn/Vtn_an*exp((Vmn+Imn*Rsn)/Vtn_an)+Rsn/Rpn);
        fprintf('f(c) = %f\n',dP_dV(ni));
        break;
    end

    ni = ni+1;
    if(dPm_dVm_1*dP_dV(ni-1)>0)
        aa = Rsn;
    else
        bb = Rsn;
    end
end


% Stop the timer 
toc
fprintf('\nParameters at STC with bisection search \n');

fprintf('  Step = %d\n',ni);
fprintf('    Gn = %f\n',Gn);
fprintf('    Tn = %f\n',Tn-273.15);
fprintf('   Rsn = %f\n',Rsn);
fprintf('   Rpn = %f\n',Rpn);
fprintf('    an = %f\n',an);
fprintf('   ILn = %f\n',ILn);
fprintf('   Ion = %g\n',Ion);
fprintf('   Imn = %g\n',Imn);
fprintf('   Vmn = %g\n',Vmn);
fprintf('   Pmn = %g\n',Pmn);
fprintf('\nCalculated parameters under operation condiditons (u.o.c)\n');

%calculated parameters under operation conditions
fprintf('     G = %f\n',Go);
fprintf('     T = %f\n',To-273.15);
Rs_ = Rsn;                      %Eq. (13)
fprintf('    Rs = %f\n',Rs_);  
Rp_ = (Gn/Go)*Rpn;              %Eq. (14)
fprintf('    Rp = %f\n',Rp_);
dT_ = To-Tn;                     

a = an; 

Vt = k*Ns*To/q;
Vt_a_ = a*Vt; % modified diode ideality factor 

Voc = Vocn + Kv*dT_ + Vt_a_*log(Go/Gn); %Eq. (15)
Isc = (Iscn + Ki * dT_) * (Go/Gn);      %(Villalva et. al, 2009 Eq.(4) )

Impp = Imn*(Go/Gn);
Vmpp = Vmn+Kv*dT_;
IL_ = (ILn + Ki*dT_)*(Go/Gn);           %Eq. (11)
Pmpp = Impp*Vmpp;
%NONEED
Io_ = (Iscn + Ki*dT_)/(exp((Vocn + Kv*dT_)/Vt_a_)-1);   %Eq. (12)

fprintf('    IL = %f\n',IL_);
fprintf('    Io = %g\n',Io_);
fprintf('   Voc = %f\n',Voc);
fprintf('   Isc = %f   (Villalva et. al, 2009 Eq.(4) )\n',Isc);
fprintf('     a = %f\n',a);



i = 0;
idx=1;
Iout = zeros(1, length(0:Voc/NumPoints:Voc));
for V=0:Voc/NumPoints:Voc
    Iout(idx)= IL_ - Io_.*(exp((V+(i.*Rs_))./Vt_a_)-1)-((V+(i.*Rs_))/Rp_);
    i = Iout(idx);                     %Update Current
    idx=idx+1;
end

V=0:Voc/NumPoints:Voc;
P = Iout.*V;
%find the maximum power value like MPPT devices do.
maxx = find(P(1,:)==max(P(1,:)));
last = find(P(P(1,:)>0), 1, 'last' );

fprintf('\nModule output\n');
fprintf('     I = %f\n',Iout(maxx));
fprintf('     V = %f\n',V(maxx));
fprintf('     P = %f\n',max(P));

% I-V curve
figure(1) 
grid on
hold on 
title('');
xlabel('Voltage (V)');
ylabel('Current (A)');
xlim([0 max(V)*1.1]);
ylim([0 max(Iout)*1.1]);
plot(V,Iout,'LineWidth',2,'Color','k') %
%markers on the curve
plot([0 V(maxx) V(last) ],[max(Iout) Iout(maxx) 0 ],'o','LineWidth',2,'MarkerSize',4,'Color','k') 
title('I-V curve', 'Color', 'k');  

% P-V curve
figure(2) 
grid on
hold on 
title('');
xlabel('Voltage (V)');
ylabel('Power (W)');
xlim([0 Voc*1.1]);
ylim([0 max(P(1,:))*1.1]);  
plot(V,P,'LineWidth',2,'Color','k') %
%markers on the curve
plot([0 V(maxx) V(last) ],[0 max(P(1,:)) 0 ],'o','LineWidth',2,'MarkerSize',4,'Color','k')
title('P-V curve', 'Color', 'k'); 