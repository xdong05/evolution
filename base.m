clear all;
rng(100580);

rain = 60;
C_P   = 0.05;
C_chi = 10*10^(-3);
yr1   = 10000; 
raindropto = 5;

n = 300;
nYrs = yr1+ceil((rain-raindropto)/C_P);  

dx = 2;          %grid length(m)
dt = 0.04/365;   %time step in years
nt = nYrs/dt;    %total time steps

W = zeros(n,nYrs); %mm
P = zeros(n,nYrs); %g/m2
O = zeros(n,nYrs); %mm
CHI = zeros(n,nYrs); 
rain_record = zeros(1,nYrs);

Gamma = 20;  %  water uptake rate -- m^2/(kg y);
c = 0.032/20;
f = 0.01;    %  infiltration contrast -- no unit
A = 40.0;    %  maximal value of infiltration rate (I) -- 1/y
Q = 0.06;    %  reference biomass at I ~ A/2 for f << 1 -- kg/m^2
L_0 = 4.0;   %  evaporation rate in bare soil -- 1/y
K_min = 0.1; %  min reference biomass for 50% growth attenuation -- kg/m^2
K_max = 0.6; %  max reference biomass for 50% growth attenuation -- kg/m^2
M_min = 0.5; %  min mortality rate -- 1/y (Default: 0.5)
M_max = 0.9; %  max mortality rate -- 1/y (Default: 0.9)
Dp = 1;      %  biomass dispersal rate -- m^2/y
Dw = 100;    %  below ground water diffusion coef -- m^2/y
Do = 10000;  %  above ground water diffusion coef -- m^2/y

% initial condition
a    = 0.1; b = 0.9;  % 0.1-0.15; 0.55-0.6
CHI(:,1) = a + (b - a).* rand(n,1);
iCov = max(rand, 0.004);
bioDensity = 0.2*rand;
bioDensity = max(bioDensity, 0.0001);
P(randsample(n,round(n*iCov)),1) = rand(round(n*iCov),1).*bioDensity;
W(:,1) = 5;
O(:,1) = 5;

w0 = W(:,1);
p0 = P(:,1);
o0 = O(:,1);
chi0 = CHI(:,1);
rain_record(1) = rain;

eta_M = M_max - M_min;
eta_K = K_max - K_min;


for t = 2:nt

    K = K_max + chi0.* (K_min - K_max);
    M = M_max + chi0.* (M_min - M_max);
    Y = 0.9;
    Lambda = (c*Gamma).* ( 1 - p0./K);
    I = A.*( Y.* p0 + f * Q )./( Y.* p0 + Q );
    
    
    p1 = p0 + ((Lambda.* w0).* p0 - M.*p0).* dt + (Dp.* my_laplacian(p0, n, dx)).* dt;
    w1 = w0 + (I.* o0 - L_0.* w0 - (Gamma.* w0).* p0).* dt + (Dw.* my_laplacian(w0, n, dx)).* dt;
    o1 = o0 + (rain - I.* o0).* dt + (Do.* my_laplacian(o0, n, dx)).* dt;
    drdchi = -(((eta_K*c*Gamma).*p0).*w0)./(K_max - eta_K.*chi0).^2 + eta_M;
    chi1   = chi0 + (C_chi.*drdchi).*dt + (Dp.* my_laplacian(chi0, n, dx)).* dt +...
        ((2*Dp).*my_advection(p0,chi0,dx)).*dt;

    chi1(chi1 < 0) = 0;   
    chi1(chi1 > 1) = 1;  


    % Integrate into a time/space matrix
    if rem((t-1)*dt,1) == 0
        W(:,(t-1)*dt+1) = w1;
        P(:,(t-1)*dt+1) = p1;
        O(:,(t-1)*dt+1) = o1;
        CHI(:,(t-1)*dt+1) = chi1;

        if ((t-1)*dt+1) > (yr1-1)
            rain = rain - C_P;
        end
        
        rain_record((t-1)*dt+1) = rain;
    end

    if rem((t-1)*dt,10) == 0
        round(t/nt.*100,1) % display percent finish 
    end

    % watch out for numerical instability
    if min(o1(:))<0 | isnan(min(o1(:))) | max(chi1(:))>1 | min(chi1(:))< 0 
        break
    end

    % prepare for the next loop
    w0 = w1;
    p0 = p1;
    o0 = o1;
    chi0 = chi1;

end

