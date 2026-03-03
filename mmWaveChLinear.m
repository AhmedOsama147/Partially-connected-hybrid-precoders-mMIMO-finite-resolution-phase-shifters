function [H_n,At,Ar] = mmWaveChLinear(NT,NR,Ncl,Nray,DEG)
    gamma = sqrt((NT*NR)/(Ncl*Nray));
    alpha = sqrt(1/(2)) *  (randn(Ncl,Nray) + 1i * randn(Ncl,Nray)) ;
    phiT = (-pi) + ((2*pi)*rand(Ncl,1));
    phiR = (-pi) + ((2*pi)*rand(Ncl,1));
    sig = deg2rad(DEG); % this is the s.d. of the Laplace distribution for the scattering clusters
    At = [];
    Ar = [];
    H = zeros(NR,NT);
    for j=1:Ncl
         % The final azimuth angles from each cluster i    
         phi_AOD = randLaplacian(Nray, 1, phiT(j), sig^2);
         phi_AOA = randLaplacian(Nray, 1, phiR(j), sig^2);    
         for l=1:Nray
             at = exp(-1j*pi*sin(phi_AOD(l))*(0:1:NT-1)).'/sqrt(NT);
             ar = exp(-1j*pi*sin(phi_AOA(l))*(0:1:NR-1)).'/sqrt(NR);
             Alpha = alpha(j,l);
             H = H + Alpha*ar*at'; % channel
             At = [At at];
             Ar = [Ar ar];
         end
     end
     H_n = gamma*H;
end