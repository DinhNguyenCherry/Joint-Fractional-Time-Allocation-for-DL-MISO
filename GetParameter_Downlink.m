function [nu, xi, lambda] = GetParameter_Downlink(vartheta, mu)


nu = 2/mu*log(1+1/vartheta) + 1/((vartheta+1)*mu);

xi = -1/(vartheta*(vartheta+1)*mu);

lambda = -1/(mu^2)*log(1+1/vartheta);



% varphi = -log(1 - ( (real(h'*w)^2) ) / ( phi^2 + real(h'*w)^2 ) ) - real(h'*w)^2/(phi^2);
% 
% chi = 2*real(h'*w)/(phi^2);
% 
% varpi = real( (h'*w)^2 ) / ( phi^2 * ( phi^2 + real(h'*w)^2 ) );

end