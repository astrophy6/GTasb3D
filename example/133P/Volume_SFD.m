function [Volume_smaller] = Volume_SFD( D_max_lift )
% integrate the accumulative particle volume with size smaller than D_max

D_min = 100e-6;
alpha = -3.3;

Volume_smaller = pi/6.0 * ( D_min^(3+alpha) - D_max_lift^(3+alpha) ) * alpha/(3+alpha);

end