clear 
close all
clc
tic

% Soil Layer Parameters
H=50;
lambda=[(1E-5)/H 1/H 2/H 3/H 4/H 5/H 6/H 7/H 8/H 9/H 10/H];
alpha= 0.1 * ones(1, length(lambda));

Finhomo_Linear=zeros(length(lambda),1 );
Finhomo_Parabolic=zeros(length(lambda), 1);
Finhomo_Sinsusoidal=zeros(length(lambda), 1);


for o=1:length(lambda)


Finhomo_Linear(o) = (sqrt(6) * H * sqrt((exp(-2*H*lambda(o)) * (4*exp(H*lambda(o)) - alpha(o)) * alpha(o) + (-4 + alpha(o)) * alpha(o) + 2 * H * lambda(o)) / (H^3 * lambda(o)))) / pi;

Finhomo_Parabolic(o) = (sqrt(5/2) * H * sqrt((-48 * alpha(o)  + 3 * alpha(o) ^2 + 4 * H^3 * lambda(o) ^3 - 3 * exp(-2 * H * lambda(o) ) * alpha(o) ^2 * (1 + 2 * H * lambda(o)  * (1 + H * lambda(o) )) + 24 * exp(-H * lambda(o) ) * alpha(o)  * (2 + H * lambda(o)  * (2 + H * lambda(o) ))) / (H^5 * lambda(o) ^3))) / pi;

Finhomo_Sinsusoidal(o) = (H * sqrt((exp(-2*H*lambda(o)) * (4 * exp(H*lambda(o)) * alpha(o) * (pi^4 + 6 * H^2 * pi^2 * lambda(o)^2 + 8 * H^4 * lambda(o)^4) - alpha(o)^2 * (pi^4 + 9 * H^2 * pi^2 * lambda(o)^2 + 8 * H^4 * lambda(o)^4) + exp(2 * H * lambda(o)) * (8 * H^5 * lambda(o)^5 + pi^4 * ((-4 + alpha(o)) * alpha(o) + 2 * H * lambda(o)) + H^2 * pi^2 * lambda(o)^2 * ((-16 + alpha(o)) * alpha(o) + 10 * H * lambda(o))))) / (H^3 * lambda(o) * (pi^4 + 5 * H^2 * pi^2 * lambda(o)^2 + 4 * H^4 * lambda(o)^4)))) / sqrt(2);



end

RayleighHZ=[Finhomo_Linear, Finhomo_Parabolic, Finhomo_Sinsusoidal];


toc