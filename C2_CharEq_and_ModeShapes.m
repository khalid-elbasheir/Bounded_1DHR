% This code solve the characteristic equation using initial guesses and
% plot the first three modal shapes for each combination of inhomogeneity
% parameters
clear 
close all
clc
tic

H=50;
lambda=[3/H 1/H 3/H 8/H 3/H 1/H];
alpha=[0.9 0.9 0.8 0.5 0.5 1E-5]; % an alpha=0 generates a singularity in the hypergeometric argument

%====== In a decreasing stiffness with depth profile (alpha<0) uncomment
%the commneted lines below here and at the initial guessess =======
%lambda = [8/H 5/H 3/H 10/H 1/H];
%alpha = [-0.5 -0.9 -1 -3 1E-5];


z=[];
for k = 0:0.2:H
        z = [z k];
end

% Section for Normal Profile

mode_shapes = cell(length(lambda), 1);
C_m_Values = zeros(1, 3);
C_m_Values_lambda=cell(length(lambda), 1);
mode_shapes_for_lambda = cell(3, 1);




for o=1:length(lambda)
lam = lambda(o);
lamH = lam*H;
a = alpha(o);

    if o==1
        Iguess= [0.4317 0.9813 1.5455];
        %Iguess = [0.1977 0.6079 1.0237]; %Decreasing stiffness Case
    elseif o==2   
        Iguess=[0.7687 1.7208 2.88888];
        %Iguess = [0.3276 1.0515 1.7792]; %Decreasing Case
    elseif o==3
        Iguess=[0.4469 1.0956 1.7582];
        %Iguess = [0.5942	1.9602	3.3079]; %Decreasing Case
    elseif o==4
        Iguess=[0.1944 0.5632 0.9205];
        %Iguess = [0.1600 0.55555555 0.8805]; %Decreasing Case
    elseif o==5
        Iguess=[0.4804 1.3160 2.1607];
        %For Decreasing Case Comment the fifth initial guess above
    elseif o==6 % For Decreasing case change this to 0==5
        Iguess=[3 5 7];
    end

    for m = 1:3 
     

   char_eq = @(C_m) real((lam*exp(lamH))*((((exp(lamH) - a)^((-1 + sqrt(1 - 4*C_m^2))/2))*((1 + sqrt(1 - 4*C_m^2))/2)*hypergeom([(1 - 2*1i*C_m + sqrt(1 - 4*C_m^2))/2, (1 + 2*1i*C_m + sqrt(1 - 4*C_m^2))/2], 1 + sqrt(1 - 4*C_m^2), -((exp(lamH) - a)/a)))...
    -((exp(lamH) - a)^((-1 -sqrt(1 - 4*C_m^2))/2))*((1 - sqrt(1 - 4*C_m^2))/2)*((1 - a)^sqrt(1 - 4*C_m^2))*(hypergeom([(1 - 2*1i*C_m + sqrt(1 - 4*C_m^2))/2,  (1 + 2*1i*C_m + sqrt(1 - 4*C_m^2))/2], 1+ sqrt(1 - 4*C_m^2), -((1 - a)/a))*hypergeom([(1 - 2*1i*C_m - sqrt(1 - 4*C_m^2))/2, (1 + 2*1i*C_m - sqrt(1 - 4*C_m^2))/2], 1 - sqrt(1 - 4*C_m^2), -((exp(lamH) - a)/a))/hypergeom([(1 - 2*1i*C_m - sqrt(1 - 4*C_m^2))/2, (1 + 2*1i*C_m - sqrt(1 - 4*C_m^2))/2], 1 - sqrt(1 - 4*C_m^2), -((1 - a)/a)))) ...
      + (((exp(lamH) - a)^((1 + sqrt(1 - 4*C_m^2))/2))*((- lam*exp(lamH))/(2*a))* hypergeom([(3 - 2*1i*C_m + sqrt(1 - 4*C_m^2))/2, (3 + 2*1i*C_m + sqrt(1 - 4*C_m^2))/2], 2 + sqrt(1 - 4*C_m^2), -((exp(lamH) - a)/a))-((exp(lamH) - a)^((1 - sqrt(1 - 4*C_m^2))/2))*((1 - a)^sqrt(1 - 4*C_m^2))*(hypergeom([(1 - 2*1i*C_m + sqrt(1 - 4*C_m^2))/2,  (1 + 2*1i*C_m + sqrt(1 - 4*C_m^2))/2], 1+ sqrt(1 - 4*C_m^2), -((1 - a)/a))*((- lam*exp(lamH))/(2*a))* hypergeom([(3 - 2*1i*C_m - sqrt(1 - 4*C_m^2))/2, (3 + 2*1i*C_m - sqrt(1 - 4*C_m^2))/2], 2 - sqrt(1 - 4*C_m^2), -((exp(lamH) - a)/a))/hypergeom([(1 - 2*1i*C_m - sqrt(1 - 4*C_m^2))/2, (1 + 2*1i*C_m - sqrt(1 - 4*C_m^2))/2], 1 - sqrt(1 - 4*C_m^2), -((1 - a)/a)))));
   
    options = optimset('Display', 'iter','TolFun', 1e-6, 'TolX', 1e-6);
    C_m(m) = fzero(char_eq, Iguess(m), options); % Start with an initial guess of 1
     
    C_m_Values(m) =C_m(m);  
       
    parfor i = 1:length(z)
        ps = C_m(m);
        ct = 2*1i*ps;
        sqf = sqrt(1 - 4 * ps^2);
        s1 = 0.5 * (1 + sqf);
        s2 = 0.5 * (1 - sqf);
        a1 = 0.5 * (1 - ct + sqf);
        b1 = 0.5 * (1 + ct + sqf);
        c1 = 1 + sqf;
        a2 = 0.5 * (1 - ct - sqf);
        b2 = 0.5 * (1 + ct - sqf);
        c2 = 1 - sqf;
        Arg_Z(i) = -(exp(z(i) * lam) - a) / a;
        Arg_0 = -(1 - a) / a;
        % Define hypergeometric functions F_1 and F_2 at z and 0
        F_1_z(i) =  hypergeom([a1,b1],c1,Arg_Z(i));
        F_2_z(i)=  hypergeom([a2,b2],c2,Arg_Z(i));
        F_1_0 =  hypergeom([a1,b1],c1,Arg_0);
        F_2_0 = hypergeom([a2,b2],c2,Arg_0);
        % Define hypergeometric functions phi_1 and phi_2 at z and 0
        phi_1_z(i) =  (-lam * exp(z(i) * lam) ) / (2 * a) * hypergeom([a1+1,b1+1],c1+1,Arg_Z(i));
        phi_2_z(i) = (-lam * exp(z(i) * lam) ) / (2 * a) * hypergeom([a2+1,b2+1],c2+1,Arg_Z(i));
        phi_1_0 = -lam / (2 * a) * hypergeom([a1+1,b1+1],c1+1,Arg_0);
        phi_2_0 =  -lam / (2 * a) * hypergeom([a2+1,b2+1],c2+1,Arg_0);


        % Define E and M
        E(i) =  lam * exp(z(i) * lam) * (s1)* (exp(z(i) * lam) - a)^(-s2) * F_1_z(i) + (exp(z(i) * lam) - a)^(s1) * phi_1_z(i);

        M(i) =  lam * exp(z(i) * lam) * (s2)* (exp(z(i) * lam) - a)^(-s1) * F_2_z(i) + (exp(z(i) * lam) - a)^(s2) * phi_2_z(i);
        % Define E and M at z=0
        E_0 = lam * 1 * (s1)* (1 - a)^(-s2) * F_1_0 + (1 - a)^(s1) * phi_1_0;
        M_0 =  lam * 1 * (s2)* (1 - a)^(-s1) * F_2_0 + (1 - a)^(s2) * phi_2_0;
        Term_Inside = (1 - a)^(sqf) * (F_1_0 /F_2_0 ) ;
        % Define the mode shape Phi
        Phi(i) =  ( E(i) - (Term_Inside*M(i)) ) / (E_0 - (Term_Inside*M_0) );

    end
 mode_shapes_for_lambda{m} = (Phi);
    end
   
    mode_shapes{o} = mode_shapes_for_lambda;
    C_m_Values_lambda{o}=C_m_Values;
end

figure(1);
for m = 1:3
    subplot(1, 3, m); % Adjust the subplot position
    hold on;
    for o = 1:length(lambda)
            plot(mode_shapes{o}{m}, z/H);
    end
    xlabel(['\Phi_' num2str(m)]);
    ylabel('z/H');
    xlim([-1.5 1.5]);
    %xlim([-2 2]); %Uncomment this in case of decreasing stiffness
   
   legendEntries = arrayfun(@(o) sprintf('\\lambdaH = %.2f, \\alpha = %.2f', ...
    lambda(o)*H, alpha(o)), 1:length(lambda), 'UniformOutput', false);
   set(gca, 'YDir', 'reverse');
   legend(legendEntries, 'Location', 'Best','FontSize',10);



end
toc;