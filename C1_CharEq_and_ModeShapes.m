clear 
close all
clc
tic

H=50;
lambda=[3/H 1/H 3/H 8/H 3/H 1/H];
alpha=[0.9 0.9 0.8 0.5 0.5 1E-5];



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
    if o==1
        Iguess= [0.1864 0.9630 2.3887];
    elseif o==2   
        Iguess=[0.5909 2.9612 7.6966];
    elseif o==3
        Iguess=[0.1997 1.2003 3.0911];
    elseif o==4
        Iguess=[0.0378	0.3172	0.8474];
    elseif o==5
        Iguess=[0.2308 1.7318 4.6686];
    elseif o==6
        Iguess=[10 30 50];
    end

    for m = 1:3 
             
   char_eq = @(C_m) real((lambda(o)*exp(lambda(o)*H))*((((exp(lambda(o)*H) - alpha(o))^((-1 + sqrt(1 - 4*C_m))/2))*((1 + sqrt(1 - 4*C_m))/2)*hypergeom([(1 - 2*1i*sqrt(C_m) + sqrt(1 - 4*C_m))/2, (1 + 2*1i*sqrt(C_m) + sqrt(1 - 4*C_m))/2], 1 + sqrt(1 - 4*C_m), -((exp(lambda(o)*H) - alpha(o))/alpha(o))))...
    -((exp(lambda(o)*H) - alpha(o))^((-1 -sqrt(1 - 4*C_m))/2))*((1 - sqrt(1 - 4*C_m))/2)*((1 - alpha(o))^sqrt(1 - 4*C_m))*(hypergeom([(1 - 2*1i*sqrt(C_m) + sqrt(1 - 4*C_m))/2,  (1 + 2*1i*sqrt(C_m) + sqrt(1 - 4*C_m))/2], 1+ sqrt(1 - 4*C_m), -((1 - alpha(o))/alpha(o)))*hypergeom([(1 - 2*1i*sqrt(C_m) - sqrt(1 - 4*C_m))/2, (1 + 2*1i*sqrt(C_m) - sqrt(1 - 4*C_m))/2], 1 - sqrt(1 - 4*C_m), -((exp(lambda(o)*H) - alpha(o))/alpha(o)))/hypergeom([(1 - 2*1i*sqrt(C_m) - sqrt(1 - 4*C_m))/2, (1 + 2*1i*sqrt(C_m) - sqrt(1 - 4*C_m))/2], 1 - sqrt(1 - 4*C_m), -((1 - alpha(o))/alpha(o))))) ...
      + (((exp(lambda(o)*H) - alpha(o))^((1 + sqrt(1 - 4*C_m))/2))*((- lambda(o)*exp(lambda(o)*H))/(2*alpha(o)))* hypergeom([(3 - 2*1i*sqrt(C_m) + sqrt(1 - 4*C_m))/2, (3 + 2*1i*sqrt(C_m) + sqrt(1 - 4*C_m))/2], 2 + sqrt(1 - 4*C_m), -((exp(lambda(o)*H) - alpha(o))/alpha(o)))-((exp(lambda(o)*H) - alpha(o))^((1 - sqrt(1 - 4*C_m))/2))*((1 - alpha(o))^sqrt(1 - 4*C_m))*(hypergeom([(1 - 2*1i*sqrt(C_m) + sqrt(1 - 4*C_m))/2,  (1 + 2*1i*sqrt(C_m) + sqrt(1 - 4*C_m))/2], 1+ sqrt(1 - 4*C_m), -((1 - alpha(o))/alpha(o)))*((- lambda(o)*exp(lambda(o)*H))/(2*alpha(o)))* hypergeom([(3 - 2*1i*sqrt(C_m) - sqrt(1 - 4*C_m))/2, (3 + 2*1i*sqrt(C_m) - sqrt(1 - 4*C_m))/2], 2 - sqrt(1 - 4*C_m), -((exp(lambda(o)*H) - alpha(o))/alpha(o)))/hypergeom([(1 - 2*1i*sqrt(C_m) - sqrt(1 - 4*C_m))/2, (1 + 2*1i*sqrt(C_m) - sqrt(1 - 4*C_m))/2], 1 - sqrt(1 - 4*C_m), -((1 - alpha(o))/alpha(o))))));
   
    options = optimset('Display', 'iter','TolFun', 1e-6, 'TolX', 1e-6);
    C_m(m) = fzero(char_eq, Iguess(m), options); % Start with an initial guess of 1
     
    C_m_Values(m) =C_m(m);  
       
 parfor i = 1:length(z)
    % Define hypergeometric functions F_1 and F_2 at z and 0
F_1_z(i) =  hypergeom([0.5 * (1 - 2i * sqrt(C_m(m)) + sqrt(1 - 4 * C_m(m))), ...
                              0.5 * (1 + 2i * sqrt(C_m(m)) + sqrt(1 - 4 * C_m(m)))], ...
                             1 + sqrt(1 - 4 * C_m(m)), ...
                             -(exp(z(i) * lambda(o)) - alpha(o)) / alpha(o));

F_2_z(i)=  hypergeom([0.5 * (1 - 2i * sqrt(C_m(m)) - sqrt(1 - 4 * C_m(m))), ...
                              0.5 * (1 + 2i * sqrt(C_m(m)) - sqrt(1 - 4 * C_m(m)))], ...
                             1 - sqrt(1 - 4 * C_m(m)), ...
                             -(exp(z(i) * lambda(o)) - alpha(o)) / alpha(o));

F_1_0(i) =  hypergeom([0.5 * (1 - 2i * sqrt(C_m(m)) + sqrt(1 - 4 * C_m(m))), ...
                          0.5 * (1 + 2i * sqrt(C_m(m)) + sqrt(1 - 4 * C_m(m)))], ...
                         1 + sqrt(1 - 4 * C_m(m)), ...
                         -(1 - alpha(o)) / alpha(o));

F_2_0(i) = hypergeom([0.5 * (1 - 2i * sqrt(C_m(m)) - sqrt(1 - 4 * C_m(m))), ...
                          0.5 * (1 + 2i * sqrt(C_m(m)) - sqrt(1 - 4 * C_m(m)))], ...
                         1 - sqrt(1 - 4 * C_m(m)), ...
                         -(1 - alpha(o)) / alpha(o));

% Define hypergeometric functions phi_1 and phi_2 at z and 0
phi_1_z(i) =  -lambda(o) * exp(z(i) * lambda(o)) / (2 * alpha(o)) * ...
                      hypergeom([0.5 * (3 - 2i * sqrt(C_m(m)) + sqrt(1 - 4 * C_m(m))), ...
                                 0.5 * (3 + 2i * sqrt(C_m(m)) + sqrt(1 - 4 * C_m(m)))], ...
                                2 + sqrt(1 - 4 * C_m(m)), ...
                                -(exp(z(i) * lambda(o)) - alpha(o)) / alpha(o));

phi_2_z(i) = -lambda(o) * exp(z(i) * lambda(o)) / (2 * alpha(o)) * ...
                      hypergeom([0.5 * (3 - 2i * sqrt(C_m(m)) - sqrt(1 - 4 * C_m(m))), ...
                                 0.5 * (3 + 2i * sqrt(C_m(m)) - sqrt(1 - 4 * C_m(m)))], ...
                                2 - sqrt(1 - 4 * C_m(m)), ...
                                -(exp(z(i) * lambda(o)) - alpha(o)) / alpha(o));

phi_1_0(i) = -lambda(o) / (2 * alpha(o)) * ...
                  hypergeom([0.5 * (3 - 2i * sqrt(C_m(m)) + sqrt(1 - 4 * C_m(m))), ...
                             0.5 * (3 + 2i * sqrt(C_m(m)) + sqrt(1 - 4 * C_m(m)))], ...
                            2 + sqrt(1 - 4 * C_m(m)), ...
                            -(1 - alpha(o)) / alpha(o));

phi_2_0(i) =  -lambda(o) / (2 * alpha(o)) * ...
                  hypergeom([0.5 * (3 - 2i * sqrt(C_m(m)) - sqrt(1 - 4 * C_m(m))), ...
                             0.5 * (3 + 2i * sqrt(C_m(m)) - sqrt(1 - 4 * C_m(m)))], ...
                            2 - sqrt(1 - 4 * C_m(m)), ...
                            -(1 - alpha(o)) / alpha(o));

% Define E and M
E(i) =  lambda(o) * exp(z(i) * lambda(o)) * 0.5 * (1 + sqrt(1 - 4 * C_m(m))) .* ...
               (exp(z(i) * lambda(o)) - alpha(o)).^ (0.5 * (-1 + sqrt(1 - 4 * C_m(m)))) .* F_1_z(i) + ...
               (exp(z(i) * lambda(o)) - alpha(o)).^ (0.5 * (1 + sqrt(1 - 4 * C_m(m)))) .* phi_1_z(i);

M(i) =  lambda(o) * exp(z(i) * lambda(o)) * 0.5 * (1 - sqrt(1 - 4 * C_m(m))) .* ...
               (exp(z(i) * lambda(o)) - alpha(o)).^ (0.5 * (-1 - sqrt(1 - 4 * C_m(m)))) .* F_2_z(i) + ...
               (exp(z(i) * lambda(o)) - alpha(o)).^ (0.5 * (1 - sqrt(1 - 4 * C_m(m)))) .* phi_2_z(i);

% Define D
D(i) =  lambda(o) * sqrt(1 - 4 * C_m(m)) / (1 - alpha(o)) + ...
            (F_2_0(i) * phi_1_0(i) - F_1_0(i) * phi_2_0(i)) / (F_1_0(i) * F_2_0(i));

% Define the mode shape Phi
Phi(i) =  ((1 - alpha(o)).^ (-0.5 * (1 + sqrt(1 - 4 * C_m(m)))) / F_1_0(i)) .* E(i) / D(i) - ...
                 ((1 - alpha(o)).^ (-0.5 * (1 - sqrt(1 - 4 * C_m(m)))) / F_2_0(i)) .* M(i) / D(i);



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
        
        if o == length(lambda) % For the last one, use black color and solid line
            plot(mode_shapes{o}{m}, z/H, 'k-');
        else
            plot(mode_shapes{o}{m}, z/H);
        end  

    end
    %title(['Mode Shape ' num2str(m)]);
    xlabel(['\Phi_' num2str(m)]);
    ylabel('z/H');
    xlim([-1.5 1.5]);
   
   legendEntries = arrayfun(@(o) sprintf('\\lambdaH = %.2f, \\alpha = %.2f', ...
    lambda(o)*H, alpha(o)), 1:length(lambda), 'UniformOutput', false);
   set(gca, 'YDir', 'reverse');
   legend(legendEntries, 'Location', 'Best','FontSize',10);



end

%% Section For Profile Reversal
H=50;

z=[];
for k = 0:0.2:H
        z = [z k];
end

lambda_Reversal=[8/H 5/H 3/H 10/H 1/H];
alpha_Reversal=[-0.5 -0.9 -1 -3 1E-5];


mode_shapes_Reversal = cell(length(lambda_Reversal), 1);
C_m_Values_Reversal = zeros(1, 3);
C_m_Values_lambda_Reversal=cell(length(lambda_Reversal), 1);
mode_shapes_for_lambda_Reversal = cell(3, 1);



for o_Reversal=1:length(lambda_Reversal)
   
    if o_Reversal ==1
        Iguess_Reversal=[0.0391 0.3695 1.0480];
    elseif o_Reversal ==2   
        Iguess_Reversal=[0.1073	1.1057 3.1656];
    elseif o_Reversal ==3
        Iguess=[0.25 3.8 7.8];
    elseif o_Reversal ==4
        Iguess_Reversal=[0.0256 0.2625 0.7752];
    elseif o_Reversal ==5
        Iguess_Reversal=[10 30 50];
    end

    for m_Reversal = 1:3 
             
   char_eq_Reversal = @(C_m_Reversal) real((lambda_Reversal(o_Reversal)*exp(lambda_Reversal(o_Reversal)*H))*((((exp(lambda_Reversal(o_Reversal)*H) - alpha_Reversal(o_Reversal))^((-1 + sqrt(1 - 4*C_m_Reversal))/2))*((1 + sqrt(1 - 4*C_m_Reversal))/2)*hypergeom([(1 - 2*1i*sqrt(C_m_Reversal) + sqrt(1 - 4*C_m_Reversal))/2, (1 + 2*1i*sqrt(C_m_Reversal) + sqrt(1 - 4*C_m_Reversal))/2], 1 + sqrt(1 - 4*C_m_Reversal), -((exp(lambda_Reversal(o_Reversal)*H) - alpha_Reversal(o_Reversal))/alpha_Reversal(o_Reversal))))...
    -((exp(lambda_Reversal(o_Reversal)*H) - alpha_Reversal(o_Reversal))^((-1 -sqrt(1 - 4*C_m_Reversal))/2))*((1 - sqrt(1 - 4*C_m_Reversal))/2)*((1 - alpha_Reversal(o_Reversal))^sqrt(1 - 4*C_m_Reversal))*(hypergeom([(1 - 2*1i*sqrt(C_m_Reversal) + sqrt(1 - 4*C_m_Reversal))/2,  (1 + 2*1i*sqrt(C_m_Reversal) + sqrt(1 - 4*C_m_Reversal))/2], 1+ sqrt(1 - 4*C_m_Reversal), -((1 - alpha_Reversal(o_Reversal))/alpha_Reversal(o_Reversal)))*hypergeom([(1 - 2*1i*sqrt(C_m_Reversal) - sqrt(1 - 4*C_m_Reversal))/2, (1 + 2*1i*sqrt(C_m_Reversal) - sqrt(1 - 4*C_m_Reversal))/2], 1 - sqrt(1 - 4*C_m_Reversal), -((exp(lambda_Reversal(o_Reversal)*H) - alpha_Reversal(o_Reversal))/alpha_Reversal(o_Reversal)))/hypergeom([(1 - 2*1i*sqrt(C_m_Reversal) - sqrt(1 - 4*C_m_Reversal))/2, (1 + 2*1i*sqrt(C_m_Reversal) - sqrt(1 - 4*C_m_Reversal))/2], 1 - sqrt(1 - 4*C_m_Reversal), -((1 - alpha_Reversal(o_Reversal))/alpha_Reversal(o_Reversal))))) ...
      + (((exp(lambda_Reversal(o_Reversal)*H) - alpha_Reversal(o_Reversal))^((1 + sqrt(1 - 4*C_m_Reversal))/2))*((- lambda_Reversal(o_Reversal)*exp(lambda_Reversal(o_Reversal)*H))/(2*alpha_Reversal(o_Reversal)))* hypergeom([(3 - 2*1i*sqrt(C_m_Reversal) + sqrt(1 - 4*C_m_Reversal))/2, (3 + 2*1i*sqrt(C_m_Reversal) + sqrt(1 - 4*C_m_Reversal))/2], 2 + sqrt(1 - 4*C_m_Reversal), -((exp(lambda_Reversal(o_Reversal)*H) - alpha_Reversal(o_Reversal))/alpha_Reversal(o_Reversal)))-((exp(lambda_Reversal(o_Reversal)*H) - alpha_Reversal(o_Reversal))^((1 - sqrt(1 - 4*C_m_Reversal))/2))*((1 - alpha_Reversal(o_Reversal))^sqrt(1 - 4*C_m_Reversal))*(hypergeom([(1 - 2*1i*sqrt(C_m_Reversal) + sqrt(1 - 4*C_m_Reversal))/2,  (1 + 2*1i*sqrt(C_m_Reversal) + sqrt(1 - 4*C_m_Reversal))/2], 1+ sqrt(1 - 4*C_m_Reversal), -((1 - alpha_Reversal(o_Reversal))/alpha_Reversal(o_Reversal)))*((- lambda_Reversal(o_Reversal)*exp(lambda_Reversal(o_Reversal)*H))/(2*alpha_Reversal(o_Reversal)))* hypergeom([(3 - 2*1i*sqrt(C_m_Reversal) - sqrt(1 - 4*C_m_Reversal))/2, (3 + 2*1i*sqrt(C_m_Reversal) - sqrt(1 - 4*C_m_Reversal))/2], 2 - sqrt(1 - 4*C_m_Reversal), -((exp(lambda_Reversal(o_Reversal)*H) - alpha_Reversal(o_Reversal))/alpha_Reversal(o_Reversal)))/hypergeom([(1 - 2*1i*sqrt(C_m_Reversal) - sqrt(1 - 4*C_m_Reversal))/2, (1 + 2*1i*sqrt(C_m_Reversal) - sqrt(1 - 4*C_m_Reversal))/2], 1 - sqrt(1 - 4*C_m_Reversal), -((1 - alpha_Reversal(o_Reversal))/alpha_Reversal(o_Reversal))))));
   
    options_Reversal = optimset('Display', 'iter','TolFun', 1e-18, 'TolX', 1e-18);
    C_m_Reversal(m_Reversal) = fzero(char_eq_Reversal, Iguess_Reversal(m_Reversal), options_Reversal); % Start with an initial guess of 1
     
    C_m_Values_Reversal(m_Reversal) =C_m_Reversal(m_Reversal);  
       
 parfor i_Reversal = 1:length(z)
     % Define hypergeometric functions F_1 and F_2 at z and 0
F_1_z(i_Reversal) =  hypergeom([0.5 * (1 - 2i * sqrt(C_m_Reversal(m_Reversal)) + sqrt(1 - 4 * C_m_Reversal(m_Reversal))), ...
                              0.5 * (1 + 2i * sqrt(C_m_Reversal(m_Reversal)) + sqrt(1 - 4 * C_m_Reversal(m_Reversal)))], ...
                             1 + sqrt(1 - 4 * C_m_Reversal(m_Reversal)), ...
                             -(exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) - alpha_Reversal(o_Reversal)) / alpha_Reversal(o_Reversal));

F_2_z(i_Reversal)=  hypergeom([0.5 * (1 - 2i * sqrt(C_m_Reversal(m_Reversal)) - sqrt(1 - 4 * C_m_Reversal(m_Reversal))), ...
                              0.5 * (1 + 2i * sqrt(C_m_Reversal(m_Reversal)) - sqrt(1 - 4 * C_m_Reversal(m_Reversal)))], ...
                             1 - sqrt(1 - 4 * C_m_Reversal(m_Reversal)), ...
                             -(exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) - alpha_Reversal(o_Reversal)) / alpha_Reversal(o_Reversal));

F_1_0(i_Reversal) =  hypergeom([0.5 * (1 - 2i * sqrt(C_m_Reversal(m_Reversal)) + sqrt(1 - 4 * C_m_Reversal(m_Reversal))), ...
                          0.5 * (1 + 2i * sqrt(C_m_Reversal(m_Reversal)) + sqrt(1 - 4 * C_m_Reversal(m_Reversal)))], ...
                         1 + sqrt(1 - 4 * C_m_Reversal(m_Reversal)), ...
                         -(1 - alpha_Reversal(o_Reversal)) / alpha_Reversal(o_Reversal));

F_2_0(i_Reversal) = hypergeom([0.5 * (1 - 2i * sqrt(C_m_Reversal(m_Reversal)) - sqrt(1 - 4 * C_m_Reversal(m_Reversal))), ...
                          0.5 * (1 + 2i * sqrt(C_m_Reversal(m_Reversal)) - sqrt(1 - 4 * C_m_Reversal(m_Reversal)))], ...
                         1 - sqrt(1 - 4 * C_m_Reversal(m_Reversal)), ...
                         -(1 - alpha_Reversal(o_Reversal)) / alpha_Reversal(o_Reversal));

% Define hypergeometric functions phi_1 and phi_2 at z and 0
phi_1_z(i_Reversal) =  -lambda_Reversal(o_Reversal) * exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) / (2 * alpha_Reversal(o_Reversal)) * ...
                      hypergeom([0.5 * (3 - 2i * sqrt(C_m_Reversal(m_Reversal)) + sqrt(1 - 4 * C_m_Reversal(m_Reversal))), ...
                                 0.5 * (3 + 2i * sqrt(C_m_Reversal(m_Reversal)) + sqrt(1 - 4 * C_m_Reversal(m_Reversal)))], ...
                                2 + sqrt(1 - 4 * C_m_Reversal(m_Reversal)), ...
                                -(exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) - alpha_Reversal(o_Reversal)) / alpha_Reversal(o_Reversal));

phi_2_z(i_Reversal) = -lambda_Reversal(o_Reversal) * exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) / (2 * alpha_Reversal(o_Reversal)) * ...
                      hypergeom([0.5 * (3 - 2i * sqrt(C_m_Reversal(m_Reversal)) - sqrt(1 - 4 * C_m_Reversal(m_Reversal))), ...
                                 0.5 * (3 + 2i * sqrt(C_m_Reversal(m_Reversal)) - sqrt(1 - 4 * C_m_Reversal(m_Reversal)))], ...
                                2 - sqrt(1 - 4 * C_m_Reversal(m_Reversal)), ...
                                -(exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) - alpha_Reversal(o_Reversal)) / alpha_Reversal(o_Reversal));

phi_1_0(i_Reversal) = -lambda_Reversal(o_Reversal) / (2 * alpha_Reversal(o_Reversal)) * ...
                  hypergeom([0.5 * (3 - 2i * sqrt(C_m_Reversal(m_Reversal)) + sqrt(1 - 4 * C_m_Reversal(m_Reversal))), ...
                             0.5 * (3 + 2i * sqrt(C_m_Reversal(m_Reversal)) + sqrt(1 - 4 * C_m_Reversal(m_Reversal)))], ...
                            2 + sqrt(1 - 4 * C_m_Reversal(m_Reversal)), ...
                            -(1 - alpha_Reversal(o_Reversal)) / alpha_Reversal(o_Reversal));

phi_2_0(i_Reversal) =  -lambda_Reversal(o_Reversal) / (2 * alpha_Reversal(o_Reversal)) * ...
                  hypergeom([0.5 * (3 - 2i * sqrt(C_m_Reversal(m_Reversal)) - sqrt(1 - 4 * C_m_Reversal(m_Reversal))), ...
                             0.5 * (3 + 2i * sqrt(C_m_Reversal(m_Reversal)) - sqrt(1 - 4 * C_m_Reversal(m_Reversal)))], ...
                            2 - sqrt(1 - 4 * C_m_Reversal(m_Reversal)), ...
                            -(1 - alpha_Reversal(o_Reversal)) / alpha_Reversal(o_Reversal));

% Define E and M
E_Reversal(i_Reversal) =  lambda_Reversal(o_Reversal) * exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) * 0.5 * (1 + sqrt(1 - 4 * C_m_Reversal(m_Reversal))) .* ...
               (exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) - alpha_Reversal(o_Reversal)).^ (0.5 * (-1 + sqrt(1 - 4 * C_m_Reversal(m_Reversal)))) .* F_1_z(i_Reversal) + ...
               (exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) - alpha_Reversal(o_Reversal)).^ (0.5 * (1 + sqrt(1 - 4 * C_m_Reversal(m_Reversal)))) .* phi_1_z(i_Reversal);

M_Reversal(i_Reversal) =  lambda_Reversal(o_Reversal) * exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) * 0.5 * (1 - sqrt(1 - 4 * C_m_Reversal(m_Reversal))) .* ...
               (exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) - alpha_Reversal(o_Reversal)).^ (0.5 * (-1 - sqrt(1 - 4 * C_m_Reversal(m_Reversal)))) .* F_2_z(i_Reversal) + ...
               (exp(z(i_Reversal) * lambda_Reversal(o_Reversal)) - alpha_Reversal(o_Reversal)).^ (0.5 * (1 - sqrt(1 - 4 * C_m_Reversal(m_Reversal)))) .* phi_2_z(i_Reversal);

% Define D
D_Reversal(i_Reversal) =  lambda_Reversal(o_Reversal) * sqrt(1 - 4 * C_m_Reversal(m_Reversal)) / (1 - alpha_Reversal(o_Reversal)) + ...
            (F_2_0(i_Reversal) * phi_1_0(i_Reversal) - F_1_0(i_Reversal) * phi_2_0(i_Reversal)) / (F_1_0(i_Reversal) * F_2_0(i_Reversal));

% Define the mode shape Phi
Phi_Reversal(i_Reversal) =  ((1 - alpha_Reversal(o_Reversal)).^ (-0.5 * (1 + sqrt(1 - 4 * C_m_Reversal(m_Reversal)))) / F_1_0(i_Reversal)) .* E_Reversal(i_Reversal) / D_Reversal(i_Reversal) - ...
                 ((1 - alpha_Reversal(o_Reversal)).^ (-0.5 * (1 - sqrt(1 - 4 * C_m_Reversal(m_Reversal)))) / F_2_0(i_Reversal)) .* M_Reversal(i_Reversal) / D_Reversal(i_Reversal);


 end  


 mode_shapes_for_lambda_Reversal{m_Reversal} = (Phi_Reversal); % Store mode shape for this mode
    end
    % Store mode shapes for the current n in the cell array
    mode_shapes_Reversal{o_Reversal} = mode_shapes_for_lambda_Reversal;
    C_m_Values_lambda_Reversal{o_Reversal}=C_m_Values_Reversal;
end


figure(2);
for m_Reversal = 1:3
   
    subplot(1, 3, m_Reversal); % Adjust the subplot position
    hold on;

    for o_Reversal = 1:length(lambda_Reversal)
        if o_Reversal == length(lambda_Reversal) % For the last one, use black color and solid line
            plot(mode_shapes_Reversal{o_Reversal}{m_Reversal}, z/H, 'k-');
        else
            plot(mode_shapes_Reversal{o_Reversal}{m_Reversal}, z/H);
        end  

    end
   
    xlabel(['\Phi_' num2str(m_Reversal)]);
    ylabel('z/H');
    xlim([-1.5 1.5]);
    set(gca, 'YDir', 'reverse');
   legendEntries = arrayfun(@(o) sprintf('\\lambdaH = %.2f, \\alpha = %.2f', ...
    lambda_Reversal(o)*H, alpha_Reversal(o)), 1:length(lambda_Reversal), 'UniformOutput', false);
   legend(legendEntries, 'Location', 'Best','FontSize',10);

end




toc
