% Assignment 6 - Simulations of Isotropic Dendritic Growth
Gridsize = 300;
T_old = zeros(Gridsize,Gridsize); % Temperature matrix
p_old = zeros(Gridsize,Gridsize); % Phase field matrix

%Parameter values
dt = 0.0002;  % Time discretization
dx = 0.03; % Space discretization along x
dy = 0.03; % Space discretization along y
%t = 0;  % Total time taken for simulation
K = 1.6;  % Latent heat
eps = 0.01;
tou = 0.0003;
Te = 1; % Equilibrium temperature
gamma = 10;
alpha = 0.9;
a = 0.01; % Amplitude of noise term


% Initialization of boundaries
for i = 148:152
    for j = 148:152
    p_old(i,j) = 1; % Makes a portion of the domain on one side solid
    end 
end

T_old(:) = 0;

T_new = zeros(Gridsize,Gridsize);
p_new = p_old;

%%% Solving for each grid iteratively using explicit scheme %%%
for t = 1:50000
    
    for i = 2:Gridsize-1
        for j = 2:Gridsize-1
            m = (alpha/pi)*atan(gamma*(Te-T_old(i,j)));
            if (p_old(i,j)>=0 && p_old(i,j)<=0.5)
 % Noise Term is added here because at p = 0.5, maximum is reached. 
X1 = rand(1) - 0.5; 
p_new(i,j) = p_old(i,j) + eps^2*dt/(tou*dx^2)*(p_old(i+1,j) + ...
             p_old(i-1,j) + p_old(i,j+1) + p_old(i,j-1) - ...
             4*p_old(i,j)) + ...
             dt/tou*p_old(i,j)*(1-p_old(i,j))*(p_old(i,j)+m-0.5) + ...
             a*p_old(i,j)*(1-p_old(i,j))*X1;
            else
p_new(i,j) = p_old(i,j) + eps^2*dt/(tou*dx^2)*(p_old(i+1,j) + ...
             p_old(i-1,j) + p_old(i,j+1) +p_old(i,j-1) - ...
             4*p_old(i,j)) +...
             dt/tou*p_old(i,j)*(1-p_old(i,j))*(p_old(i,j)+m-0.5);
            end                
T_new(i,j) = T_old(i,j)+ dt/dx^2*(T_old(i+1,j)+ T_old(i-1,j)+ ...
             T_old(i,j+1)+ T_old(i,j-1) -4*T_old(i,j)) + ...
             K*(p_new(i,j)-p_old(i,j));
        end
    end
    
    % Giving no flux boundary conditions
    p_new(1,:) = p_new(2,:);
    p_new(Gridsize,:) = p_new(Gridsize-1,:);
    p_new(:,1) = p_new(:,2);
    p_new(:,Gridsize) = p_new(:,Gridsize-1);
    T_new(:,1) = T_new(:,2);
    T_new(:,Gridsize) = T_new(:,Gridsize-1);
    T_new(1,:) = T_new(2,:);
    T_new(Gridsize,:) = T_new(Gridsize-1,:) - 100;
    T_old = T_new;
    p_old = p_new;
    % Plotting at required time intervals
    if (mod(t,100)==0)
        pcolor(p_old);
title("Isotropic dendritic growth in phase field solidification"); 
        pause(1);
        str1 = ['isotropic',num2str(t)];
        print(str1,'-dpng');
    end
    
end
