data = load("AIdata22.dat", "r");
in = data(:,1);
out = data(:,2);

k = 0;
mu = 1000;
lmbd = 5*mu;
parent = zeros(lmbd,3);
pop = zeros(lmbd*2,3);
sigma = zeros(lmbd,3);
fit1 = zeros(lmbd,1);
fit2 = zeros(lmbd,1);
fit3 = zeros(lmbd*2,1);
fit3_i = zeros(lmbd*2,1);
fit3_delete = zeros(lmbd,1);
y = zeros(101,1);

tau1 = 1/(sqrt(2*6));
tau2 = 1/(sqrt(2*sqrt(6)));

% creating the initial population
for i = 1:lmbd
    parent(i,1) = -10 + (20).*rand(); % a
    parent(i,2) = -10 + (20).*rand(); % b
    parent(i,3) = -10 + (20).*rand(); % c
end

% creating random sigmas
for i = 1:lmbd
    sigma(i,1) = 0.0001 + (9.9998).*rand(); % sig_a
    sigma(i,2) = 0.0001 + (9.9998).*rand(); % sig_b
    sigma(i,3) = 0.0001 + (9.9998).*rand(); % sig_c
end

% iterations
while k < 500
    
    k = k + 1;
    
    for i = 1:lmbd
        % throwing parents into the whole population
        pop(i,1) = parent(i,1);
        pop(i,2) = parent(i,2);
        pop(i,3) = parent(i,3);
        
        % fitnesses of parents before mutation
        fit1(i) = 0;
        for j = 1:101
            ohat = parent(i,1)*(in(j)^2-parent(i,2)*cos(parent(i,3)*pi*in(j)));
            fit1(i) = fit1(i) + (out(j) - ohat)^2;
        end
    end
    
    % mutation (creating offsprings)
    for i = 1:lmbd
        a = parent(i,1) + normrnd(0, sigma(i,1));
        b = parent(i,2) + normrnd(0, sigma(i,2));
        c = parent(i,3) + normrnd(0, sigma(i,3));
        
        r1 = normrnd(0, tau2);
        sigma(i,1) = sigma(i,1)*exp(normrnd(0, tau1))*exp(r1);
        sigma(i,2) = sigma(i,2)*exp(normrnd(0, tau1))*exp(r1);
        sigma(i,3) = sigma(i,3)*exp(normrnd(0, tau1))*exp(r1);
        
        % fitnesses after mutation
        fit2(i) = 0;
        for j = 1:101
            ohat = a*(in(j)^2-b*cos(c*pi*in(j)));
            fit2(i) = fit2(i) + (out(j) - ohat)^2;
        end
        
        % throwing offsping into the whole population
        pop(i+lmbd,1) = a;
        pop(i+lmbd,2) = b;
        pop(i+lmbd,3) = c;
    end
    
    % fitness of the whole population
    for i = 1:lmbd*2
        fit3(i) = 0;
        for j = 1:101
            ohat = pop(i,1)*(in(j)^2-pop(i,2)*cos(pop(i,3)*pi*in(j)));
            fit3(i) = fit3(i) + (out(j) - ohat)^2;
        end
    end
    
    % find the best fitness of the whole population
    [fit3, fit3_i] = sort(fit3);
    fit3_delete = fit3_i(1:lmbd,:);
    
    % creating new parents
    for i = 1:lmbd
        parent(i,1) = pop(fit3_delete(i),1);
        parent(i,2) = pop(fit3_delete(i),2);
        parent(i,3) = pop(fit3_delete(i),3);
    end
    
    fit1 = sort(fit1);
    fit2 = sort(fit2);
    best_f1 = fit1(1);
    best_f2 = fit2(1);
    
    % stop criterion
    if abs(best_f1 - best_f2) < 0.00001
        break
    end
     
end

% find a function of minimum
for i=1:101
    y(i)= parent(1,1)*(in(i)^2-parent(1,2)*cos(parent(1,3)*pi*in(i)));
end

% printing outputs
fprintf('k=%1d err=%1.4f\n', k, fit3(1)/101);
fprintf('a=%1.4f b=%1.4f c=%1.4f\n', parent(1,1), parent(1,2), parent(1,3));

% plotting
plot(in, y, in, out);
legend('Optimization result', 'Accurate result');