
clear

% Compute the quantiles of the pivotal statistic W
Q = 25;
nu0 = 1/2;
dt = (1-nu0) / Q;
nn = 10000;
ww = zeros(nn, 1);
nuq = nu0 + (1-nu0)/Q*(1:Q)';
for i = 1:nn
    b = zeros(Q, 1);
    b(1) = sqrt(nu0+dt) * randn;
    for j = 2:Q
        b(j) = b(j-1)+sqrt(dt)*randn;
    end
    %     temp = nuq .* (nuq .* b - b(end));
    temp = nuq .* b - b(end);
    temp = norm(temp)^2;
    temp = sqrt(temp/Q*(1-nu0));
    ww(i) = b(end) / temp;
end
Wq = quantile(ww, 0.95);
clear ww


Delta = [0.7 1 1.3];
delta = 0.1:0.1:2;

% Compute empirical rejection probabilities
SIMU = 50;
rj = zeros(length(Delta), length(delta));

gridx = 0:0.01:1;
lg = length(gridx);
dgrid = gridx(2)-gridx(1);
tildebeta = ones(1,lg);


% beta setting S1
f = sqrt(2)*cos((1:49)'.*gridx*pi);
for k=1:49
    tildebeta = tildebeta+4*f(k,:)/(k+1)^2*(-1)^k;
end
% standardize
tildebeta = tildebeta / sqrt(sum(tildebeta.^2) * dgrid);

n = 200;    % number of curves
for a = 1:length(delta)
    for d = 1:length(Delta)
        beta = tildebeta * sqrt(delta(a));
        
        for simu = 1:SIMU
            % Generate samples
            % X setting 1
            x = zeros(lg,n);
            theta = unifrnd(-1/sqrt(2),1/sqrt(2),n,1);
            z = normrnd(0,1,lg,n+1);
            eta = ones(lg, n+1).*z(1,:);
            for k=1:49
                eta = eta+f(k,:)'/(k+1).*z;
            end
            for i=1:n
                x(:,i) = eta(:,i) + theta(i)*eta(:,i+1);
            end

            
            % error setting
            xi = normrnd(0,1,n+2,1);
            epsilon = normrnd(0,1,n,1);
            upsilon = unifrnd(-1/sqrt(2),1/sqrt(2),n,2);
            for i=1:n
                epsilon(i) = xi(i) + upsilon(i,1)*xi(i+1) + upsilon(i,2)*xi(i+2);
            end
            % rescale errors
            c = sqrt(0.3*var(sqrt(sum(x'*beta') * dgrid))/var(epsilon));
            epsilon = c * epsilon;
            
            x = x';
            y = x * beta' *dgrid + epsilon;
            
            [T,V] = TandV(x,y,gridx,Q,nu0);
            rj(d,a) = rj(d,a) + (T>Wq*V+Delta(d));
        end
        rj(d,a) = rj(d,a)/SIMU;
        
    end
    
end


rj

