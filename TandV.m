function [T,V] = TandV(x,y,gridx,Q,nu0)
% function x is measured on the equally spaced grid gridx
% nrow(x)=n

%%%%%%%%%%%%%%%%
% Please source package chebfun
%%%%%%%%%%%%%%%%
% addpath('/Users/jiajunt2/Dropbox/RKHS/rkhs_2021_05_25/chebfun-master')
nuq = nu0 + (1:Q) * (1-nu0) / Q;
nq = round(size(x,1)*nuq);

xori = x;
yori = y;

dgrid = gridx(2) - gridx(1);
lg = length(gridx);
xcov = x' * x /size(x,1);

%%%%%%%% compute eigenfunction
% num_max = min(ceil(n^(2/5)),ceil(log(n)^2));       %number of candidate eigenfunctions

v = floor(size(x,1)^0.4);
vv = v:(2*v);
CV = zeros(1,v+1);

for kk = 1:length(vv)
    num = vv(kk);
    
    phi = zeros(num,lg);
    % differential operator
    L = chebop(0,1);
    L.op = @(v) diff(v,4);
    L.lbc = @(v) [diff(v,2);diff(v,3)];       % left boundary condition
    L.rbc = @(v) [diff(v,2);diff(v,3)];        % right boundary condition
    
    B = chebop(0,1);
    K = @(x,s) F_cheb(x,s,xcov);                 % integral operator
    B.op = @(v) fred(K,v);   % fred:  Fredholm integral operator (Chebfun@blockFunction)
    [e,valtemp] = eigs(L,B,num);
    val= abs(real(diag(valtemp)));
    
    
    % Fetch eigenfunctions of the integro-differential equations
    % from Chebfun object "e"
    for k = 1:num
        fun = chebfun(e(:,k));
        temp = real(fun(gridx));
        phi(k,:) = temp'/sqrt(temp*xcov*temp')/dgrid;
    end
    
    Phi = zeros(num, num);
    for i = 1:num
        for j = 1:i
            Phi(i,j) = sum(phi(i,:).*phi(j,:))*dgrid;
            Phi(i,j) = Phi(j,i);
        end
    end
    Lambda = diag(val);
    
    % Five-fold cross-validation to select r
    for ii = 1:5
        foldsize = round(size(x,1)/5);
        foldid = (foldsize *(ii-1)+1):(foldsize*ii);
        xtrain = xori(setdiff(1:size(x,1),foldid),:);
        xtest = xori(foldid,:);
        ytrain = yori(setdiff(1:size(x,1),foldid),:);
        ytest = yori(foldid,:);
        
        % % GCV to select lambda
        temp = 1e-7;
        % Layout a grid for choosing lambda
        lambda_cand = temp*[0.1:0.1:3];
        gcv = zeros(1, length(lambda_cand));
        
        for q = 1:Q
            n = nq(q);
            x = xori(1:n,:);
            y = yori(1:n);
            Omega = x*phi'*dgrid;
            for k = 1:length(lambda_cand)
                tr_Hlambda = 0;
                Lambda = diag(val);
                temp2 = (Omega'*Omega+n*lambda_cand(k)*Lambda)\Omega';
                hatbeta = phi'*temp2*y;
                haty = x * hatbeta*dgrid;
                tr_Hlambda = tr_Hlambda+trace(Omega*temp2);
                gcv(k) = gcv(k)+sum((haty-y).^2)*dgrid/n/(1-tr_Hlambda/n)^2;
            end
        end
        [na,Index]=min(gcv);
        lambda = lambda_cand(Index);
        
        Omega = xtrain*phi'*dgrid;
        temp2 = (Omega'*Omega+n*lambda*Lambda)\Omega';
        hatbeta = phi'*temp2*ytrain;
        hatytest = xtest * hatbeta*dgrid;
        CV(ii) = CV(ii) + sum((ytest-hatytest).^2)*dgrid;
    end
    
end

[na,Index]=min(CV);
num = vv(Index);


phi = zeros(num,lg);
% differential operator
L = chebop(0,1);
L.op = @(v) diff(v,4);
L.lbc = @(v) [diff(v,2);diff(v,3)];       % left boundary condition
L.rbc = @(v) [diff(v,2);diff(v,3)];        % right boundary condition

B = chebop(0,1);
K = @(x,s) F_cheb(x,s,xcov);                 % integral operator
B.op = @(v) fred(K,v);   % fred:  Fredholm integral operator (Chebfun@blockFunction)
[e,valtemp] = eigs(L,B,num);
val= abs(real(diag(valtemp)));
for k = 1:num
    fun = chebfun(e(:,k));
    temp = real(fun(gridx));
    phi(k,:) = temp'/sqrt(temp*xcov*temp')/dgrid;
end

Phi = zeros(num, num);
for i = 1:num
    for j = 1:i
        Phi(i,j) = sum(phi(i,:).*phi(j,:))*dgrid;
        Phi(i,j) = Phi(j,i);
    end
end
Lambda = diag(val);

temp = 1e-7;

lambda_cand = temp*[0.1:0.1:3];
gcv = zeros(1, length(lambda_cand));

for q = 1:Q
    n = nq(q);
    x = xori(1:n,:);
    y = yori(1:n);
    Omega = x*phi'*dgrid;
    for k = 1:length(lambda_cand)
        tr_Hlambda = 0;
        Lambda = diag(val);
        temp2 = (Omega'*Omega+n*lambda_cand(k)*Lambda)\Omega';
        hatbeta = phi'*temp2*y;
        haty = x * hatbeta*dgrid;
        tr_Hlambda = tr_Hlambda+trace(Omega*temp2);
        gcv(k) = gcv(k)+sum((haty-y).^2)*dgrid/n/(1-tr_Hlambda/n)^2;
    end
end
[na,Index]=min(gcv);
lambda = lambda_cand(Index);

% Compute T and V
hatb = zeros(Q,num);

for q = 1:Q
    n = nq(q);
    x = xori(1:n,:);
    y = yori(1:n);
    Omega = x*phi'*dgrid;
    hatb(q,:) = (Omega'*Omega+n*lambda*Lambda)\Omega'*y;
end

T = hatb(Q,:)* Phi *hatb(Q,:)';
V = 0;
for q = 1:Q
    V = V+nuq(q)^4*(hatb(q,:)* Phi *hatb(q,:)'-T)^2;
end
V = sqrt((1-nu0)/Q*V);



end