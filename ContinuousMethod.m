function X = ContinuousMethod(A, N, K)
%% Initialization
X = ones(N,K)/K;
eta = 1e-6;
% Scale = N*K;
Scale = 1; 
Fpre = inf;
zeta = -1;
minStep = 1e-3;
maxStep = 1e-1;
stepL = minStep;
Const = abs(fl(X, A, zeta));
while zeta <= 1
%     zeta
    for FW_i = 1:100
        D = GF(X, A, zeta);
        if sum(sum(abs(D))) < Scale*eta
            continue;
        end
        
        %line search     X = X + alpha*(Y-X);
        Y = MaxY(D);
        if sum(sum((Y-X).*D)) < 0
            break;
        end
        
        XPre = X;
        %divided method
        Xleft = X;
        Xright = Y;
        deltaX = (Xright - Xleft);
        for i = 1:10
            X1 = Xleft + 0.5*deltaX;
            X2 = Xleft + (0.5 + eta)*deltaX;
            J1 = fl(X1, A, zeta);
            J2 = fl(X2, A, zeta);
            if J1 < J2
               Xleft = X1;
            else
               Xright = X1;
            end
            deltaX = Xright - Xleft;
        end 
        X = Xleft;
        if sum(sum(abs(X-XPre)))< Scale*eta
           break; 
        end
        
        Jnew = fl(X, A, zeta);
        Dnew = GF(X, A, zeta);
        if sum(sum((Y-X).*Dnew))< Scale*eta*abs(Jnew-sum(sum((Y-X).*Dnew)))
            break;
        end       
     
    end
    if sum(sum(abs(X - MaxY(X)))) < Scale*eta
        break;
    end
    
    Fcur = fl(X, A, zeta);
    if (abs(Fcur - Fpre) < Scale*Const*eta) | (Fcur < Fpre)
%     if abs(Fcur - Fpre) < Scale*Const*eta
        stepL = stepL * 2;
    else
        stepL = stepL / 2;
    end
    if stepL < minStep
        stepL = minStep;
    elseif stepL>maxStep
        stepL = maxStep;
    end
    Fpre = Fcur;
    if zeta ==1
        break;
    end
%     stepL
    zeta = zeta + stepL;
    if zeta > 1
        zeta = 1;
    end
end
X = MaxY(X);
end
function f = fl(X, B, gamma)
    if gamma > 0
        f = gamma * trace(X'*X) + (1-gamma) * trace(X'*B*X*pinv(X'*X)) ;
    else
        f = gamma * trace(X'*X) + (1+gamma) * trace(X'*B*X*pinv(X'*X));
    end
end
function g = GF(X, B, gamma)
    if gamma > 0
%         g = gamma*2*X + (1-gamma) * (eye(size(X,1))-X*pinv(X'*X)*X')*(B+B')*X*pinv(X'*X);
        g = gamma*2*X + (1-gamma) * (sparse(1:size(X,1), 1:size(X,1), 1, size(X,1), size(X,1))*(B+B')*X*pinv(X'*X)-X*pinv(X'*X)*(X'*(B+B')*X*pinv(X'*X))); % Sparse
    else
%         g = gamma*2*X + (1+gamma) * (eye(size(X,1))-X*pinv(X'*X)*X')*(B+B')*X*pinv(X'*X);
        g = gamma*2*X + (1+gamma) * (sparse(1:size(X,1), 1:size(X,1), 1, size(X,1), size(X,1))*(B+B')*X*pinv(X'*X)-X*pinv(X'*X)*(X'*(B+B')*X*pinv(X'*X))); % Sparse
    end
end