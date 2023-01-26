function Synthetic_exp

% 2D Demo
for ExpCircle = 1:1
    for iter = 1:1
        [ExpCircle iter]
        N = 200; 
        K = 2;
%         KNN = 40;
        KNN = 20;
        [Data X_GT] = GeneDataCircle(N,K);
               
        X = Data;
        Y = X;
        % MdlES = ExhaustiveSearcher(X)
        % IdxES = knnsearch(MdlES,Y,'K', 3)
        MdlKDT = KDTreeSearcher(X);
        IdxKDT = knnsearch(MdlKDT,Y,'K',KNN+1);
        if abs(IdxKDT(:,1) - (1:N)')~=0
           error('kNN Graph!') 
        end
        E = [];
        for i = 1:KNN
           insE = [IdxKDT(:,1) IdxKDT(:,i+1)];
           E= [E;insE];
        end
        G = sparse(E(:,1),E(:,2),ones(size(E,1),1),N,N);
        Dis = sqrt(sum((X(E(:,1),:) - X(E(:,2),:)).^2,2));
        Aff = exp(-Dis.^2/max(Dis.^2)/1);

        A = sparse(E(:,1),E(:,2),Aff,N,N);
        A = A+A';

        XX = ContinuousMethod(A, N, K);

        if 1
            figure;
            hold on;
            plot(Data(find(XX(:,1)),1),Data(find(XX(:,1)),2),'r.');
            plot(Data(find(XX(:,2)),1),Data(find(XX(:,2)),2),'b.');
            axis equal
            hold off;
        end

    end

end

% 3D Demo
for ExpGauss = 2:2
    for iter = 1:1
        [ExpGauss iter]
        N = 700; 
        K = 7;
        D = 3;
        KNN = 70;
        [Data X_GT] = GeneDataGauss(N, K, D);
 
        X = Data;
        Y = X;
        % MdlES = ExhaustiveSearcher(X)
        % IdxES = knnsearch(MdlES,Y,'K', 3)
        MdlKDT = KDTreeSearcher(X);
        IdxKDT = knnsearch(MdlKDT,Y,'K',KNN+1);
        if abs(IdxKDT(:,1) - (1:N)')~=0
           error('kNN Graph!') 
        end
        E = [];
        for i = 1:KNN
           insE = [IdxKDT(:,1) IdxKDT(:,i+1)];
           E= [E;insE];
        end
        G = sparse(E(:,1),E(:,2),ones(size(E,1),1),N,N);
        Dis = sqrt(sum((X(E(:,1),:) - X(E(:,2),:)).^2,2));
        Aff = exp(-Dis.^2/max(Dis.^2)/1);
        % Aff = -Dis.^2;
        A = sparse(E(:,1),E(:,2),Aff,N,N);
        A = A+A';

        XX = ContinuousMethod(A, N, K);

        if 1
            figure
            plot3(Data(find(XX(:,1)),1),Data(find(XX(:,1)),2),Data(find(XX(:,1)),3),'r.');
            hold on;
            plot3(Data(find(XX(:,2)),1),Data(find(XX(:,2)),2),Data(find(XX(:,2)),3),'b.');
            plot3(Data(find(XX(:,3)),1),Data(find(XX(:,3)),2),Data(find(XX(:,3)),3),'g.');
            plot3(Data(find(XX(:,4)),1),Data(find(XX(:,4)),2),Data(find(XX(:,4)),3),'c.');
            plot3(Data(find(XX(:,5)),1),Data(find(XX(:,5)),2),Data(find(XX(:,5)),3),'k.');
            plot3(Data(find(XX(:,6)),1),Data(find(XX(:,6)),2),Data(find(XX(:,6)),3),'y.');
            plot3(Data(find(XX(:,7)),1),Data(find(XX(:,7)),2),Data(find(XX(:,7)),3),'m.');
            axis equal
            hold off;
            
        end
        
    end

end

clc
clear
close all
end

function [Data GTX] = GeneDataCircle(N, K)
Data = [];
GTX = zeros(N,K);
SumN = 0;
for k = 1:K
    if k==K
        Nk = N-SumN; 
        SumN = N;
    else
        Nk = round(N*k/sum(1:K)); 
        SumN = SumN+Nk;
    end
    while 1
        rk = 1*k+0.1*randn(Nk,1)-0.7;
        if sum(rk<=0)==0 %prunning
            break
        end
    end
    
    thetak = rand(Nk,1)*2*pi;
    Pk = [rk.*cos(thetak) rk.*sin(thetak)];
    Data = [Data;Pk];
    GTX(SumN-Nk+1:SumN,k) = 1;
end
end

function [Data GTX] = GeneDataGauss(N, K, D)
% K <= 1331 
% D >= 3
Data = [];
GTX = zeros(N,K);
sumK = 0;

[X,Y,Z] = meshgrid(-5:5,-5:5,-5:5);
Grids = [X(:) Y(:) Z(:)];

[~,I] = sort(sum(Grids.^2,2),'ascend');
Centers = Grids(I,:)*6;

for k = 1:K
   Nk = round(N/K);
   if k == K
       Nk = N - sumK;
   end
   Centerk = [Centers(k,:) zeros(1,D-3)];
   Pk = ones(Nk,1)*Centerk + randn(Nk,D);
   Data = [Data; Pk];
   GTX(sumK+1:sumK+Nk,k)=1;
   sumK = sumK + Nk;
end

end
