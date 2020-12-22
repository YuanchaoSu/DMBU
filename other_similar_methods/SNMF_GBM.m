function [A,B,M,Gamma,RMSE] = SNMF_GBM(X,E,ite_max,A_init,B_init)

%------------------------------------------------------------------
% SEMI-NONNEGATIVE MATRIX FACTORIZATION FOR HS UNMIXING
%
% [A,B,M,Gamma] = SNMF_GBM(X,E,delta,ite_max,A_init,B_init)
% 
% INPUT
%       X               : HS image (L, N)
%       E               : Endmember matrix (L, R)
%       delta           : Parameter of sum to one constraint
%       ite_max         : Maximum number of iteration
%       A_init          : Initial abundance matrix (option)
%       B_init          : Initial interaction matrix (option)
%
% OUTPUT
%       A               : Abundance matrix (R, N)
%       B               : Interaction matrix (R*(R-1)/2, N)
%       M               : 2 times reflection endmember spectra (L, R*(R-1)/2)
%       Gamma           : Interaction coefficient matrix (R*(R-1)/2, N)
%
%------------------------------------------------------------------
initpara=0.001;
eps = 0.001;
delta = 0.1;
% delta = 10;
L = size(X,1);
N = size(X,2);
R = size(E,2);
M = zeros(L,R*(R-1)/2);
B_max = zeros(R*(R-1)/2,N);
Gamma = zeros(R*(R-1)/2,N);
RMSE = zeros(1,ite_max);
% Initialize M
u = 1;
for i = 1:R-1
    for j = i+1:R
        M(:,u) = E(:,i).*E(:,j);
        u = u+1;
    end
end

if nargin==4
    [A] = hyperFcls(X,E);
    % Initialize B_max
    for i = 1:N
        a = A(:,i);
        b = (triu(a*a')-diag(diag(a*a'))+tril(ones(size(a,1)))*(-1))';
        B_max(:,i) = b(b~=-1);
    end
    B = B_max*initpara;
end
if nargin==5
    A = A_init;
    % Initialize B_max
    for i = 1:N
        a = A(:,i);
        b = (triu(a*a')-diag(diag(a*a'))+tril(ones(size(a,1)))*(-1))';
        B_max(:,i) = b(b~=-1);
    end
    B = B_max*initpara;
end
if nargin==6
    A = A_init;
    % Initialize B_max
    for i = 1:N
        a = A(:,i);
        b = (triu(a*a')-diag(diag(a*a'))+tril(ones(size(a,1)))*(-1))';
        B_max(:,i) = b(b~=-1);
    end
    B = B_init;
end

E_1 = [E; delta*ones(1,R)];

% NMF
%disp('Iteration');
for i=1:ite_max
    B(:,sum(B)==0) = eps;
    % Update A (NMF update rule)
    X_1 = [X-M*B; delta*ones(1,N)]; % not always nonnegative
    XtEplus = (abs(X_1'*E_1)+X_1'*E_1)/2;
    XtEminus = (abs(X_1'*E_1)-X_1'*E_1)/2;
    EtEplus = (abs(E_1'*E_1)+E_1'*E_1)/2;
    EtEminus = (abs(E_1'*E_1)-E_1'*E_1)/2;
    A = (A'.*(((XtEplus+A'*EtEminus)./(XtEminus+A'*EtEplus)).^0.5))';
    % Update B (Semi-NMF update rule)
    M_2 = M;
    X_2 = X-E*A; % not always nonnegative
    XtMplus = (abs(X_2'*M_2)+X_2'*M_2)/2;
    XtMminus = (abs(X_2'*M_2)-X_2'*M_2)/2;
    MtMplus = (abs(M_2'*M_2)+M_2'*M_2)/2;
    MtMminus = (abs(M_2'*M_2)-M_2'*M_2)/2;
    B = (B'.*(((XtMplus+B'*MtMminus)./(XtMminus+B'*MtMplus)).^0.5))';
    B(B>B_max) = B_max(B>B_max);
    % Initialize B_max
    for k = 1:N
        a = A(:,k);
        b = (triu(a*a')-diag(diag(a*a'))+tril(ones(size(a,1)))*(-1))';
        B_max(:,k) = b(b~=-1);
    end
    RMSE(1,i) = (sum(sum((X-E*A-M*B).^2))/(L*N)).^0.5;
end
Gamma(B_max>0) = B(B_max>0)./B_max(B_max>0);


