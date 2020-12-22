function [alpha gam]=GBM_gradient(y,MPlus)

%------------------------------------------------------------------
% NONLINEAR SPECTRAL UNMIXING ASSUMING THE GENERALIZED BILINEAR MODEL
% by A. HALIMI, Y. ALTMANN, N. DOBIGEON and J.-Y. TOURNERET 
% IRIT/ENSEEIHT/TéSA - France - 2011
% 
% INPUT
%       y        : L x 1 pixel to be unmixed where L is the number 
%                  of spectral band
%       MPlus    : L x R endmember matrix involved in the mixture
%                  where R is the number of endmembers
%
% OUTPUT
%       alpha    : abundance vector estimated by the procedure
%       gam      : gamma vector estimated by the procedure
%
%------------------------------------------------------------------

% number of spectral bands x number of endmember spectra
[L R]=size(MPlus);


%% Initialization
Niter=100;
% Initialization with the FCLS algorithm (assuming the linear mixing model) 
% for the abundance vector
theta0=[hyperFcls(y,MPlus)' 0.01*ones(1,R*(R-1)/2)]';

THETA(:,1)=theta0;
ER(1)=20;
Gam_sq=0.01*ones(R,R);

%% Iterations
for t=2:Niter
y0=gene_Gamma(MPlus,theta0);

% Derivatives with respect to the R-1 first abundances
for r=1:R-1
d_alpha(:,r)= MPlus(:,r)-MPlus(:,R);
    for i=1:R-1
        if i~=r
           d_alpha(:,r) = d_alpha(:,r) + Gam_sq(r,i)*theta0(i)*MPlus(:,r).*MPlus(:,i)...
           - Gam_sq(r,R)*theta0(i)*MPlus(:,r).*MPlus(:,R)- Gam_sq(i,R)*theta0(i)*MPlus(:,i).*MPlus(:,R);
        end
    end
d_alpha(:,r) = d_alpha(:,r) + Gam_sq(r,R)*MPlus(:,r).*MPlus(:,R) ...
                - 2* Gam_sq(r,R)*theta0(r)*MPlus(:,r).*MPlus(:,R);     
end            
  

% Derivative with respect to gamma 
d_gamma=zeros(L,R,R);
for i=1:R-1
    for j=i+1:R
        d_gamma(:,i,j)=theta0(i)*theta0(j)*MPlus(:,i).*MPlus(:,j);
        df_gamma(i,j)=(y0-y)'*squeeze(d_gamma(:,i,j)); 
    end
end  

df_alpha=(y0-y)'*d_alpha;

% Modification of the gradient vector for move optimization
if sum(theta0(1:R-1))>0.99
    n=ones(R-1,1);
    if (df_alpha*n<0)
        df_alpha(R-1)=-sum(df_alpha(1:R-2));
    end
end

for r=1:R-1
if theta0(r)< 0.01
    if df_alpha(r)>0
        df_alpha(r)=0;
    end
end
end
      
% Determination of the bounds for the line search parameter
 if(df_alpha==zeros(1,R-1))
     lb=0;ub=0;
 else
     lb=min(theta0(1)/df_alpha(1)*(df_alpha(1)~=0),(theta0(1)-1)/df_alpha(1))*(df_alpha(1)~=0);
     ub=max(theta0(1)/df_alpha(1)*(df_alpha(1)~=0),(theta0(1)-1)/df_alpha(1))*(df_alpha(1)~=0);
     for r=1:R-1
         lb=max(lb,min(theta0(r)/df_alpha(r)*(df_alpha(r)~=0),(theta0(r)-1)/df_alpha(r))*(df_alpha(r)~=0));
         ub=min(ub,max(theta0(r)/df_alpha(r)*(df_alpha(r)~=0),(theta0(r)-1)/df_alpha(r))*(df_alpha(r)~=0));
         for j=r+1:R
            lb=max(lb, min(Gam_sq(r,j)/df_gamma(r,j)*(df_gamma(r,j)~=0),(Gam_sq(r,j)-1)/df_gamma(r,j))*(df_gamma(r,j)~=0));
            ub=min(ub, max(Gam_sq(r,j)/df_gamma(r,j)*(df_gamma(r,j)~=0),(Gam_sq(r,j)-1)/df_gamma(r,j))*(df_gamma(r,j)~=0));
         end
     end       

 end
lb=max(lb,min(sum(theta0(1:R-1))/sum(df_alpha(1:R-1))*(sum(df_alpha(1:R-1))~=0),(sum(theta0(1:R-1))-1)/sum(df_alpha(1:R-1)))*(sum(df_alpha(1:R-1))~=0));
ub=min(ub,max(sum(theta0(1:R-1))/sum(df_alpha(1:R-1))*(sum(df_alpha(1:R-1))~=0),(sum(theta0(1:R-1))-1)/sum(df_alpha(1:R-1)))*(sum(df_alpha(1:R-1))~=0));  

lb = max(lb,0);
ub = max(ub,0);
df=df_alpha;
u=1;
for i=1:R-1
    for j=i+1:R
        df(R-1+u)=df_gamma(i,j);
        u=u+1;
    end
end

% Line search procedure using the Golden section method
[lambda1]=golden_section( theta0,df,MPlus,y,lb,ub);

% Unknown parameters update
theta0(1:R-1)=theta0(1:R-1)-lambda1(1)*df(1:R-1)';
theta0(R)=1-sum(theta0(1:R-1));

for r=R+1:R*(R+1)/2
theta0(r)=theta0(r)-lambda1(1)*df(r-1);
end


THETA(:,t)=theta0;
DALPHA(t)=norm(THETA(:,t)-THETA(:,t-1));
ER(t)=norm(y-gene_Gamma(MPlus,THETA(:,t)))^2;
DER(t)=ER(t)-ER(t-1);
u=1;
for i=1:R-1
    for j=i+1:R
        Gam_sq(i,j)=theta0(R+u);
        u=u+1;
        Gam_sq(j,i)=Gam_sq(i,j);
    end
end

if(abs(DER(t)) < 10^(-5)) break;
end

end
alpha=theta0(1:R);
gam=theta0(R+1:R*(R+1)/2);