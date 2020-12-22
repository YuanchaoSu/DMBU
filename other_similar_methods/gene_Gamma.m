function [y] = gene_Gamma(MPlus,theta)
R=size(MPlus,2);
a=theta(1:R);
u=1;
y=MPlus*a;

for i=1:R-1
    for j=i+1:R
        y=y+ theta(R+u)*a(i)*a(j)*MPlus(:,i).*MPlus(:,j);
        u=u+1;
    end
end
