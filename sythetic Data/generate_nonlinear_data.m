function [Y,M,X] = generate_nonlinear_data(varargin)
% [Y,M,X] = generate_nonlinear_data('MAXPURITY',MaxPurity,'SNR',SNR,'SHOW_LINEAR','yes','SHOW_NONLINEAR','yes');
% Input----------------------------------------------
% Recommended SNR = 100 to 20 dB.
% The data is design based on bilinear unmixied model.
% MaxPurity aim at removing pure pixels.
% Output---------------------------------------------
% Y: The generated nonlinear data
% M: extracted endmembers
% X: The estimated abundances
%% Start
show_linear = 'no';
show_nonlinear = 'no';
if (nargin-length(varargin)) ~= 0
    error('Wrong number of required parameters');
elseif (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
for i=1:2:(length(varargin)-1)
   switch upper(varargin{i})
          case 'MAXPURITY'
                MaxPurity = varargin{i+1};% m is the number of endmembers
          case 'SNR'
                SNR = varargin{i+1};% maxiterations is the maximum iteration.
          case 'SHOW_LINEAR'
                show_linear = varargin{i+1};%display
          case 'SHOW_NONLINEAR'
                show_nonlinear = varargin{i+1};%display
          otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
   end;
end

m = 3;
[Y_LMM,X,M] = sythetic_data(MaxPurity,SNR,m);
[y1,Up,~,~] = dataProj(Y_LMM,m,'proj_type','affine');
y1=Up'*y1;
I = 1;
J = 2;
K = 3;
E_I = eye(m);
v1 = E_I(:,I);
v2 = E_I(:,J);
v3 = E_I(:,K);
y1 = [v1 v2 v3]'*y1;

if strcmp(show_linear, 'yes')
    figure
    scatter3(y1(1,:),y1(2,:),y1(3,:),50,[0.5 0.5 0.5],'.');
    title('Linear Model');
end

Y=Y_LMM;
[~, N]=size(Y);
% runing nolinearity
for  num=1:N
     u=1;
     a=X(:,num);
     theta=[0.01 0.05 0.01]'; 
     R=m;
        for i=1:R-1
            for j=i+1:R
                Y(:,num)=Y(:,num) + theta(u)*a(i)*a(j)*M(:,i).*M(:,j);
                u=u+1;
            end
        end  
end
 [~, PCAresults] = pca(Y');
 IList = PCAresults(:, 1:3);
if strcmp(show_nonlinear, 'yes')
    figure;
    scatter3(IList(1:1:end,1), IList(1:1:end,2), IList(1:1:end,3),10);
    xlabel('PC 1','fontsize',15);
    ylabel('PC 2','fontsize',15);
    zlabel('PC 3','fontsize',15);
%     title('Nonlinear Model');
end
end
end