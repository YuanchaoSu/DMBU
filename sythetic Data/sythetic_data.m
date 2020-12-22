function [Y,X,M] = sythetic_data(MaxPurity,SNR,p)
% Linear unmixing model
% rand('seed',5);
load('USGS_1995_Library')
load True_M
[bands,n_materiais]=size(datalib);
% select randomly
sel_mat = 4+randperm(n_materiais-4);
sel_mat = sel_mat(1:p);
% M = datalib(:,sel_mat);
M = True_M;

dimension = 32;
label = ones((dimension/8)^2,1);
num = floor(length(label)/p);
  for i=1:p-1
      label((i-1)*num+1:i*num) = (i+1); 
  end
  ridx = randperm(length(label));
  label = label(ridx)';
  label = reshape(label,dimension/8,dimension/8);
  X = zeros(dimension,dimension,p);
  img = zeros(dimension,dimension);
  for i=1:dimension
      for j=1:dimension
          for cls = 1:p
              if label(floor((i-1)/8)+1,floor((j-1)/8)+1) == cls
                 tmp = zeros(p,1);
                 tmp(cls) = 1;
                 X(i,j,:) = tmp;
                 img(i,j) = p;
              end
          end
      end 
  end
  win = 7;
  H = ones(win,win)/(win*win);
  for i=1:p
      X(:,:,i) = filter2(H,X(:,:,i));
  end
  X = X(ceil(win/2):end-floor(win/2),ceil(win/2):end-floor(win/2),:);
  [m,n,p] = size(X);
  X = reshape(X,m*n,p)';
 %¡¡purity
  Index = ceil(find(X>MaxPurity)/p);
  X(:,Index) = 1/p*ones(p,1)*ones(1,length(Index));

   HIM = reshape((M*X)',m,n,bands);    
   Y = reshape(HIM,[m*n,bands,1]);
 if SNR <100
% add noise
   variance = sum(Y(:).^2)/10^(SNR/10)/m/m/bands;
   n = sqrt(variance)*randn([bands m*n]);
   Y = Y' + n;
 else
   Y=Y';
 end
end