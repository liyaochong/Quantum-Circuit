function CaFluxplot(Ej,Ec,alpha,beta,kappa,sigma,Cc,Cr,Csh,FluxBias,nk,nl,nm,nlevels,Var)
%     example:CaFluxplot(148e9,3.29e9,0.613,0,0,0,4.19e-15,150e-15,0,0.5,5,10,2,20,'Ej');
tic;
hbar=1.054560652926899e-034;
h = hbar*2*pi;
e = 1.60217662e-19;        
                 
%    parpool(3);
%%
if strcmpi(Var,'FluxBias') 
  w1 = zeros(1,41);w2 = zeros(1,41);w3 = zeros(1,41);w4 = zeros(1,41);
  parfor i=1:41
      FluxBias = 0.48+(i-1)*0.001;
       [EL,~] = CaFluxQubit(Ej,Ec,alpha,beta,kappa,sigma,Cc,Cr,Csh,FluxBias,nk,nl,nm,nlevels);
       w1(i)=EL(3)-EL(1);w2(i)=EL(5)-EL(3);w3(i)=EL(5)-EL(1);w4(i)=(EL(5)-EL(1))/2;
%         w1(i)=(EL(3)-EL(1))/(2*pi);w2(i)=(EL(5)-EL(3))/(2*pi);w3(i)=(EL(5)-EL(1))/(2*pi);w4(i)=(EL(5)-EL(1))/(2*pi)/2;
   end
   I = zeros(1,40);
   for i=1:40
       I(i) = (w1(i+1)-w1(i))*(2*e);
   end   
%    delete(gcp('nocreate'));%关闭Pool
   x=-0.02:0.1/100:0.02;
   figure;
   plot(x,w1);hold on;plot(x,w2,'--');hold on;plot(x,w3,'-.');hold on;plot(x,w4,':');hold on;
   legend('\omega_{01}','\omega_{12}','\omega_{02}','\omega_{02}/2')
   xlabel('FluxBias');ylabel('频率v/Hz');
   figure;
   plot(x(1:40),I);hold on;
   xlabel('FluxBias');ylabel('Ip');
%%
elseif strcmpi(Var,'Ej') 
      w1 = zeros(1,41);w2 = zeros(1,41);w3 = zeros(1,41);w4 = zeros(1,41);
  parfor i=1:41
      Ej = (68+(i-1)*4)*10^9;
       [EL,~] = CaFluxQubit(Ej,Ec,alpha,beta,kappa,sigma,Cc,Cr,Csh,FluxBias,nk,nl,nm,nlevels);
       w1(i)=EL(3)-EL(1);w2(i)=EL(5)-EL(3);w3(i)=EL(5)-EL(1);w4(i)=(EL(5)-EL(1))/2;
%         w1(i)=(EL(3)-EL(1))/(2*pi);w2(i)=(EL(5)-EL(3))/(2*pi);w3(i)=(EL(5)-EL(1))/(2*pi);w4(i)=(EL(5)-EL(1))/(2*pi)/2;
  end
%    delete(gcp('nocreate'));%关闭Pool
   x=68:4:228;
   figure;
   plot(x,w1);hold on;plot(x,w2,'--');hold on;plot(x,w3,'-.');hold on;plot(x,w4,':');hold on;
   legend('\omega_{01}','\omega_{12}','\omega_{02}','\omega_{02}/2')
   xlabel('Ej');ylabel('频率v/Hz');
%%
elseif strcmpi(Var,'Ec') 
    w1 = zeros(1,41);w2 = zeros(1,41);w3 = zeros(1,41);w4 = zeros(1,41);
    parfor i=1:41
        Ec = (1.19+(i-1)*0.1)*10^9;
       [EL,~] = CaFluxQubit(Ej,Ec,alpha,beta,kappa,sigma,Cc,Cr,Csh,FluxBias,nk,nl,nm,nlevels);
       w1(i)=EL(3)-EL(1);w2(i)=EL(5)-EL(3);w3(i)=EL(5)-EL(1);w4(i)=(EL(5)-EL(1))/2;
%         w1(i)=(EL(3)-EL(1))/(2*pi);w2(i)=(EL(5)-EL(3))/(2*pi);w3(i)=(EL(5)-EL(1))/(2*pi);w4(i)=(EL(5)-EL(1))/(2*pi)/2;
    end
%    delete(gcp('nocreate'));%关闭Pool
    x=1.19:0.1:5.19;
    figure;
    plot(x,w1);hold on;plot(x,w2,'--');hold on;plot(x,w3,'-.');hold on;plot(x,w4,':');hold on;
    legend('\omega_{01}','\omega_{12}','\omega_{02}','\omega_{02}/2')
    xlabel('Ec');ylabel('频率v/Hz');
%%
elseif strcmpi(Var,'alpha') 
          w1 = zeros(1,41);w2 = zeros(1,41);w3 = zeros(1,41);w4 = zeros(1,41);
  parfor i=1:41
      alpha = 0.4+(i-1)*0.0125;
       [EL,~] = CaFluxQubit(Ej,Ec,alpha,beta,kappa,sigma,Cc,Cr,Csh,FluxBias,nk,nl,nm,nlevels);
       w1(i)=EL(3)-EL(1);w2(i)=EL(5)-EL(3);w3(i)=EL(5)-EL(1);w4(i)=(EL(5)-EL(1))/2;
%         w1(i)=(EL(3)-EL(1))/(2*pi);w2(i)=(EL(5)-EL(3))/(2*pi);w3(i)=(EL(5)-EL(1))/(2*pi);w4(i)=(EL(5)-EL(1))/(2*pi)/2;
  end
%    delete(gcp('nocreate'));%关闭Pool
   x=0.4:0.0125:0.9;
   figure;
   plot(x,w1);hold on;plot(x,w2,'--');hold on;plot(x,w3,'-.');hold on;plot(x,w4,':');hold on;
   legend('\omega_{01}','\omega_{12}','\omega_{02}','\omega_{02}/2')
   xlabel('alpha');ylabel('频率v/Hz');
%%
elseif strcmpi(Var,'beta') 
          w1 = zeros(1,41);w2 = zeros(1,41);w3 = zeros(1,41);w4 = zeros(1,41);
  parfor i=1:41
      beta = 1e-6+(i-1)*1e-5;
       [EL,~] = CaFluxQubit(Ej,Ec,alpha,beta,kappa,sigma,Cc,Cr,Csh,FluxBias,nk,nl,nm,nlevels);
       w1(i)=EL(3)-EL(1);w2(i)=EL(5)-EL(3);w3(i)=EL(5)-EL(1);w4(i)=(EL(5)-EL(1))/2;
%         w1(i)=(EL(3)-EL(1))/(2*pi);w2(i)=(EL(5)-EL(3))/(2*pi);w3(i)=(EL(5)-EL(1))/(2*pi);w4(i)=(EL(5)-EL(1))/(2*pi)/2;
  end
%    delete(gcp('nocreate'));%关闭Pool
   x=1e-6:1e-5:4.01e-4;
   figure;
   plot(x,w1);hold on;plot(x,w2,'--');hold on;plot(x,w3,'-.');hold on;plot(x,w4,':');hold on;
   legend('\omega_{01}','\omega_{12}','\omega_{02}','\omega_{02}/2')
   xlabel('beta');ylabel('频率v/Hz');
%%   
elseif strcmpi(Var,'Cc') 
  w1 = zeros(1,41);w2 = zeros(1,41);w3 = zeros(1,41);w4 = zeros(1,41);
  parfor i=1:41
      Cc = (1.19+(i-1)*0.4)*10^(-15);
       [EL,~] = CaFluxQubit(Ej,Ec,alpha,beta,kappa,sigma,Cc,Cr,Csh,FluxBias,nk,nl,nm,nlevels);
       w1(i)=EL(3)-EL(1);w2(i)=EL(5)-EL(3);w3(i)=EL(5)-EL(1);w4(i)=(EL(5)-EL(1))/2;
%         w1(i)=(EL(3)-EL(1))/(2*pi);w2(i)=(EL(5)-EL(3))/(2*pi);w3(i)=(EL(5)-EL(1))/(2*pi);w4(i)=(EL(5)-EL(1))/(2*pi)/2;
  end
%    delete(gcp('nocreate'));%关闭Pool
   x=1.19:0.4:17.19;
   figure;
   plot(x,w1);hold on;plot(x,w2,'--');hold on;plot(x,w3,'-.');hold on;plot(x,w4,':');hold on;
   legend('\omega_{01}','\omega_{12}','\omega_{02}','\omega_{02}/2')
   xlabel('Cc');ylabel('频率v/Hz');   


%%    
elseif strcmpi(Var,'Csh') 
  w1 = zeros(1,41);w2 = zeros(1,41);w3 = zeros(1,41);w4 = zeros(1,41);
  parfor i=1:41
      Csh = (i-1)*1.25*10^(-15);
       [EL,~] = CaFluxQubit(Ej,Ec,alpha,beta,kappa,sigma,Cc,Cr,Csh,FluxBias,nk,nl,nm,nlevels);
       w1(i)=EL(3)-EL(1);w2(i)=EL(5)-EL(3);w3(i)=EL(5)-EL(1);w4(i)=(EL(5)-EL(1))/2;
%         w1(i)=(EL(3)-EL(1))/(2*pi);w2(i)=(EL(5)-EL(3))/(2*pi);w3(i)=(EL(5)-EL(1))/(2*pi);w4(i)=(EL(5)-EL(1))/(2*pi)/2;
  end
%    delete(gcp('nocreate'));%关闭Pool
   x=0:1.25:50;
   figure;
   plot(x,w1);hold on;plot(x,w2,'--');hold on;plot(x,w3,'-.');hold on;plot(x,w4,':');hold on;
   legend('\omega_{01}','\omega_{12}','\omega_{02}','\omega_{02}/2')
   xlabel('Csh');ylabel('频率v/Hz'); 
%%
end
toc;    
end