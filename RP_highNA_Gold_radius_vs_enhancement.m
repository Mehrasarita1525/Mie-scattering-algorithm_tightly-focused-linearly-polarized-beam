clc
close all
clear all
 
tic
%%%%%%%%%%%%%%%%%%%%% Gold NP %%%%%%%%%%%%%%%%%%%%%%%%%%
n_m=1.52;                   %optical index of medium(glass/air)
 
N=15;                     %order of multipole


alpha=pi/3;

for n=1
    syms t;
    at=cos(t)^(1/2);
    Pn=legendreP(n,cos(t));
    dP=diff(Pn,t); 
    PL=at.*dP.*sin(t);
    Palpha=int(PL,t,0, alpha);
    B1=(-1i)^n.*Palpha.*3./4;
    B1_radial=abs(B1)^2;
    B1_radial=double(B1_radial);
    
end
    C1=0;
for a=0:1:100
    C1=C1+1;
    C2=0;
for Inter_lambda=300:1:900
    C2=C2+1;
   
    lambda=[301 311 320 331 342 354 368 381 397 413 430 451 471 496 521 548 582 616 659 704 755 821 891 900];    %lambda=wavelength of incident light(nm)  
    R=[1.53 1.53 1.54 1.48 1.48 1.50 1.48 1.46 1.47 1.46 1.45 1.38 1.31 1.04 0.62 0.43 0.29 0.21 0.14 0.13 0.14 0.16 0.17 0.18]; 
    I=[1.889 1.893 1.898 1.883 1.871 1.866 1.895 1.933 1.952 1.958 1.948 1.914 1.849 1.833 2.081 2.455 2.863 3.272 3.697 4.103 4.542 5.083 5.663 5.720];
    
    Aur=interp1(lambda,R,Inter_lambda);
    Aui=interp1(lambda,I,Inter_lambda);
    
   Kr_R=0;
   Kt_R=0;

% Input Parameters  
    n_Au=Aur+(1i*Aui);          %optical index(complex) of particle(gold/silver)
   k=2*pi*n_m/Inter_lambda;          %wavenumber
   m=n_Au/n_m;                 %m=relative refractive index
   x=k*a;                      %x=size parameter
   z=m*x;
 
for n=1:N
   sx=sqrt(pi*x/2);
   sz=sqrt(pi*z/2);

   Phinx=sx.*besselj(n+0.5,x);                 
   Phin1x=sx.*besselj(n-0.5,x);   
   Phinz=sz.*besselj(n+0.5,z);                       
   Phin1z=sz.*besselj(n-0.5,z);
     
   Xinx=sx.*((besselj(n+0.5,x))+1i*(bessely(n+0.5,x)));       
   Xin1x=sx.*((besselj(n-0.5,x))+1i*(bessely(n-0.5,x))); 
      
   Phindx=Phin1x-(n.*Phinx./x);
   Phindz=Phin1z-(n.*Phinz./z);
   Xindx=Xin1x-(n.*(Xinx)./x);
   
an=((m*Phinz.*Phindx)-(Phinx.*Phindz))./((m*Phinz.*Xindx)-(Xinx.*Phindz));
bn=((Phinz.*Phindx)-(m*Phinx.*Phindz))./((Phinz.*Xindx)-(m*Xinx.*Phindz));

W=((abs(an))^2.*(abs(Xinx))^2)+(abs(Phinx)^2)+(2*real(conj(an).*conj(Xinx)).*Phinx);
X=((abs(an))^2.*(abs(Xindx))^2)+(abs(Phindx)^2)+(2*real(conj(an).*conj(Xindx)).*Phindx);
Y=((n^2)*(n+1)^2)./((2*n)+1);
Z=(n*(n+1))./((2*n)+1);

 
    
    syms t;
    at=cos(t)^(1/2);
    Pn=legendreP(n,cos(t));
    dP=diff(Pn,t); 
    PL=at.*dP.*sin(t);
    Palpha=int(PL,t,0, alpha);
    Bn=(-1i)^n.*Palpha.*(2*n+1)./(2*n*(n+1));
    Bn_radial=abs(Bn)^2;
    Bn_radial=double(Bn_radial);
        
    Bn_R=Bn_radial./B1_radial;

    KrR=Bn_R.*Y.*W;
    Kr_R=KrR+Kr_R;
   
end
    Kr_Radial=9.*(Kr_R)./(4.*x^4);
  
d(C1,C2)=Kr_Radial;
end
end

a=0:1:100;
Inter_lambda=300:100:900;
figure(1);imagesc(Inter_lambda,a,d);
colormap(hot)
xlabel('Wavelength(nm)')
ylabel('Radius of nanoparticle(nm)')
set(gca,'fontweight','bold','FontSize',12)
axis square

toc

 