function [Kfax,Kmfx,Ksfx,Kfmx,Cxla,Cxlv,Cxdatal,relerr,nit,crit,relerrp] = ...
    iterate_liver_data(glnodes,glweights,t,alpha,V,Ca,Ksgx,Cv,Cdatal,r,...
    Kfax,Kmfx,Ksfx,Kfmx,Cxla,Cxlv,Cxdatal,relerr,nit,crit)

% crit is a 5x5 matrix, each row for the relative errors between successive
% iterations, ie for example, if nit=8:
% first row: relative errors between step 3 and step 4
% second row: relative errors between step 4 and step 5
% third row: relative errors between step 5 and step 6
% fourth row: relative errors between step 6 and step 7
% fifth row: relative errors between step 7 and step 8
% then, after the function, nit=9 and crit becomes
% first row: relative errors between step 4 and step 5
% second row: relative errors between step 5 and step 6
% third row: relative errors between step 6 and step 7
% fourth row: relative errors between step 7 and step 8
% fifth row: relative errors between step 8 and step 9

nit=nit+1;
Ax=[[-(Kmfx+Ksfx);Kmfx],[Kfmx;-Kfmx]];
x=[Kfax;Kmfx;Ksfx;Kfmx];

M=(1-V)*(mat_derivative(alpha,Ax,Kfax,Ca,t,Cxla,glnodes,glweights)+...
    mat_derivative(alpha,Ax,Ksgx,Cv,t,Cxlv,glnodes,glweights,'A'));

M=[M(:,1),M(:,3)-M(:,2),-M(:,2),M(:,4)-M(:,5)];

h=(r*diag([1,1,1,1])+M.'*M)\(M.'*(Cdatal-Cxdatal));
xph=x+h;
if any(xph<=0)
    h(xph<=0)=0;
end
x=x+h;

% keep trace of previous values
[Kfaxp,Kmfxp,Ksfxp,Kfmxp,relerrp]=deal(Kfax,Kmfx,Ksfx,Kfmx,relerr);

Kfax=x(1);Kmfx=x(2);Ksfx=x(3);Kfmx=x(4);
Ax=[[-(Kmfx+Ksfx);Kmfx],[Kfmx;-Kfmx]];

Cxla=concentration_K1(Ax,Ca,0,[0;0],t,glnodes,glweights);
Cxlv=concentration_K1(Ax,Cv,0,[0;0],t,glnodes,glweights);

Cxdatal=(1-V)*(alpha*(Kfax*Cxla+Ksgx*Cxlv)).'+ V/100 * (15*Ca(t)+ 85*Cv(t)).';
relerr=norm(Cxdatal-Cdatal)/norm(Cdatal);

crit(1,:)=[];
crit(end+1,:)=[abs(Kfax-Kfaxp)/abs(Kfax),...
    abs(Kmfx-Kmfxp)/abs(Kmfx),...
    abs(Ksfx-Ksfxp)/abs(Ksfx),...
    abs(Kfmx-Kfmxp)/abs(Kfmx),...
    abs(relerr-relerrp)/abs(relerr)];

end