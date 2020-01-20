function [Kfbx,Kmfx,Kbfx,Kfmx,Cx,Cxdata,relerr,nit,crit,relerrp]=...
    iterate_tumor_data(glnodes,glweights,t,alpha,Ca,Cdata,r,...
    Kfbx,Kmfx,Kbfx,Kfmx,Cx,Cxdata,relerr,nit,crit)

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
Ax=[[-(Kmfx+Kbfx);Kmfx],[Kfmx;-Kfmx]];
x=[Kfbx;Kmfx;Kbfx;Kfmx];

M=mat_derivative(alpha,Ax,Kfbx,Ca,t,Cx,glnodes,glweights);
M=[M(:,1),M(:,3)-M(:,2),-M(:,2),M(:,4)-M(:,5)];

h=(r*diag([1,1,1,1])+M.'*M)\(M.'*(Cdata-Cxdata));

xph=x+h;
if any(xph<=0)
    h(xph<=0)=0;
end
x=x+h;

% keep trace of previous values
[Kfbxp,Kmfxp,Kbfxp,Kfmxp,relerrp]=deal(Kfbx,Kmfx,Kbfx,Kfmx,relerr);

Kfbx=x(1);Kmfx=x(2);Kbfx=x(3);Kfmx=x(4);
Ax=[[-(Kmfx+Kbfx);Kmfx],[Kfmx;-Kfmx]];
Cx=concentration_K1(Ax,Ca,0,[0;0],t,glnodes,glweights);
Cxdata=Kfbx*(alpha*Cx).';
relerr=norm(Cxdata-Cdata)/norm(Cdata);

crit(1,:)=[];
crit(end+1,:)=[abs(Kfbx-Kfbxp)/abs(Kfbx),...
    abs(Kmfx-Kmfxp)/abs(Kmfx),...
    abs(Kbfx-Kbfxp)/abs(Kbfx),...
    abs(Kfmx-Kfmxp)/abs(Kfmx),...
    abs(relerr-relerrp)/abs(relerr)];

end