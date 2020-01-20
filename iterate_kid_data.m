function [Kfax,Kmax,Kafx,Kmfx,Kfmx,Ktmx,Kutx,Cx1,Cx2,Cxdata,relerr,nit,crit,relerrp]=...
    iterate_kid_data(glnodes,glweights,t,alpha,Ca,Cdata,r,...
    Kfax,Kmax,Kafx,Kmfx,Kfmx,Ktmx,Kutx,Cx1,Cx2,Cxdata,relerr,nit,crit)

% crit is a 8x8 matrix, each row for the relative errors between successive
% iterations, ie for example, if nit=8:
% first row: relative errors between step 3 and step 4
% second row: relative errors between step 4 and step 5
% third row: relative errors between step 5 and step 6
% fourth row: relative errors between step 6 and step 7
% fifth row: relative errors between step 7 and step 8
% ...to eighth row
% then, after the function, nit=9 and crit becomes
% first row: relative errors between step 4 and step 5
% second row: relative errors between step 5 and step 6
% third row: relative errors between step 6 and step 7
% fourth row: relative errors between step 7 and step 8
% fifth row: relative errors between step 8 and step 9
% ...to eighth row

nit=nit+1;

Ax=[[-(Kafx+Kmfx);Kmfx;0],[Kfmx;-(Kfmx+Ktmx);Ktmx],[0;0;-Kutx]];
x=[Kfax;Kmax;Kafx;Kmfx;Kfmx;Ktmx;Kutx];

M1=mat_derivative_kid(alpha,Ax,Kfax,Ca,[1;0;0],t,Cx1,glnodes,glweights); % M1=[dalphaCdK1,dalphaCdA] on e1
M2=mat_derivative_kid(alpha,Ax,Kmax,Ca,[0;1;0],t,Cx2,glnodes,glweights); % M2=[dalphaCdK2,dalphaCdA] on e2
M=M1+M2;
M=[M1(:,1),M2(:,1),-M(:,2),M(:,3)-M(:,2),M(:,5)-M(:,6),M(:,7)-M(:,6),-M(:,10)];

h=(r*diag([1,1,1,1,1,1,1])+M.'*M)\(M.'*(Cdata-Cxdata));

xph=x+h;
if any(xph<=0)
    h(xph<=0)=0;
end
x=x+h;

% keep trace of previous values
[Kfaxp,Kmaxp,Kafxp,Kmfxp,Kfmxp,Ktmxp,Kutxp,relerrp]=deal(Kfax,Kmax,Kafx,Kmfx,Kfmx,Ktmx,Kutx,relerr);

Kfax=x(1);Kmax=x(2);Kafx=x(3);Kmfx=x(4);Kfmx=x(5);Ktmx=x(6);Kutx=x(7);

Ax=[[-(Kafx+Kmfx);Kmfx;0],[Kfmx;-(Kfmx+Ktmx);Ktmx],[0;0;-Kutx]];
Cx1=concentration_K1_kid(Ax,Ca,[1;0;0],0,[0;0;0],t,glnodes,glweights); %[Ca;0;0] on e1
Cx2=concentration_K1_kid(Ax,Ca,[0;1;0],0,[0;0;0],t,glnodes,glweights); %[0;Ca;0] on e2
Cx=Kfax*Cx1+Kmax*Cx2;
Cxdata=(alpha*Cx).';

relerr=norm(Cxdata-Cdata)/norm(Cdata);

crit(1,:)=[];
crit(end+1,:)=[abs(Kfax-Kfaxp)/abs(Kfax),...
    abs(Kmax-Kmaxp)/abs(Kmax),...
    abs(Kafx-Kafxp)/abs(Kafx),...
    abs(Kmfx-Kmfxp)/abs(Kmfx),... 
    abs(Kfmx-Kfmxp)/abs(Kfmx),...
    abs(Ktmx-Ktmxp)/abs(Ktmx),...
    abs(Kutx-Kutxp)/abs(Kutx),...
    abs(relerr-relerrp)/abs(relerr)];

end