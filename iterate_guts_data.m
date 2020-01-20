function [Kgax,Ktgx,Ksgx,Kgtx,Cxg,Cxdatag,relerr,nit,crit,relerrp] = ...
    iterate_guts_data(glnodes,glweights,t,alpha,Ca,Cdatag,r,...
    Kgax,Ktgx,Ksgx,Kgtx,Cxg,Cxdatag,relerr,nit,crit)

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
Ax=[[-(Ktgx+Ksgx);Ktgx],[Kgtx;-Kgtx]];
x=[Kgax;Ktgx;Ksgx;Kgtx];

M=mat_derivative(alpha,Ax,Kgax,Ca,t,Cxg,glnodes,glweights);
M=[M(:,1),M(:,3)-M(:,2),-M(:,2),M(:,4)-M(:,5)];
h=(r*diag([1,1,1,1])+M.'*M)\(M.'*(Cdatag-Cxdatag));
xph=x+h;
if any(xph<=0)
    h(xph<=0)=0;
end
x=x+h;

% keep trace of previous values
[Kgaxp,Ktgxp,Ksgxp,Kgtxp,relerrp]=deal(Kgax,Ktgx,Ksgx,Kgtx,relerr);

Kgax=x(1);Ktgx=x(2);Ksgx=x(3);Kgtx=x(4);
Ax=[[-(Ktgx+Ksgx);Ktgx],[Kgtx;-Kgtx]];
Cxg=concentration_K1(Ax,Ca,0,[0;0],t,glnodes,glweights);
Cxdatag=Kgax*(alpha*Cxg).';
relerr=norm(Cxdatag-Cdatag)/norm(Cdatag);

crit(1,:)=[];
crit(end+1,:)=[abs(Kgax-Kgaxp)/abs(Kgax),...
    abs(Ktgx-Ktgxp)/abs(Ktgx),...
    abs(Ksgx-Ksgxp)/abs(Ksgx),...
    abs(Kgtx-Kgtxp)/abs(Kgtx),...
    abs(relerr-relerrp)/abs(relerr)];

end