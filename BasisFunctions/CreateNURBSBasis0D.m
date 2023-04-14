function [phi,dphi,ddphi]=CreateNURBSBasis0D(mesh_file,p,yi,type_nurbs,isclosed,nband,sizeim)

if nargin < 4 , type_nurbs=1;end
if nargin < 5 , isclosed=0;end
if nargin < 6 , nband=[];end
load(mesh_file,'yo','vo');
if ~isempty(nband)
    assert(nargin==7);
    yi=yi(nband);
end
npt=length(yi(:));
indp=zeros((p+1)*npt,1);
indn=zeros((p+1)*npt,1);
val=zeros((p+1)*npt,1);
dval=zeros((p+1)*npt,1);
ddval=zeros((p+1)*npt,1);
if type_nurbs
Nny=length(yo)+p-1;
else
Nny=length(yo)*p-(p-1);    
end
nel=0;
toto=1;
Nelems=length(yo)-1;

for iy=1:Nelems
    if iy==Nelems
    found=find((yi(:)>=yo(iy))&(yi(:)<=yo(iy+1)));
    else
    found=find((yi(:)>=yo(iy))&(yi(:)<yo(iy+1)));
    end
    yp=yi(found);
    if isempty(yp)
        figure
plot(sort(yi),sort(yi),'b+');
hold on
plot(sort(yo),sort(yo),'ro')
    yo(iy)
    yo(iy+1)
    pause
    end
	if type_nurbs
    [N]=NURBSBasisFunc(iy+p,p,yp',vo,2);
    else
    toto=toto+p;
    [N]=NURBSBasisFunc(toto,p,yp',vo,2);
    end

	Sel=length(yp);
    for ip=1:(p+1)
        if type_nurbs
        indn(nel+(1:Sel))=iy+ip-1;
        else
        indn(nel+(1:Sel))=iy*p-p+ip;
        end            
        indp(nel+(1:Sel))=found;
		val(nel+(1:Sel))=N(:,ip,1);
		dval(nel+(1:Sel))=N(:,ip,2);
		ddval(nel+(1:Sel))=N(:,ip,3);
        nel=nel+Sel;
    end
end
found=find(indp);
indp=indp(found);
indn=indn(found);
val=val(found);
dval=dval(found);
ddval=ddval(found);
if isclosed
    for ip=1:p
    found=find(indn==Nny-ip+1);
    indn(found)=p-ip+1;
    end
    Nny=Nny-p;
end
if isempty(nband)
    size1=npt;
else
    indp=nband(indp);
    size1=prod(sizeim);
end
phi=sparse(indp,indn,val,size1,Nny);
dphi=sparse(indp,indn,dval,size1,Nny);
ddphi=sparse(indp,indn,ddval,size1,Nny);

end