function [phi,dphix,dphiy]=CreateNURBSBasis2D(mesh_file,p,yi,xi,type_nurbs,isclosed,nband,sizeim)

if nargin < 5 , type_nurbs=1;end
if nargin < 6 , isclosed=[0,0];end
if nargin < 7 , nband=[];end
load(mesh_file,'xo','yo','uo','vo');
if ~isempty(nband)
    assert(nargin==8);
    yi=yi(nband);
    xi=xi(nband);
end
if numel(p)==1
    p=[p,p];
end
npt=length(yi(:));
indp=zeros((p(1)+1)*npt,1);
indn=zeros((p(1)+1)*npt,1);
val=zeros((p(1)+1)*npt,1);
dval=zeros((p(1)+1)*npt,1);
if type_nurbs
Nny=length(yo)+p(1)-1;
else
Nny=length(yo)*p(1)-(p(1)-1);    
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
	if type_nurbs
    [N]=NURBSBasisFunc(iy+p(1),p(1),yp',vo,2);
    else
    toto=toto+p;
    [N]=NURBSBasisFunc(toto,p(1),yp',vo,2);
    end

	Sel=length(yp);
    for ip=1:(p(1)+1)
        if type_nurbs
        indn(nel+(1:Sel))=iy+ip-1;
        else
        indn(nel+(1:Sel))=iy*p(1)-p+ip;
        end            
        indp(nel+(1:Sel))=found;
		val(nel+(1:Sel))=N(:,ip,1);
		dval(nel+(1:Sel))=N(:,ip,2);
        nel=nel+Sel;
    end
end
found=find(indp);
indp=indp(found);
indn=indn(found);
val=val(found);
dval=dval(found);
if isclosed(1)
    for ip=1:p(1)
    found=find(indn==Nny-ip+1);
    indn(found)=p(1)-ip+1;
    end
    Nny=Nny-p(1);
end
if isempty(nband)
    size1=npt;
else
    indp=nband(indp);
    size1=prod(sizeim);
end
phiy=sparse(indp,indn,val,size1,Nny);
dphiy1=sparse(indp,indn,dval,size1,Nny);


npt=length(xi(:));
indp=zeros((p(2)+1)*npt,1);
indn=zeros((p(2)+1)*npt,1);
val=zeros((p(2)+1)*npt,1);
dval=zeros((p(2)+1)*npt,1);
if type_nurbs
Nny=length(xo)+p(2)-1;
else
Nny=length(xo)*p(2)-(p(2)-1);    
end
nel=0;
toto=1;
Nelems=length(xo)-1;
for iy=1:Nelems
    if iy==Nelems
    found=find((xi(:)>=xo(iy))&(xi(:)<=xo(iy+1)));
    else
    found=find((xi(:)>=xo(iy))&(xi(:)<xo(iy+1)));
    end
    xp=xi(found);
	if type_nurbs
    [N]=NURBSBasisFunc(iy+p(2),p(2),xp',uo,2);
    else
    toto=toto+p;
    [N]=NURBSBasisFunc(toto,p(2),xp',uo,2);
    end

	Sel=length(xp);
    for ip=1:(p(2)+1)
        if type_nurbs
        indn(nel+(1:Sel))=iy+ip-1;
        else
        indn(nel+(1:Sel))=iy*p(2)-p(2)+ip;
        end            
        indp(nel+(1:Sel))=found;
		val(nel+(1:Sel))=N(:,ip,1);
		dval(nel+(1:Sel))=N(:,ip,2);
        nel=nel+Sel;
    end
end
indp((nel+1):end)=[];
indn((nel+1):end)=[];
val((nel+1):end)=[];
dval((nel+1):end)=[];


if isclosed(2)
    for ip=1:p(2)
    found=find(indn==Nny-ip+1);
    indn(found)=p(2)-ip+1;
    end
    Nny=Nny-p(2);
end
if isempty(nband)
    size1=npt;
else
    indp=nband(indp);
    size1=prod(sizeim);
end
phix=sparse(indp,indn,val,size1,Nny);
dphix1=sparse(indp,indn,dval,size1,Nny);
[indj,indi]=meshgrid(1:size(phix,2),1:size(phiy,2));

phi=phix(:,indj).*phiy(:,indi);
dphix=dphix1(:,indj).*phiy(:,indi);
dphiy=phix(:,indj).*dphiy1(:,indi);

% phi=phix(:)*(phiy(:)');
% 
% [indx,indy,val]=find(phi);
% indp=mod(indx,sizeim(1))+sizeim(1)*(mod(indy,sizeim(2))-1)+sizeim(1)*sizeim(2)*(mod(indy,sizeim(2))-1);
% indn=ceil(indx/(sizeim(1)))+Nnx*(ceil(indy/(sizeim(2)))-1);
% phi=sparse(indp,indn,val,prod(sizeim),Nnx*Nny);



end