function [phi]=CreateElementwiseBasis(mesh_file,sizeim,pscale,selected_elts,GaussPts,full_size)

load(mesh_file,'xo','yo','Nnodes','Nelems','Smesh','mesh_size');
if nargin < 3 , pscale=1;end
if nargin < 4 , selected_elts=[];end
if nargin < 5 , GaussPts='pixels';end
if nargin < 6 , full_size=false;end

if isempty(selected_elts)||full_size
    ny=prod(Nelems);
else
    ny=length(selected_elts);
end

xo=(xo-0.5)*pscale+0.5;
yo=(yo-0.5)*pscale+0.5;
    xo=reshape(xo,Nnodes)+0.5;
    yo=reshape(yo,Nnodes)+0.5;
    xo=xo(1:Nnodes(1),1)';
    yo=yo(1,1:Nnodes(2));

ince(1)=1;ince(2)=ince(1)*Nelems(1);
incp(1)=1;incp(2)=incp(1)*sizeim(1)*pscale;
incg(1)=4;incg(2)=incg(1)*Nelems(1);

 if strcmp(GaussPts,'pixels') , npts=prod(Smesh)*pscale^2;nx=prod(sizeim*pscale);
 elseif strcmp(GaussPts,'Gauss_points'), npts=4*prod(Nelems);nx=4*prod(Nelems);
 end

indn=zeros(npts,1);
indp=zeros(npts,1);
val=zeros(npts,1);
nel=0;

    
for i1=1:Nelems(1)
      for i2=1:Nelems(2)

            ielt=1+ince(1)*(i1-1)+ince(2)*(i2-1);
            ienr=find(selected_elts==ielt);
            if (~isempty(ienr))||(isempty(selected_elts))
                  
                  if (~isempty(ienr))&&(~full_size)
                        ielt=ienr;
                  end
                  indxn=i1+(0:1);
                  indyn=i2+(0:1);

            xn=(xo(indxn)-xo(i1));
            yn=(yo(indyn)-yo(i2));

            xmax=max(abs(xn(:)));
            ymax=max(abs(yn(:)));

            xn=xn/xmax;
            yn=yn/ymax;
            if strcmp(GaussPts,'pixels')
                xp=xo(indxn(1)):(xo(indxn(length(indxn)))-1);
                yp=yo(indyn(1)):(yo(indyn(length(indyn)))-1);
                [Yp,Xp]=meshgrid(yp,xp);
                ipix=1+incp(1)*(Xp-1)+incp(2)*(Yp-1);
            elseif strcmp(GaussPts,'Gauss_points')
                dg=1/(2*sqrt(3));
                xp=sort([xn(1:length(xn)-1)+0.5-dg;xn(2:length(xn))-0.5+dg]);
                yp=sort([yn(1:length(yn)-1)+0.5-dg;yn(2:length(yn))-0.5+dg]);
                [Yp,Xp]=meshgrid(yp,xp);
                groot=1+incg(1)*(indxn(1)-1)+incg(2)*(indyn(1)-1);
                ipix(1:2)=groot+(1:2)-1;
                if length(indxn)==3
                ipix(2+(1:2))=ipix(1:2)+incg(1);
                end
                ipix(2*(length(indxn)-1)+(1:2*(length(indxn)-1)))=ipix(1:2*(length(indxn)-1))+2;
                if length(indyn)==3
                ipix(4*(length(indxn)-1)+(1:4*(length(indxn)-1)))=ipix(1:4*(length(indxn)-1))+incg(2);
                end
            end

                  Sel=length(ipix(:));
                  indn(nel+(1:Sel))=ielt;
                  indp(nel+(1:Sel))=ipix(:);
                  val(nel+(1:Sel))=1;
                  nel=nel+Sel;
            clear ipix

            end
      end
end    
    
%end


if nel<length(indn)
    indn=indn(1:nel);
    indp=indp(1:nel);
    val=val(1:nel);
end



phi=sparse(indp,indn,val,nx,ny);

end
