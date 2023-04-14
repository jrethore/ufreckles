function [xc,yc,conn,elt,zc]=CleanMesh(xi,yi,conni,elti,zi)

if nargin<5, zi=0*xi;end

conn=conni;
elt=elti;
ll=max(max(max(xi)-min(xi),max(yi)-min(yi)),max(zi)-min(zi));
xc=[];yc=[];zc=[];
for ie=1:length(elt)
   inods=conn(ie,:);
   for in=1:elt(ie)
      xn=xi(inods(in));
      yn=yi(inods(in));
      zn=zi(inods(in));
      found=find(sqrt((xc-xn).^2+(yc-yn).^2+(zc-zn).^2)<1.e-9*ll);
      if isempty(found)||isempty(xc)
          xc=[xc;xn];
          yc=[yc;yn];
          zc=[zc;zn];
          conn(ie,in)=length(xc);
      else
          conn(ie,in)=found;
      end
   end
end

end