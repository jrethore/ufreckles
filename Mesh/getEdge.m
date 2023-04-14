function getEdge(meshfile, xo, yo, zo, conn, elt, select,keepall)
if nargin<8,keepall=1;end
conn1D = [];
elt1D = [];
oldtonew=[];
facet=[1,2;2,3;3,1];
faceq=[1,2;2,3;3,4;4,1];
for i1=1:size(conn,1)
    
    if elt(i1) == 3
        
        inods=conn(i1,1:elt(i1));
        
        is = [sum(select==inods(1)), sum(select==inods(2)), sum(select==inods(3))];
        
        if any(is>0)
            
            for iface=1:3
                onface=is(facet(iface,:));
                if sum(onface)==2
                    conn1D(end+1,1:4) = [inods(facet(iface,:)) 0 0];
                    elt1D(end+1) = 2;
                    for in=1:2
                        inod=inods(facet(iface,in));
                        if ~any(oldtonew==inod)
                            oldtonew=[oldtonew,inod];
                        end
                        
                    end
                end
                
            end
        end
        
        %         if sum(is) == 2
        %             found=find(is==1);
        %             if ~(diff(found)==1)
        %                 found=flipdim(found,2);
        %             end
        %             conn1D(end+1,1:4) = [inods(found) 0 0];
        %             elt1D(end+1) = 2;
        %
        %         end
        
    elseif elt(i1) == 4
        
        inods=conn(i1,1:elt(i1));
        
        is = [sum(select==inods(1)), sum(select==inods(2)), sum(select==inods(3)), sum(select==inods(4))];
        if any(is>0)
            
            for iface=1:4
                onface=is(faceq(iface,:));
                if sum(onface)==2
                    conn1D(end+1,1:4) = [inods(faceq(iface,:)) 0 0];
                    elt1D(end+1) = 2;
                    for in=1:2
                        inod=inods(faceq(iface,in));
                        if ~any(oldtonew==inod)
                            oldtonew=[oldtonew,inod];
                        end
                        
                    end
                end
                
            end
        end
        
%         if sum(is) == 2
%             
%             found=find(is);
%             conn1D(end+1,1:4) = [found 0 0];
%             elt1D(end+1) = 2;
%         end
        
    end
    
end

conn = conn1D;
elt = elt1D;
if ~keepall
    for ie=1:length(elt)
        for in=1:elt(ie)
            conn(ie,in)=find(odltonew==conn(ie,in));
        end
    end
    xo=xo(oldtonew);
    yo=yo(oldtonew);
    zo=zo(oldtonew);
end
Xo = xo;
Yo = yo;
Zo = zo;
Nnodes = [length(Xo),1,1];
Nelems = [length(elt),1,1];
rint=0;
ng=1;
ns=ones(2,1);
Smesh=[];

save(meshfile,'rint','Xo','Yo','Zo','xo','yo','zo','Nnodes','Nelems','Smesh','conn','ng','ns','elt');

end