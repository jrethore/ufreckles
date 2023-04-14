function getSurf(meshfile, xo, yo, zo, conn, elt, select,keepall)
if nargin<8,keepall=1;end

conn2D = [];
elt2D = [];
oldtonew=[];
facet=[1,2;2,3;3,1];
faceq=[4,3,2,1;...
    5,6,7,8;...
    1,2,6,5;...
    5,8,4,1;...
    3,4,8,7;...
    7,6,2,3];
for i1=1:size(conn,1)
    
    if elt(i1) == 4

%         inods=conn(i1,1:elt(i1));
% 
%         is = [sum(select==inods(1)), sum(select==inods(2)), sum(select==inods(3)), sum(select==inods(4))];
% 
%         if sum(is) == 3
% found=find(is);
%             for in=1:3
%                 inod=inods(found(in));
%                if ~any(oldtonew==inod)
%                    oldtonew=[oldtonew,inod];
%                end
%                 
%                 
%             end
%             
%             
%             conn2D(end+1,1:8) = [inods(find(is==1)) 0 0 0 0 0];
%             elt2D(end+1) = 3;
% 
%     %         inods2D=conn2D(end,1:elt2D(end));
%     %         
%     %         xn=xo(inods2D);
%     %         yn=yo(inods2D);
%     %         zn=zo(inods2D);
%     %         
%     %         hold on;
%     %         fill3(xn([1 2 3 1]),yn([1 2 3 1]),zn([1 2 3 1]),'y')
%     %         hold off;
%     %         pause(0.001)
% 
%         end
    error('NOT CODED YET');
    elseif elt(i1) == 8
        
        inods=conn(i1,1:elt(i1));

        is = [sum(select==inods(1)), sum(select==inods(2)), sum(select==inods(3)), sum(select==inods(4)), sum(select==inods(5)), sum(select==inods(6)), sum(select==inods(7)), sum(select==inods(8))];

        if any(is>0)
            
            for iface=1:6
                onface=is(faceq(iface,:));
                if sum(onface)==4
                    conn2D(end+1,1:4) = [inods(faceq(iface,:))];
                    elt2D(end+1) = 4;
                    for in=1:4
                        inod=inods(faceq(iface,in));
                        if ~any(oldtonew==inod)
                            oldtonew=[oldtonew,inod];
                        end
                        
                    end
                end
                
            end
        end
%         if sum(is) == 4
% 
%  found=find(is);
%             for in=1:4
%                 inod=inods(found(in));
%                if ~any(oldtonew==inod)
%                    oldtonew=[oldtonew,inod];
%                end
%                 
%                 
%             end
%            conn2D(end+1,1:8) = [inods(find(is==1)) 0 0 0 0];
%             elt2D(end+1) = 4;
%             
%         end
        
    end
        
end

conn = conn2D;
elt = elt2D;
if ~keepall
    for ie=1:length(elt)
        for in=1:elt(ie)
            conn(ie,in)=find(select==conn(ie,in));
        end
    end
    xo=xo(select);
    yo=yo(select);
    zo=zo(select);
end

Xo = xo;
Yo = yo;
Zo = zo;


Nnodes = [length(Xo),1,1];
Nelems = [length(elt),1,1];
rint=0;
ng=0;
ns=ones(2,1);
Smesh=[];

selected=ones(length(xo),1);
if keepall
    selected(selected)=0;
end

save(meshfile,'rint','Xo','Yo','Zo','xo','yo','zo','Nnodes','Nelems','Smesh','conn','ng','ns','elt','selected');

end