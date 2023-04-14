function operation=GetGluingParameters(nmod,xo,yo,zo,conno,elt)
load(fullfile('TMP','params'),'param');
roi=param.roi;
if nargin<4,zo=0*xo;end
if nargin<5
    conn=[];
else
%     if size(conn,2)==3
%         conn=[conn,zeros(numel(elt),1)];
%     end
% conn=conn(:,[1,2,2,3,3,4,4,1]);
% conn=reshape(conn',2,4*numel(elt));
% ind=find(conn==0);
% ind=ind2sub(size(conn),ind);
% conn(ind,:)=[];
    
%     conn=[];
%     for ie=1:numel(elt)
%         switch elt(ie)
%             case 4
%                 conn(end+1,1:3)=conno(ie,[1,2,3]);
%                 conn(end+1,1:3)=conno(ie,[3,4,1]);
%             case 3
%                 conn(end+1,1:3)=conno(ie,[1,2,3]);
%         end
%         
%     end
        seg3=reshape(conn(elt==3,[1,2,1,3,2,3])',2,3*sum(elt==3))';
    seg4=reshape(conn(elt==4,[1,2,2,3,3,4,4,1])',2,4*sum(elt==4))';
    
    conn=unique(sort([seg3;seg4],2),'rows');

end

%default parameters
load(fullfile('TMP',sprintf('sample0')),'im0');
load(fullfile('TMP',sprintf('sample0_%d',1-1)),'sizeim');
mxo=mean(xo);
myo=mean(yo);
lx=max(xo)-min(xo);
ly=max(yo)-min(yo);
operation{1,:}={'translate',[-mxo;-myo;0]};
% if lx>ly
%     if sizeim(1)>sizeim(2)
%         operation{2}={'rotate',{[0;0;1]},{mxo;myo;0},0};
%         operation{3}={'scale',{mxo;myo;0},mean([diff(roi(1:2))/lx,diff(roi(3:4))/ly};
%     else
%         operation{2}={'rotate',{[0;0;1]},{mxo;myo;0},pi/2};
%         operation{3}={'scale',{mxo;myo;0},mean([diff(roi(1:2))/ly,diff(roi(3:4))/lx};
%     end
%     
% else
%     if sizeim(1)<sizeim(2)
%         operation{2}={'rotate',{[0;0;1]},{mxo;myo;0},0};
%         operation{3}={'scale',{mxo;myo;0},mean([diff(roi(1:2))/lx,diff(roi(3:4))/ly};
%     else
%         operation{2}={'rotate',{[0;0;1]},{mxo;myo;0},pi/2};
%         operation{3}={'scale',{mxo;myo;0},mean([diff(roi(1:2))/ly,diff(roi(3:4))/lx};
%     end
% end
%operation{4}={'translate',{roi(1)+0.5*sizeim(1);roi(3)+0.5*sizeim(2);0}};
tx=roi(1)+0.5*sizeim(1);
ty=roi(3)+0.5*sizeim(2);
scale=0.5*(diff(roi(1:2))/lx+diff(roi(3:4))/ly);
scaleo=scale;
angl=0;
operation{2,:}={'scale',[0;0;0],scale};
operation{3,:}={'rotate',[0;0;1],[0;0;0],angl};
operation{4,:}={'translate',[tx;ty;0]};

ff=figure;
hax = axes('Units','normalized');
colormap(gray)

imagesc(im0','Parent',hax,'CDataMapping','scaled')
hold on;
axis equal
axis xy;
axis off
hmesh=0;
plotmesh(operation);

% uicontrol('Style', 'edit','String',num2str(tx),...
%     'Position', [182 14 120 20],...
%     'Callback',  @setTx);
uicontrol('Style', 'slider',...
    'Min',tx-20,'Max',tx+20,'Value',tx,...
    'Position', [82 74 20 20],...
    'Callback', @setTx);
uicontrol('Style','text',...
    'Position',[22 74 40 20],...
    'String','Tx');
uicontrol('Style', 'slider',...
    'Min',ty-20,'Max',ty+20,'Value',ty,...
    'Position', [82 44 20 20],...
    'Callback', @setTy);
uicontrol('Style','text',...
    'Position',[22 44 40 20],...
    'String','Ty');
uicontrol('Style', 'slider',...
    'Min',log10(0.1*scaleo),'Max',log10(10*scaleo),'Value',log10(scale),...
    'Position', [82 104 20 20],...
    'Callback', @setScale);
uicontrol('Style','text',...
    'Position',[22 104 40 20],...
    'String','Scale');
uicontrol('Style', 'slider',...
    'Min',angl*180/pi-10,'Max',angl*180/pi+10,'Value',angl*180/pi,...
    'Position', [82 134 20 20],...
    'Callback', @setAngl);
uicontrol('Style','text',...
    'Position',[22 134 40 20],...
    'String','Angle');

uicontrol('Style', 'pushbutton', 'String', 'OK',...
    'Position', [40 14 50 20],...
    'Callback', @finish);
uiwait
    function finish(hObj,event)
        close(ff)
        pause(0.1)
        uiresume
        
    end
    function setTy(hObj,event)
        ty = (get(hObj,'Value'));
        set(hObj,'Min',ty-20,'Max',ty+20,'Value',ty);
        operation{4,:}={'translate',[tx;ty;0]};
        delete(hmesh)
plotmesh(operation)
        
    end
    function setTx(hObj,event)
        tx = (get(hObj,'Value'));
        set(hObj,'Min',tx-20,'Max',tx+20,'Value',tx);
        operation{4,:}={'translate',[tx;ty;0]};
        delete(hmesh)
plotmesh(operation)
        
    end
    function setScale(hObj,event)
        scale = 10^(get(hObj,'Value'));
operation{2,:}={'scale',[0;0;0],scale};
        delete(hmesh)
plotmesh(operation)
    end
    function setAngl(hObj,event)
        angl = (get(hObj,'Value'))*pi/180;
        set(hObj,'Min',angl*180/pi-10,'Max',angl*180/pi+10,'Value',angl*180/pi);
operation{3,:}={'rotate',[0;0;1],[0;0;0],angl};
        delete(hmesh)
plotmesh(operation)
        
    end
    function plotmesh(ops)
        [xi,yi,zi]=GlueMesh(ops,xo,yo,zo);
        if isempty(conn)
        hmesh=plot(xi,yi,'b+','LineWidth',0.5,'MarkerSize',10);
        else
        hmesh=plot(xi(conn)',yi(conn)','b+','LineWidth',0.5);
        end
%                 hmesh=trimesh(conn(:,1:3),xi,yi,zi,zi);
%    set(hmesh,'EdgeColor','r',...
%         'FaceColor','none',...
%         'Marker','none');

    end
end