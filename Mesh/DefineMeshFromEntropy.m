function smin=DefineMeshFromEntropy(meshfile,roi,im00,smin)
if nargin<4,smin=0;end
im0=im00(roi(1):roi(2),roi(3):roi(4));
load(meshfile,'xo','yo','Nnodes','Nelems','conn','elt','selected')
Xe=zeros(prod(Nelems),1);
Ye=zeros(prod(Nelems),1);
for ie=1:prod(Nelems)
   Xe(ie)=mean(xo(conn(ie,1:elt(ie))))+roi(1)-1; 
   Ye(ie)=mean(yo(conn(ie,1:elt(ie))))+roi(3)-1; 
end
% Xe=mean(xo(conn),2)+roi(1)-1;
% Ye=mean(yo(conn),2)+roi(3)-1;
ff=figure;
hax = axes('Units','pixels');
colormap(gray)

imagesc(im00')
hold on;
axis equal
axis xy;
axis off
[S]=GetEntropy(meshfile,im0);
%plot(xo,yo,'r.')
S=(S-min(S))/(max(S)-min(S))*255;
hmesh=0;
plotmesh(S,smin);

title('Element selection from image entropie')
uicontrol('Style', 'slider',...
    'Min',0,'Max',255,'Value',smin,...
    'Position', [182 14 120 20],...
    'Callback', @setSmin);
uicontrol('Style', 'pushbutton', 'String', 'OK',...
    'Position', [373 13 50 20],...
    'Callback', @savemesh);
uiwait

    function savemesh(hObj,event)
        
        ekeep=find(S>=smin);
        
        
        
        conn=conn(ekeep,:);
        elt=elt(ekeep,:);
        [keep,conn]=GetNodesFromElts(conn,elt);
        xo=xo(keep);yo=yo(keep);
        Nnodes=[length(xo),1,1];
        Nelems=[size(conn,1),1,1];
        selected=selected(keep);
        save(meshfile,'xo','yo','conn','elt','Nnodes','Nelems','-append');
        close(ff)
        pause(0.1)
        uiresume
    end
    function setSmin(hObj,event)
        
        smin = get(hObj,'Value');
        delete(hmesh);
        plotmesh(S,smin);
        
    end
    function plotmesh(F,fmin)
        
        hmesh=scatter((Xe)',(Ye)',10+0*Ye',255*(F>=fmin),'filled');
        
    end
end