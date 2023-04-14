function T=UTcalib(path0,fil0)

cd(path0);
[filc,pathres]=uiputfile({'*.cal','Ufreckles calibration file (*.cal)'},'Save as a new calibration data set');

fils=SortImageFiles(fil0);

        fo=0;
        for iim=(length(fils)-9):length(fils)
            tmp=double(readim(fils{iim}));
            fo=fo+tmp;
        end
        fo=fo/10;
T.fo=fo;

            scrsz = get(0,'ScreenSize');
fig=figure(100);
set(fig,'NumberTitle','off','Name','IR speckle calibration','Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2],'Renderer','painters');
ax=axes('Parent',fig,'drawmode','fast','FontSize',16);
gim=imagesc(fo');
colormap('gray')
hold on
axis equal
axis xy;
axis off


title('Select the calibration area')
sizeim=size(fo);
zone=[];
while isempty(zone)
    waitforbuttonpress
    if strcmp(get(fig,'selectiontype'),'normal')
        point1=get(gca,'CurrentPoint');
        point1=point1(1,1:2);
            rbbox;
            point2=get(gca,'CurrentPoint');
            point2=point2(1,1:2);
            ly=abs((point1(1)-point2(1)));
            lx=abs((point1(2)-point2(2)));
    end
            if ~(max(lx,ly)<1)
                zone=[sort(round([point1(1),point2(1)])),sort(round([point1(2),point2(2)]))];
                zone(1)=max(zone(1),1);
                zone(2)=min(zone(2),sizeim(1));
                zone(3)=max(zone(3),1);
                zone(4)=min(zone(4),sizeim(2));
            end
end
plot(zone([1,2,2,1]),zone([3,3,4,4]),'r-','LineWidth',2)
            pause(1.0)
            title('Running calibration...')

fo=fo(zone(1):zone(2),zone(3):zone(4));
to=mean(fo(:));
tts=zeros(size(fil0,2),1);
txs=zeros(size(fil0,2),1);
%tx=0;
for iim=1:(size(fil0,2))
    tmp=double(readim(fil0{iim}));
    tmp=tmp(zone(1):zone(2),zone(3):zone(4));
    tt=mean(tmp(:));
    tmp=((-to*(tmp-fo)+fo*(tt-to))./(tt-to-(tmp-fo)));
%    tx=tx+mean(tmp(:));
    txs(iim)=mean(tmp(:));
    tts(iim)=tt;
end
%tx=tx/size(fil0,2)
dtts=gradient(tts);
tx=sum(txs.*dtts)/sum(dtts);


T.to=to;
T.tx=tx;

 save(filc,'T');
            title('End running calibration...')
            pause(1.0)

end