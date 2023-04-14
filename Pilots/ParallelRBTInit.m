function [Ue,Ve,We]=ParallelRBTInit(roi,h,Nelems,conn,xo,yo,zo,im0,jm3,dec,sizeim)
        Ue=zeros(Nelems);
        Ve=zeros(Nelems);
        We=zeros(Nelems);
            hh=ceil(-h/2):floor(h/2);
    [Urbt]=rbt3(im0(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6)),jm3(roi(1):roi(2),roi(3):roi(4),roi(5):roi(6)));
Ug=round(Urbt(1));Vg=round(Urbt(2));Wg=round(Urbt(3));

xn=zeros(size(conn));
yn=zeros(size(conn));
zn=zeros(size(conn));
im0e=zeros([length(hh),length(hh),length(hh),prod(Nelems)]);
im1e=zeros([length(hh),length(hh),length(hh),prod(Nelems)]);
for i1=1:prod(Nelems)
                   found=find(conn(i1,:)>0);
            inods=conn(i1,found);
            xn(i1,:)=xo(inods);yn(i1,:)=yo(inods);zn(i1,:)=zo(inods);
           xg=round(mean(xn(i1,:)));yg=round(mean(yn(i1,:)));zg=round(mean(zn(i1,:)));
            xp=dec(1)+xg+hh;
            found=find((xp>=1)&(xp<=sizeim(1))&(xp+Ug>=1)&(xp+Ug<=sizeim(1)));
            xp=xp(found);
            yp=dec(2)+yg+hh;
            found=find((yp>=1)&(yp<=sizeim(2))&(yp+Vg>=1)&(yp+Vg<=sizeim(2)));
            yp=yp(found);
            zp=dec(3)+zg+hh;
           found=find((zp>=1)&(zp<=sizeim(3))&(zp+Wg>=1)&(zp+Wg<=sizeim(3)));
            zp=zp(found);
            im0e(:,:,:,i1)=im0(xp,yp,zp);
            im1e(:,:,:,i1)=jm3(xp+Ug,yp+Vg,zp+Wg);

end



parfor i1=1:prod(Nelems)
             im00=im0e(:,:,:,i1);
            im11=im1e(:,:,:,i1);
            [Urbt]=rbt3(im00,im11);
            Ue(i1)=Ug+Urbt(1);
            Ve(i1)=Vg+Urbt(2);
            We(i1)=Wg+Urbt(3);

   
end




end