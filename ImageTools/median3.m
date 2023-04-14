function U=median3(Ui)
U=0*Ui;
for ix=1:size(Ui,1)
    indx=max(1,ix-1):min(size(Ui,1),ix+1);
    for iy=1:size(Ui,2)
    indy=max(1,iy-1):min(size(Ui,2),iy+1);
        for iz=1:size(Ui,3)
    indz=max(1,iz-1):min(size(Ui,3),iz+1);
    tmp=Ui(indx,indy,indz);
    U(ix,iy,iz)=median(tmp(:));
        end
    end
end
end