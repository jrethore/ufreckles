function gp2cell=ElementAverage(elt)

ngq=4;

indi=repmat((1:length(elt))',1,ngq);
val=repmat((elt==4)/ngq,1,ngq)+[(elt==3),repmat(0,length(elt),ngq-1)];
indj=ones(size(indi'));
found=find(val'>0);
indj(found)=1:length(found);
indj=indj';
foundt3=find(elt==3);
foundq4=find(elt==4);
foundt=[foundt3;foundq4];
val=val(foundt,:);
indi=indi(foundt,:);
indi=indi(foundt,:);
gp2cell=sparse(indi,indj,val);

end