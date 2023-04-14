function [sfils]=SortImageFiles(fils)


sfils=sort(fils);

nim=length(fils);

icam=zeros(nim,1);
nums=1:nim;

for iim=1:nim
   [pp,c,ext]=fileparts(fils{iim});
if strcmp(c(end-1),'_')
    icam(iim)=str2num(c(end));

else
    ids=textscan(c,'%*s %d %d','delimiter','-');
    if isempty(ids{1})|| isempty(ids{2})
    icam(iim)=1;
    else
    icam(iim)=ids{2};
    end 
end
end
 icam=icam-min(icam)+1;
 if ~(max(icam)==min(icam))
     ncam=max(icam)-min(icam)+1;
     sfils=reshape(sfils',ncam,nim/ncam);
 end
if nim==1&&ncam>2
    sfils=sfils';
end
% for iim=1:nim
%    [pp,c,ext]=fileparts(fils{iim});
% if strcmp(c(end-1),'_')
%     icam(iim)=str2num(c(end));
%     dec=2;
% elseif strcmp(c(1:3),'Cam')
%      icam(iim)=str2num(c(4));
%      dec=0;
% elseif strcmp(c(end-2),'_')&&strcmp(c(end-5),'_')
%     icam(iim)=1;
%     for ic=1:length(c)
%         if strcmp(c(ic),'_')
%             dec=length(c)-ic+1;
%             break
%         end
%     end
% else
%     icam(iim)=1;dec=0;
% end
% 
% ii=0;num=0;
% while ii<length(c)
%     cnum=str2num(c(end-ii-dec));
%     if ~isempty(cnum)
%         num=num+cnum*10^(ii);
%     else
%         break
%     end
%     ii=ii+1;
% end
% nums(iim)=num; 
%     
% end
% icam=icam-min(icam)+1;
% if max(icam)==min(icam)
%     [smun,ind]=sort(nums);
%     sfils=fils(ind);
% else
%     ncam=max(icam)-min(icam)+1;
%     sfils=cell(ncam,nim/ncam+1);
%     snums=unique(sort(nums));
%     for iim=1:nim
%         ijm=find(snums==nums(iim));
%         sfils(icam(iim),ijm+1)=fils(iim);
%         
%     end
%     sfils(:,1)=sfils(:,2);
% end


end