function U=MedianFilter(Ui,mesh_file,lc)
persistent cond
check=0;
typ=1;
U=Ui;
if lc>0
load(mesh_file,'-mat','xo','yo','zo','Nnodes','conn','elt');

% seg=[1:size(conn,2)-1;2:size(conn,2)];
% conn=conn(:,seg(:));
% conn=reshape(conn',2,size(seg,2)*size(conn,1))';
% conn=unique(conn,'rows','stable');

Ux=Ui((1:prod(Nnodes)));
dflag=size(Ui,1)==3*prod(Nnodes);
uflag=size(Ui,1)>=2*prod(Nnodes);
if uflag
    Uy=Ui(prod(Nnodes)+(1:prod(Nnodes)));
end
if dflag
    Uz=Ui(2*prod(Nnodes)+(1:prod(Nnodes)));
end
if isempty(cond)
    if 0
        one=sparse(ones(prod(Nnodes),1));
        cx=abs(sparse(xo)*one'-one*sparse(xo'))<=lc;
        cy=abs(sparse(yo)*one'-one*sparse(yo'))<=lc;
        cond=cx&cy;
    else
        if 1
            cond=sparse(prod(Nnodes),prod(Nnodes));
            indi=conn(:);
            indi(conn(:)==0)=[];
            for ii=1:size(conn,2)
                indj=repmat(conn(:,ii),size(conn,2),1);
                indj(conn(:)==0)=[];
                cond=cond+sparse(indi(indj>0),indj(indj>0),1,prod(Nnodes),prod(Nnodes));
            end
            %             for il=1:lc-1
            %                 condo=cond;
            %                 for in=1:prod(Nnodes);
            %                     lin=sum(condo(:,condo(:,in)>0),2);
            %                     cond(lin>0,in)=1;
            %                 end
            %             end
            condo=cond;
            for il=1:lc-1
                cond=cond*condo;
            end
            if typ==0
                
                for in=1:prod(Nnodes)
                    inods=cond(:,in)>0;
                    distn=abs(xo(inods)-xo(in)+1i*(yo(inods)-yo(in))).^2;
                    lcc=sqrt(max(distn)/2);
                    distn=exp(-distn/(2*(lcc/3)^2));
                    distn=distn/sum(distn);
                    %                    cond(inods,in)=exp(-abs(xo(inods)-xo(in)+1i*(yo(inods)-yo(in))).^2/(2*(lc/3)^2))/(sqrt(2*pi)*lc/3);
                    
                    cond(inods,in)=distn;
                    
                end
            end
            
            if check
                for in=1:prod(Nnodes)
                    figure
                    plot(xo,yo,'+')
                    hold on
                    plot(xo(cond(:,in)>0),yo(cond(:,in)>0),'ro')
                    plot(xo(in),yo(in),'ks')
                    pause
                    close all
                end
            end
            %            cond=cond-diag(diag(cond));
        else
            cond=sparse(prod(Nnodes),prod(Nnodes));
            
            %indi=zeros(prod(Nnodes)*10);
            %indj=zeros(prod(Nnodes)*10);
            %nel=0;
            for in=1:prod(Nnodes)
                new=1;
                inods=in;
                while new
                    ielts=GetEltsFromNodes(conn,elt,inods);
                    inew=GetNodesFromElts(conn(ielts,:),elt(ielts));
                    %found=abs(xo(inew)-xo(in)+1i*(yo(inew)-yo(in)))<=lc;
                    found=(abs(xo(inew)-xo(in))<=lc)&(abs(yo(inew)-yo(in))<=lc);
                    new=sum(found)>length(inods);
                    if new, inods=inew;end
                end
                %            if in/10==round(in/10)
                %                figure
                %  plot(xo,yo,'x')
                %  hold on
                %  plot(xo(inods),yo(inods),'ro')
                %  plot(xo(in),yo(in),'ks')
                %  pause
                %            end
                if typ
                    %             indi(nel+(1:length(inods)))=inods;
                    cond(inods,in)=1;
                else
                    %             indi(nel+(1:length(inods)))=exp(-abs(xo(inods)-xo(in)+1i*(yo(inods)-yo(in))).^2/(2*(lc/3)^2))/(sqrt(2*pi)*lc/3);
                    cond(inods,in)=exp(-abs(xo(inods)-xo(in)+1i*(yo(inods)-yo(in))).^2/(2*(lc/3)^2))/(sqrt(2*pi)*lc/3);
                end
                %             indj(nel+(1:length(inods)))=in;
                %             nel=nel+length(inods);
            end
            %     indi((nel+1):end)=[];
            %      indj((nel+1):end)=[];
            %    cond=sparse(indi,indj,1,prod(Nnodes),prod(Nnodes));
        end
    end
end
for in=1:prod(Nnodes)
    %    new=1;
    %    inods=in;
    %     while new
    %         ielts=GetEltsFromNodes(conn,elt,inods);
    %         inew=GetNodesFromElts(conn(ielts,:),elt(ielts));
    %             found=abs(xo(inew)-xo(in)+1i*(yo(inew)-yo(in)))<=lc;
    % new=sum(found)>length(inods);
    % if new, inods=inew;end
    %     end
    %    inods=(abs(xo-xo(in)+1i*(yo-yo(in)))<=lc);
    %    inods=(abs(xo-xo(in))<=lc)&(abs(yo-yo(in))<=lc);
    inods=cond(:,in)>0;
    %         sum(inods)
    %                         figure
    %                 plot3(xo,yo,zo,'+')
    %                 hold on
    %                 plot3(xo(cond(:,in)>0),yo(cond(:,in)>0),zo(cond(:,in)>0),'ro')
    %                 pause
    %                 close all
    
    %                                if in/10==round(in/10)
    %                            figure
    %              plot(xo,yo,'x')
    %              hold on
    %              plot(xo(inods),yo(inods),'ro')
    %              plot(xo(in),yo(in),'ks')
    %              pause
    %              close all
    %
    %                        end
    
    
    
    %    inods(in)=0;
    if typ
        
     
        U(in)=median(Ux(inods));
        
        
%        U(in)=0.5*median(Ux(inods))+0.5*mean(Ux(inods));
%         tmp=Ux(inods);
%         if U(in)>Ux(in)
%             
%         U(in)=median(tmp(tmp>U(in)));
%         else
%         U(in)=median(tmp(tmp<U(in)));
%         end
            
       if uflag
            U(prod(Nnodes)+in)=median(Uy(inods));
%            U(prod(Nnodes)+in)=0.5*median(Uy(inods))+0.5*mean(Uy(inods));
%          tmp=Uy(inods);
%          if U(prod(Nnodes)+in)>Uy(in)
%             
%         U(prod(Nnodes)+in)=median(tmp(tmp>U(prod(Nnodes)+in)));
%         else
%          U(prod(Nnodes)+in)=median(tmp(tmp<U(prod(Nnodes)+in)));
%         end
       end
        if dflag
            U(2*prod(Nnodes)+in)=median(Uz(inods));
        end
    else
        U(in)=(Ux)'*cond(:,in);
        if uflag
            U(prod(Nnodes)+in)=Uy'*cond(:,in);
        end
        if dflag
            U(2*prod(Nnodes)+in)=(Uz)'*cond(:,in);
        end
    end
end
end
end