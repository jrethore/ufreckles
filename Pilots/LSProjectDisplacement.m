function [Xf]=LSProjectDisplacement(Xg,iscale,nmod)
load(fullfile('TMP','params'));
onflight=0;
if isfield(param,'onflight')
    onflight=param.onflight;
end
    persistent M phi
if ~onflight
    persistent Xi Yi wdetJ
else
    global phiy phix Xi Yi wdetJ inde on
end
load(fullfile('TMP',sprintf('%d_params',nmod)),'param');

if strcmp(param.basis,'nurbs')||strcmp(param.basis,'btri')
    load(fullfile('TMP',sprintf('%d_phio_%d',nmod,iscale-1)),'phio');
    phi=phio;
    Ng=size(Xg,1)/2;
    Ug=Xg((1:Ng),:);
    Vg=Xg(Ng+(1:Ng),:);
    M=phi'*phi;
    Fu=phi'*Ug;
    Uf=M\Fu;
    Fv=phi'*Vg;
    Vf=M\Fv;
    Xf=[Uf;Vf];
    % load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo');
    % figure
    % scatter(xo,yo,10+0*xo,Ug,'o')
    % figure
    % scatter(xo,yo,10+0*xo,phi*Uf,'o')
    % load(fullfile('TMP',sprintf('%d_phix_%d',nmod,iscale-1)),'phix','Xi','Yi');
    % figure
    % scatter(Xi,Yi,10+0*Xi,phix*[Uf;Vf],'o')
    
else
    load(fullfile('TMP',sprintf('%d_mesh_%d',nmod,iscale-1)),'xo','yo','Nnodes');
    Ng=length(Xg)/2;
    Ug=reshape(Xg((1:Ng)),Nnodes);
    Vg=reshape(Xg(Ng+(1:Ng)),Nnodes);
    % figure
    % scatter(xo,yo,'bo')
    % hold on
    
    xo=reshape(xo,Nnodes);
    yo=reshape(yo,Nnodes);
    xo=xo(1:Nnodes(1),1)';
    yo=yo(1,1:Nnodes(2));
    %
    %
    %
    % [Yi,Xi]=meshgrid(yo(1):yo(length(yo))-1,xo(1):xo(length(xo))-1);
    % Ug=interp2(yo,xo,Ug,Yi(:),Xi(:),'*linear');
    % Vg=interp2(yo,xo,Vg,Yi(:),Xi(:),'*linear');
    % sizeg=size(Xi);
    % Xg=Ug+i*Vg;
    % clear Xi Yi Ug Vg
    if isempty(M)
        if onflight
        load(fullfile('TMP',sprintf('%d_phix_%d',nmod,10*(iscale-1))),'Nddl_tot','sizeim');
        else
            load(fullfile('TMP',sprintf('%d_phix_%d',nmod,10*(iscale-1))),'wdetJ','phix','Nddl_tot','sizeim');
        load(fullfile('TMP',sprintf('%d_phiy_%d',nmod,10*(iscale-1))),'phiy');
        end
        %scatter(Xi,Yi,'rx')
        if length(wdetJ)==1
            xi=1:10:sizeim(1);
            yi=1:10:sizeim(2);
            [Yi,Xi]=meshgrid(yi,xi);
            [indj]=sub2ind(sizeim,Xi(:),Yi(:));
            P=sparse(1:length(indj),indj,1,numel(Xi),size(phix,1));
            wdetJ=1;
            phi=P*(phix+i*phiy);
            M=real(phi'*wdetJ*phi);
            wdetJ=1;
        else
            if ~onflight
                load(fullfile('TMP',sprintf('%d_phix_%d',nmod,10*(iscale-1))),'Xi','Yi');
            end
            phi=(phix+i*phiy);
            M=real(phi'*wdetJ*phi);
            
            
        end
    end
    Ug=interp2(yo,xo,Ug,min(max(yo),max(min(yo),Yi(:))),min(max(xo),max(min(xo),Xi(:))),'*linear');
    Vg=interp2(yo,xo,Vg,min(max(yo),max(min(yo),Yi(:))),min(max(xo),max(min(xo),Xi(:))),'*linear');
    Xg=Ug+i*Vg;
    
    % select=ones(sizeim);
    % select((sizeg(1)+1):sizeim(1),:)=0;
    % select(:,(sizeg(2)+1):sizeim(2))=0;
    % select=find(select);
    % select=sparse(1:prod(sizeg),select,1,prod(sizeg),prod(sizeim));
    % phi=select*phi;
    Xf=M\real(phi'*wdetJ*Xg);
    %
    % figure
    %  imagesc(reshape(Ug,sizeim))
    %  colorbar
    %  title('Ug')
    % figure
    %  imagesc(reshape(real(phix*Xf),sizeim))
    %  colorbar
    %  title('Uf')
    % figure
    % imagesc(reshape(real(phi*Xf)-Ug,sizeg))
    % colorbar
    % title('DIFF')
    %
    
end
end