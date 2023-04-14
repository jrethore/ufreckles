function CreateArlequinMasks(iscale,nmod1,nmod2)

load(fullfile('TMP','params'),'param');
param0=param;
load(fullfile('TMP',sprintf('%d_params',nmod2)),'param');
check=1;
    %%
    if isfield(param,'crack_id')
        ic=param.crack_id;
    else
        ic=1;
    end
    filfis=fullfile('TMP',sprintf('%d_levelsets',ic));
    rmin=param0.mask_radius;
    mesh_file=fullfile('TMP',sprintf('%d_mesh_%d',nmod1,iscale-1));
    load(mesh_file,'elt','conn','xo','yo','Nnodes','Smesh','mesh_size');
    load(fullfile('TMP',sprintf('%d_mask_%d',nmod1,iscale-1)),'maskn','unmasked_nodes','cut_nodes');
    if isfield(param0,'mask_radius')
        rmin=param0.mask_radius;
    else
        rmin=mean(mesh_size);
    end
    if isfield(param0,'coupling_width')
        nelem=param0.coupling_width;
    else
        nelem=1;
    end
    if isfield(param0,'distance_type')
        if strcmp(param0.distance_type,'chemical')
            typ=1;
        else
            typ=0;
        end
    else
        typ=0;
    end

    if ~typ
        Yo=yo;Xo=xo;
        load(filfis,'front');
        sizeim=size(front);
        front=interp2(front,Yo,Xo,'*linear');
        load(filfis,'crack');
        crack=interp2(crack,Yo,Xo,'*linear');
        dist=abs(front+i*crack);

        outside=(dist<=rmin);

        if isfield(param0,'mask_width')
            dmin=param0.mask_width;
            outside=outside|((abs(crack)<=dmin)&(front<0));
        end

        if isfield(param,'mask_length')||isfield(param,'cz_length')
            lmin=0;
            if isfield(param,'mask_length')
                lmin=lmin+param.mask_length;
            end
            if isfield(param,'cz_length')
                lmin=lmin+param.cz_length;
            end
            outside=outside|((front<=lmin)&(front>0)&(abs(crack)<=rmin))|(abs(front-lmin+i*crack)<=rmin);
        end
    else
        Yo=sum(yo(conn),2)./elt;
        Xo=sum(xo(conn),2)./elt;

        load(filfis,'front');
        sizeim=size(front);
        front=interp2(front,Yo,Xo,'*linear');
        load(filfis,'crack');
        crack=interp2(crack,Yo,Xo,'*linear');
        dist=abs(front+i*crack);
        outside=zeros(Nnodes);
        tip_elem=find(dist==min(dist(:)),1,'first');
        outside_nodes=conn(tip_elem,:);

        for ilayer=1:rmin-1
            outside_nodes=AddOneNodeLayer(conn,outside_nodes);
        end
        outside(outside_nodes)=1;

        if isfield(param0,'mask_width')
            dmin=param0.mask_width;
            crack_elems=find(((abs(crack)<=mean(mesh_size))&(front<0)));
            is_crack_nodes=zeros(prod(Nnodes),1);
            crack_nodes=conn(crack_elems,:);
            keep=find(crack_nodes>0);
            is_crack_nodes(crack_nodes(keep))=1;
            crack_nodes=find(is_crack_nodes);
            for ilayer=1:dmin-1
                crack_nodes=AddOneNodeLayer(conn,crack_nodes);
            end
            outside(crack_nodes)=1;

        end

        if isfield(param,'mask_length')||isfield(param,'cz_length')
            lmin=0;
            if isfield(param,'mask_length')
                lmin=lmin+param.mask_length;
            end
            if isfield(param,'cz_length')
                lmin=lmin+param.cz_length;
            end
            cz_elems=[find((front<=lmin)&(front>0)&(abs(crack)<=mean(mesh_size)));...
                find(abs(front-lmin+i*crack)==min(abs(front(:)-lmin+i*crack(:))),1,'first')];
            is_cz_nodes=zeros(prod(Nnodes),1);
            cz_nodes=conn(cz_elems,:);
            keep=find(cz_nodes>0);
            is_cz_nodes(cz_nodes(keep))=1;
            cz_nodes=find(is_cz_nodes);
            for ilayer=1:rmin-1
                cz_nodes=AddOneNodeLayer(conn,cz_nodes);
            end
            outside(cz_nodes)=1;
        end



    end
    ln=double(~outside);

    ln0=ln;


    
    outside_nodes=find(outside(:));
    for ilayer=1:nelem
        outside_nodes=AddOneNodeLayer(conn,outside_nodes);
    end
    outside(outside_nodes)=1;


    inside=ln&outside;



    inside_nodes=find(inside(:));
    selected_nodes=AddOneNodeLayer(conn,inside_nodes);

    isselected=zeros(prod(Nnodes),1);
    isselected(selected_nodes)=1;

    selected_nodes=find(isselected(:));
    boundary_nodes=find(isselected(:)&(~inside(:)));

    if check
        figure
        hold on;
        plot(xo(inside_nodes),yo(inside_nodes),'bs','MarkerSize',7);
        plot(xo(find(isselected&(~inside))),yo(find(isselected&(~inside))),'ro','MarkerSize',7);
        plot(xo(find(isselected)),yo(find(isselected)),'kx');
        legend({'Inside nodes','Selected but not inside','Selected'})
plot(xo,yo,'k+')
        axis equal
        axis off;
        title('Coupling','FontSize',24);
    end



    %%
    Kp=AssemblePenaltyMatrix(iscale,nmod1 ,selected_nodes,prod(Nnodes));
    C=sparse(1:length(boundary_nodes),boundary_nodes,1,length(boundary_nodes),prod(Nnodes));
    Kp=Kp(selected_nodes,selected_nodes);
    C=C(:,selected_nodes);
    K=[Kp,C';C,sparse(length(boundary_nodes),length(boundary_nodes))];
    Tp=ln(boundary_nodes);
    Q=[zeros(length(selected_nodes),1);Tp];
    T=K\Q;
    T=T(1:length(selected_nodes));

    ln(selected_nodes)=T;
    ln(boundary_nodes)=Tp;
    ln=ln.*double(ln0);
   
    clear T Q K Kp
    maskn=(maskn).*(ln);
    
    
    
    
    
    phi=CreateFiniteElementBasis(mesh_file,sizeim);
    if check
        PlotMap(reshape(phi*maskn,sizeim),'Model 1: modified mask',[],[param0.result_file,'-mask1'],0);
    end

    phis=CreateFiniteElementBasis(mesh_file,sizeim,1,[],'sub_cells');
    mask1s=phis*ln(:);
    lps=diag(sparse(double(mask1s)));
    phig=CreateFiniteElementBasis(mesh_file,sizeim,1,[],'Gauss_points');
    mask1=phi*ln(:);
    mask1g=phig*ln(:);
    lp=diag(sparse(double(mask1)));
    lpg=diag(sparse(double(mask1g)));
    load(fullfile('TMP',sprintf('%d_mask_%d',nmod1,iscale-1)),'mask','masks','maskg');
    mask=(mask*lp);
    maskg=(maskg*lpg);
    masks=(masks*lps);
     tmp=maskn>0;
    unmasked_nodes=find(tmp(:));
        unmasked_nodes=AddOneNodeLayer(conn,unmasked_nodes);
%     figure
%     scatter(phis*xo,phis*yo,10+0*phis*yo,diag(masks));
    

    save(fullfile('TMP',sprintf('%d_mask_%d',nmod1,iscale-1)),'mask','maskg','masks','unmasked_nodes','cut_nodes','maskn');




    load(fullfile('TMP',sprintf('%d_mask_%d',nmod2,iscale-1)),'mask');
    mask2=phi*(1-ln(:));
    lp=diag(sparse(double(mask2)));
    mask=mask*lp;
    maskg=mask;
    
    save(fullfile('TMP',sprintf('%d_mask_%d',nmod2,iscale-1)),'mask','maskg');



    coupling_nodes=selected_nodes;
    save(fullfile('TMP',sprintf('alrequin_%d',nmod2)),'coupling_nodes','inside_nodes');


end