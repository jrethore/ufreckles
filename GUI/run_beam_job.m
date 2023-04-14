function run_beam_job(paramo,modelo)
global phiy phix Xi Yi phidf on wdetJ
nmod=0;
roio=paramo.roi;
param=paramo;
filres=paramo.result_file;
filreso=strrep(filres,'.res','');
gage=modelo.zone{2,1};
param.onflight=1;
model=modelo;
if ~isfield(model,'beam_type');
    model.beam_type='euler';
end
if isfield(model,'beam_type')
    beamtype=strcmp(model.beam_type,'timoshenko');
else
    beamtype=0;
end

l1=abs(diff(gage(1:2,:)*[1;1i]));
l2=abs(diff(gage(2:3,:)*[1;1i]));
if l2>l1
    x1=gage(2:3,1);
    y1=gage(2:3,2);
    x2=gage([5,4],1);
    y2=gage([5,4],2);
else
    x1=gage(1:2,1);
    y1=gage(1:2,2);
    x2=gage([3,4],1);
    y2=gage([3,4],2);
end
x=[x1;x2];y=[y1;y2];
x=sort(round(x));y=sort(round(y));
roi(1) = min(x);
roi(2) = max(x);
roi(3) = min(y);
roi(4) = max(y);
param.roi=roi;
sizeim=[roi(2)-roi(1),roi(4)-roi(3)];

nmod=0;
LoadParameters(param);
LoadParameters(model,nmod);
ReferenceImage(nmod);
LoadMask(nmod);
nscale=1;

filres=param.result_file;
t1=[diff(x1);diff(y1)];
t1=sign(t1(abs(t1)==max(abs(t1))))*t1/norm(t1);
t2=[diff(x2);diff(y2)];
t2=t2/norm(t2);
t2=t2*sign(t1'*t2);
t=0.5*(t1+t2);
n=[-t(2);t(1)];
xc1=mean(x1);
xc2=mean(x2);
xc=0.5*(xc1+xc2);
yc1=mean(y1);
yc2=mean(y2);
yc=0.5*(yc1+yc2);

x=[x1;x2];y=[y1;y2];
x=sort(round(x));y=sort(round(y));
U=[];
for iscale=nscale:-1:1
    
    disp(sprintf('Pre-processing scale %d...',iscale));
    load(fullfile('TMP',sprintf('sample%d_%d',0,iscale-1)),'sizeim');
    pscale=2^(iscale-1);
    [Yi,Xi]=meshgrid(1:sizeim(2),1:sizeim(1));
    X=(Xi*pscale-1)+1+roi(1)-1;
    Y=(Yi*pscale-1)+1+roi(3)-1;
    S=X*t(1)+Y*t(2);
    s=sort(x*t(1)+y*t(2));
    T=X*n(1)+Y*n(2);
    sc1=xc1*t(1)+yc1*t(2);
    sc2=xc2*t(1)+yc2*t(2);
    tc1=xc1*n(1)+yc1*n(2);
    tc2=xc2*n(1)+yc2*n(2);
    sc=0.5*(sc1+sc2);
    tc=0.5*(tc1+tc2);
    mask=(T<max(tc1,tc2))&(T>min(tc1,tc2))&(S>s(2))&(S<s(3));
    mask=diag(sparse(double(mask(:))));
    T=2*(T(:)-min(tc1,tc2))/abs(tc1-tc2)-1;
    S=(S(:)-s(2))/(s(3)-s(2))-0.5;
    dsdx=1/(s(3)-s(2));
    dtdx=2/abs(tc1-tc2);
    ST=dsdx/dtdx;
    T=T*ST;
    switch model.basis
        case 'beam'
            phis=zeros(prod(sizeim), 4+double(model.exx));
            phit=zeros(prod(sizeim), 4+double(model.exx));
            phis(:,1)=1;
            
            phit(:,2)=1;
            phit(:,3)=S;
            phis(:,3)=-T;
            phit(:,4)=0.5*S.^2;
            phis(:,4)=-(T.*S);
            if model.exx
                phis(:,5)=S;
            end
            
            
            %                         phit(:,6)=(1/6)*S.^3;
            %                         phis(:,6)=-0.5*(T.*(S.^2));
            phis=mask*phis;
            phit=mask*phit;
            phix=t(1)*phis+n(1)*phit;
            phiy=t(2)*phis+n(2)*phit;
            
            phias=sparse(prod(sizeim), 4+double(model.exx));
            phits=sparse(prod(sizeim), 4+double(model.exx));
            phiss=sparse(prod(sizeim), 4+double(model.exx));
            
            
            
            phias(:,4)=-T*dsdx;
            if model.exx, phias(:,5)=dsdx;end
            phias=mask*phias;
            phits=mask*phits;
            phiss=mask*phiss;
            
            
            
        case 'beam-nurbs'
            si=linspace(0,1,model.nb_element+1);
            u=si;
            u1=u;
            for id=1:model.degree
                u=[0,u,1];
            end
            for id=1:model.degree-1
                u1=[0,u1,1];
            end
            u=u-0.5;
            u1=u1-0.5;
            si=si-0.5;
            p=model.degree;
            Smesh=numel(S);
            indp=zeros((p+1)*Smesh,1);
            indn=zeros((p+1)*Smesh,1);
            val=zeros((p+1)*Smesh,1);
            dval=zeros((p+1)*Smesh,1);
            ddval=zeros((p+1)*Smesh,1);
            Nny=length(si)+p-1;
            
            nel=0;
            for iy=1:model.nb_element
                
                if iy==model.nb_element
                    found=find((S(:)>=si(iy))&(S(:)<=si(iy+1)));
                else
                    found=find((S(:)>=si(iy))&(S(:)<si(iy+1)));
                end
                sp=S(found);
                tp=T(found);
                Sel=length(sp);
                [N]=NURBSBasisFunc(iy+p,p,sp,u,3);
                for ip=1:(p(1)+1)
                    indn(nel+(1:Sel))=iy+ip-1;
                    indp(nel+(1:Sel))=found;
                    val(nel+(1:Sel))=N(:,ip,1);
                    dval(nel+(1:Sel))=-N(:,ip,2).*tp;
                    ddval(nel+(1:Sel))=-N(:,ip,3).*tp;
                    nel=nel+Sel;
                end
                
            end
            indp((nel+1):end)=[];
            indn((nel+1):end)=[];
            val((nel+1):end)=[];
            dval((nel+1):end)=[];
            ddval((nel+1):end)=[];
            phit=sparse(indp,indn,val,prod(sizeim),Nny);
            phis=sparse(indp,indn,dval,prod(sizeim),Nny);
            
            phias=sparse(indp,indn,dsdx*ddval,prod(sizeim),Nny);
            phits=sparse(prod(sizeim),Nny);
            phiss=sparse(prod(sizeim),Nny);
            if beamtype
                phit=[phit,sparse(prod(sizeim),Nny)];
                phis=[sparse(prod(sizeim),Nny),phis];
                phias=[sparse(prod(sizeim),Nny),phias];
                phits=sparse(prod(sizeim),2*Nny);
                phiss=sparse(prod(sizeim),2*Nny);
            end
            %                         Nelt=1;
            %                         if ~(Nelt==model.nb_element)
            %
            %                          si=linspace(0,1,1+1);
            %                         u=si;
            %                         u1=[0,u,1];
            %                         si=si-0.5;
            %
            %                         end
            if model.exx==1
                p1=p-1;
                indp=zeros((p1+1)*Smesh,1);
                indn=zeros((p1+1)*Smesh,1);
                val=zeros((p1+1)*Smesh,1);
                dval=zeros((p1+1)*Smesh,1);
                Nny=length(si)+p1-1;
                nel=0;
                for iy=1:model.nb_element
                    
                    if iy==model.nb_element
                        found=find((S(:)>=si(iy))&(S(:)<=si(iy+1)));
                    else
                        found=find((S(:)>=si(iy))&(S(:)<si(iy+1)));
                    end
                    sp=S(found);
                    Sel=length(sp);
                    [N]=NURBSBasisFunc(iy+p1,p1,sp,u1,2);
                    for ip=1:(p1+1)
                        indn(nel+(1:Sel))=iy+ip-1;
                        indp(nel+(1:Sel))=found;
                        val(nel+(1:Sel))=N(:,ip,1);
                        dval(nel+(1:Sel))=N(:,ip,2);
                        nel=nel+Sel;
                    end
                    
                end
                indp((nel+1):end)=[];
                indn((nel+1):end)=[];
                val((nel+1):end)=[];
                dval((nel+1):end)=[];
                
                phis=[phis,sparse(indp,indn,val,prod(sizeim),Nny)];
                phit=[phit,sparse(prod(sizeim),Nny)];
                
                phias=[phias,sparse(indp,indn,dsdx*dval,prod(sizeim),Nny)];
                phits=[phits,sparse(prod(sizeim),Nny)];
                phiss=[phiss,sparse(prod(sizeim),Nny)];
            else
                
                phis=[phis,sparse(1:prod(sizeim),1,1,prod(sizeim),1)];
                phit=[phit,sparse(prod(sizeim),1)];
                
                phias=[phias,sparse(prod(sizeim),1)];
                phits=[phits,sparse(prod(sizeim),1)];
                phiss=[phiss,sparse(prod(sizeim),1)];
                
            end
            
            phis=mask*phis;
            phit=mask*phit;
            phix=t(1)*phis+n(1)*phit;
            phiy=t(2)*phis+n(2)*phit;
            
            phias=mask*phias;
            phits=mask*phits;
            phiss=mask*phiss;
            
    end
    Xi=Xi(:);
    Yi=Yi(:);
    Nddl_tot=size(phix,2);
    wdetJ=1;
    on=sum(abs(phix+1i*phiy),2)>0;
    %                save(fullfile('TMP',sprintf('%d_phiy_%d',nmod,iscale-1)),'phiy','sizeim','Nddl_tot','-v7.3');
    %                save(fullfile('TMP',sprintf('%d_phix_%d',nmod,iscale-1)),'phix','Xi','Yi','wdetJ','sizeim','Nddl_tot','on','-v7.3');
    mphi=double(abs(mean(phix,1))>abs(mean(phiy,1)))';
    save(fullfile('TMP',sprintf('%d_mphi_%d',nmod,iscale-1)),'mphi');
    ComputeGradFPhi(iscale,nmod);
    AssembleCorrelationOperator(iscale,nmod);
    Uini=InitializeSolution(U,nscale,nmod);
    if iscale==nscale
        switch model.basis
            case 'beam'
                Uy=Uini(1,:);
                Ux=Uini(3,:);
                Us=Ux*t(1)+Uy*t(2);
                Ut=Ux*n(1)+Uy*n(2);
                Uini=0*Uini;
                Uini(1,:)=Us;
                Uini(2,:)=Ut;
                Uini=Uini(1:(4+double(model.exx)),:);
        end
    end
    [U]=Solve(Uini,iscale,nmod);
end
copyfile(fullfile('TMP','0_error_0.mat'), [filreso,'-error.res']);

delete([fullfile('VTK',[filreso,'-error']),'-0*.vtk']);
images=1:size(U,2);
for iim=1:size(U,2)
    movefile(fullfile('VTK',sprintf('camr-1-camd-1-scale-1-%d-error.vtk',iim)),fullfile('VTK',sprintf('%s-error-%04d.vtk',filreso,images(iim))));
end
VTKExportScalarMap([filreso,sprintf('-%d',0)],'error',roi(1)-1+(1:sizeim(1)),roi(3)-1+(1:sizeim(2)),zeros(sizeim),1);
movefile(fullfile('VTK',[filreso,sprintf('-%d',0),'-error.vtk']),fullfile('VTK',sprintf('%s-error-%04d.vtk',filreso,0)));
delete(fullfile('VTK',['camr*','-error.vtk']));


%%
TT=-1:0.01:1;
T=1;
%S=(-0.5:0.01:0.5)';
switch model.basis
    case 'beam'
        S=linspace(-0.5,0.5,100);
    case 'beam-nurbs'
        S=linspace(-0.5,0.5,max(100,model.degree*model.nb_element));
end
ss=sort(x*t(1)+y*t(2));
s=ss(2):(ss(3)-ss(2))/(length(S)-1):ss(3);
TT=TT*ST;
sizeim=size(S);
switch model.basis
    case 'beam'
        phis=zeros(prod(sizeim), 4+double(model.exx));
        phit=zeros(prod(sizeim), 4+double(model.exx));
        phis(:,1)=1;
        
        phit(:,2)=1;
        phit(:,3)=S;
        phis(:,3)=-T;
        phit(:,4)=0.5*S.^2;
        phis(:,4)=-(T.*S);
        
        fleche=phit*U;
        curv=U(4,:)*dsdx*dsdx;
        strain=-TT'*U(4,:)*dsdx;
        rot=[];
        
        naxis=zeros(1,size(U,2));
        if model.exx
            strain=strain+(1+0*TT')*U(5,:)*dsdx;
            naxis=naxis+0.5*(U(5,:)./U(4,:))/ST;
            %%-0.5*strain/dsdx/(curv/(dsdx*dsdx)))*dtdx/dsdx
            %%-0.5*strain/curv*dtdx
        end
        
        fid=fopen([filreso,'-deflection.csv'],'w');
        fprintf(fid,'Result file;%s\n',param.result_file);
        fprintf(fid,'Reference image;%s\n',param.reference_image);
        fprintf(fid,'Beam;Constant curvature\n');
        fprintf(fid,'Position s [pixel];Deflection [pixel]\n');
        fprintf(fid,'Image file');
        if size(U,2)>1
            for iim=1:size(U,2)
                fprintf(fid,';"%s"',param.deformed_image{iim});
            end
        else
            fprintf(fid,';"%s"\n',param.deformed_image);
        end
        fprintf(fid,'\n');
        fprintf(fid,'Image number');
        for iim=1:size(U,2)
            fprintf(fid,';%d',iim);
        end
        fprintf(fid,'\n');
        for ix=1:length(s)
            fprintf(fid,'%f',s(ix));
            for iim=1:size(U,2)
                fprintf(fid,';%.3e',fleche(ix,iim));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        
        
        
        
        
        fid=fopen([filreso,'-curvature.csv'],'w');
        fprintf(fid,'Result file;%s\n',param.result_file);
        fprintf(fid,'Reference image;%s\n',param.reference_image);
        fprintf(fid,'Beam;Constant curvature\n');
        fprintf(fid,'File;Image;Curvature [1/pixel]\n');
        if size(U,2)>1
            for iim=1:size(U,2)
                fprintf(fid,'"%s";%d;%.3e\n',param.deformed_image{iim},iim,curv(iim));
            end
        else
            fprintf(fid,'"%s";%d;%.3e\n',param.deformed_image,1,curv);
        end
        fclose(fid);
        
        
        
        
        
        
        fid=fopen([filreso,'-strain.csv'],'w');
        fprintf(fid,'Result file;%s\n',param.result_file);
        fprintf(fid,'Reference image;%s\n',param.reference_image);
        fprintf(fid,'Beam;Constant curvature\n');
        fprintf(fid,'Position t/h [];Axial strain []\n');
        fprintf(fid,'Image file');
        if size(U,2)>1
            for iim=1:size(U,2)
                fprintf(fid,';"%s"',param.deformed_image{iim});
            end
        else
            fprintf(fid,';"%s"\n',param.deformed_image);
        end
        
        fprintf(fid,'\n');
        fprintf(fid,'Image number');
        for iim=1:size(U,2)
            fprintf(fid,';%d',iim);
        end
        fprintf(fid,'\n');
        for ix=1:length(TT)
            fprintf(fid,'%f',0.5*TT(ix)/ST);
            for iim=1:size(U,2)
                fprintf(fid,';%.3e',strain(ix,iim));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        fid=fopen([filreso,'-neutral-axis.csv'],'w');
        fprintf(fid,'Result file;%s\n',param.result_file);
        fprintf(fid,'Reference image;%s\n',param.reference_image);
        fprintf(fid,'Beam;Constant curvature\n');
        fprintf(fid,'File;Image;Normalized neutral axis position[]\n');
        if size(U,2)>1
            for iim=1:size(U,2)
                fprintf(fid,'"%s";%d;%.3e\n',param.deformed_image{iim},iim,naxis(iim));
            end
        else
            fprintf(fid,'"%s";%d;%.3e\n',param.deformed_image,1,naxis);
        end
        fclose(fid);
        
    case 'beam-nurbs'
        si=linspace(0,1,model.nb_element+1);
        u=si;
        u1=u;
        for id=1:model.degree
            u=[0,u,1];
        end
        for id=1:model.degree-1
            u1=[0,u1,1];
        end
        u=u-0.5;
        u1=u1-0.5;
        si=si-0.5;
        p=model.degree;
        Smesh=numel(S);
        indp=zeros((p+1)*Smesh,1);
        indn=zeros((p+1)*Smesh,1);
        val=zeros((p+1)*Smesh,1);
        dval=zeros((p+1)*Smesh,1);
        ddval=zeros((p+1)*Smesh,1);
        Nny=length(si)+p-1;
        
        nel=0;
        for iy=1:model.nb_element
            
            if iy==model.nb_element
                found=find((S(:)>=si(iy))&(S(:)<=si(iy+1)));
            else
                found=find((S(:)>=si(iy))&(S(:)<si(iy+1)));
            end
            sp=S(found);
            Sel=length(sp);
            [N]=NURBSBasisFunc(iy+p,p,sp,u,2);
            for ip=1:(p(1)+1)
                indn(nel+(1:Sel))=iy+ip-1;
                indp(nel+(1:Sel))=found;
                val(nel+(1:Sel))=N(:,ip,1);
                dval(nel+(1:Sel))=N(:,ip,2);
                ddval(nel+(1:Sel))=N(:,ip,3);
                nel=nel+Sel;
            end
            
        end
        indp((nel+1):end)=[];
        indn((nel+1):end)=[];
        val((nel+1):end)=[];
        dval((nel+1):end)=[];
        ddval((nel+1):end)=[];
        
        phit=sparse(indp,indn,val,prod(sizeim),Nny);
        dphit=sparse(indp,indn,dval,prod(sizeim),Nny);
        ddphit=sparse(indp,indn,ddval,prod(sizeim),Nny);
        
        
        if beamtype
            phit=[phit,sparse(prod(sizeim),Nny),sparse(prod(sizeim),size(U,1)-2*Nny)];
            dphit=[sparse(prod(sizeim),Nny),dphit,sparse(prod(sizeim),size(U,1)-2*Nny)];
            ddphit=[sparse(prod(sizeim),Nny),ddphit,sparse(prod(sizeim),size(U,1)-2*Nny)];
            phixx=sparse(prod(sizeim),size(U,1));
        else
            phit=[phit,sparse(prod(sizeim),size(U,1)-Nny)];
            dphit=[dphit,sparse(prod(sizeim),size(U,1)-Nny)];
            ddphit=[ddphit,sparse(prod(sizeim),size(U,1)-Nny)];
            phixx=sparse(prod(sizeim),size(U,1));
        end
        fleche=phit*U;
        rot=dsdx*dphit*U;
        curv=dsdx*dsdx*ddphit*U;
        
        if model.exx
            p1=p-1;
            indp=zeros((p1+1)*Smesh,1);
            indn=zeros((p1+1)*Smesh,1);
            dval=zeros((p1+1)*Smesh,1);
            Nny=length(si)+p1-1;
            nel=0;
            for iy=1:model.nb_element
                
                if iy==model.nb_element
                    found=find((S(:)>=si(iy))&(S(:)<=si(iy+1)));
                else
                    found=find((S(:)>=si(iy))&(S(:)<si(iy+1)));
                end
                sp=S(found);
                Sel=length(sp);
                [N]=NURBSBasisFunc(iy+p1,p1,sp,u1,2);
                for ip=1:(p1+1)
                    indn(nel+(1:Sel))=iy+ip-1;
                    indp(nel+(1:Sel))=found;
                    dval(nel+(1:Sel))=N(:,ip,2);
                    nel=nel+Sel;
                end
                
            end
            indp((nel+1):end)=[];
            indn((nel+1):end)=[];
            dval((nel+1):end)=[];
            phixx(:,(end-(Nny-1)):end)=sparse(indp,indn,dsdx*dval,prod(sizeim),Nny);
        end
        naxis=0.5*(phixx*U)./curv*dtdx;
        
        
        
        fid=fopen([filreso,'-deflection.csv'],'w');
        fprintf(fid,'Result file;%s\n',param.result_file);
        fprintf(fid,'Reference image;%s\n',param.reference_image);
        fprintf(fid,'Beam;NURBS degree;%d;NURBS elements;%d\n',model.degree,model.nb_element);
        fprintf(fid,'Beam type;%s\n',model.beam_type);
        fprintf(fid,'Position s [pixel];Deflection [pixel]\n');
        fprintf(fid,'Image file');
        if size(U,2)>1
            for iim=1:size(U,2)
                fprintf(fid,';"%s"',param.deformed_image{iim});
            end
        else
            fprintf(fid,';"%s"\n',param.deformed_image);
        end
        fprintf(fid,'\n');
        fprintf(fid,'Image number');
        for iim=1:size(U,2)
            fprintf(fid,';%d',iim);
        end
        fprintf(fid,'\n');
        for ix=1:length(s)
            fprintf(fid,'%f',s(ix));
            for iim=1:size(U,2)
                fprintf(fid,';%.3e',fleche(ix,iim));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        fid=fopen([filreso,'-rotation.csv'],'w');
        fprintf(fid,'Result file;%s\n',param.result_file);
        fprintf(fid,'Reference image;%s\n',param.reference_image);
        fprintf(fid,'Beam;NURBS degree;%d;NURBS elements;%d\n',model.degree,model.nb_element);
        fprintf(fid,'Beam type;%s\n',model.beam_type);
        fprintf(fid,'Position s [pixel];Rotation []\n');
        fprintf(fid,'Image file');
        if size(U,2)>1
            for iim=1:size(U,2)
                fprintf(fid,';"%s"',param.deformed_image{iim});
            end
        else
            fprintf(fid,';"%s"\n',param.deformed_image);
        end
        fprintf(fid,'\n');
        fprintf(fid,'Image number');
        for iim=1:size(U,2)
            fprintf(fid,';%d',iim);
        end
        fprintf(fid,'\n');
        for ix=1:length(s)
            fprintf(fid,'%f',s(ix));
            for iim=1:size(U,2)
                fprintf(fid,';%.3e',rot(ix,iim));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        fid=fopen([filreso,'-curvature.csv'],'w');
        fprintf(fid,'Result file;%s\n',param.result_file);
        fprintf(fid,'Reference image;%s\n',param.reference_image);
        fprintf(fid,'Beam;NURBS degree;%d;NURBS elements;%d\n',model.degree,model.nb_element);
        fprintf(fid,'Beam type;%s\n',model.beam_type);
        fprintf(fid,'Position s [pixel];Curvature [1/pixel]\n');
        fprintf(fid,'Image file');
        if size(U,2)>1
            for iim=1:size(U,2)
                fprintf(fid,';"%s"',param.deformed_image{iim});
            end
        else
            fprintf(fid,';"%s"\n',param.deformed_image);
        end
        fprintf(fid,'\n');
        fprintf(fid,'Image number');
        for iim=1:size(U,2)
            fprintf(fid,';%d',iim);
        end
        fprintf(fid,'\n');
        for ix=1:length(s)
            fprintf(fid,'%f',s(ix));
            for iim=1:size(U,2)
                fprintf(fid,';%.3e',curv(ix,iim));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        fid=fopen([filreso,'-neutral-axis.csv'],'w');
        fprintf(fid,'Result file;%s\n',param.result_file);
        fprintf(fid,'Reference image;%s\n',param.reference_image);
        fprintf(fid,'Beam;NURBS degree;%d;NURBS elements;%d\n',model.degree,model.nb_element);
        fprintf(fid,'Beam type;%s\n',model.beam_type);
        fprintf(fid,'Position s [pixel];Normalized neutral axis position []\n');
        fprintf(fid,'Image file');
        if size(U,2)>1
            for iim=1:size(U,2)
                fprintf(fid,';"%s"',param.deformed_image{iim});
            end
        else
            fprintf(fid,';"%s"\n',param.deformed_image);
        end
        fprintf(fid,'\n');
        fprintf(fid,'Image number');
        for iim=1:size(U,2)
            fprintf(fid,';%d',iim);
        end
        fprintf(fid,'\n');
        for ix=1:length(s)
            fprintf(fid,'%f',s(ix));
            for iim=1:size(U,2)
                fprintf(fid,';%.3e',naxis(ix,iim));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        
end
switch model.basis
    case 'beam'
        so=linspace(-0.5,0.5,10);
    case 'beam-nurbs'
        so=linspace(-0.5,0.5,model.degree*model.nb_element);
end
to=[-1:0.25:1];
Nnodes=[length(so),length(to),1];
Nelems=max(Nnodes-1,1);
[to,so]=meshgrid(to,so);
so=so(:);to=to(:);
soi=(so+0.5)*(ss(3)-ss(2))+ss(2);
toi=0.5*(to+1)*abs(tc1-tc2)+min(tc1,tc2);
xo=soi*t(1)+toi*n(1);
yo=soi*t(2)+toi*n(2);
zo=1;
incn=[0,1,Nnodes(1)+1,Nnodes(1)];
nroot=reshape(1:prod(Nnodes),Nnodes);
nroot=nroot(1:max(1,Nnodes(1)-1),:);
nroot=nroot(:,1:max(1,Nnodes(2)-1));
conn=repmat(nroot(:),1,4)+repmat(incn,prod(Nelems),1);
elt=repmat(4,prod(Nelems),1);

to=to*ST;

sizeim=size(xo);
switch model.basis
    case 'beam'
        phis=zeros(prod(sizeim), 4+double(model.exx));
        phit=zeros(prod(sizeim), 4+double(model.exx));
        phis(:,1)=1;
        
        phit(:,2)=1;
        phit(:,3)=so;
        phis(:,3)=-to;
        phit(:,4)=0.5*so.^2;
        phis(:,4)=-(to.*so);
        phixs=t(1)*phis+n(1)*phit;
        phiys=t(2)*phis+n(2)*phit;
        
        phias=sparse(prod(sizeim), 4+double(model.exx));
        phits=sparse(prod(sizeim), 4+double(model.exx));
        phiss=sparse(prod(sizeim), 4+double(model.exx));
        
        
        
        phias(:,4)=-to*dsdx;
        if model.exx, phias(:,5)=dsdx;end
        
    case 'beam-nurbs'
        S=so;T=to;
        sizeim=size(S);
        
        si=linspace(0,1,model.nb_element+1);
        u=si;
        u1=u;
        for id=1:model.degree
            u=[0,u,1];
        end
        for id=1:model.degree-1
            u1=[0,u1,1];
        end
        u=u-0.5;
        u1=u1-0.5;
        si=si-0.5;
        p=model.degree;
        Smesh=numel(S);
        indp=zeros((p+1)*Smesh,1);
        indn=zeros((p+1)*Smesh,1);
        val=zeros((p+1)*Smesh,1);
        dvalo=zeros((p+1)*Smesh,1);
        dval=zeros((p+1)*Smesh,1);
        ddval=zeros((p+1)*Smesh,1);
        Nny=length(si)+p-1;
        
        nel=0;
        for iy=1:model.nb_element
            
            if iy==model.nb_element
                found=find((S(:)>=si(iy))&(S(:)<=si(iy+1)));
            else
                found=find((S(:)>=si(iy))&(S(:)<si(iy+1)));
            end
            sp=S(found);
            tp=T(found);
            Sel=length(sp);
            [N]=NURBSBasisFunc(iy+p,p,sp,u,3);
            for ip=1:(p(1)+1)
                indn(nel+(1:Sel))=iy+ip-1;
                indp(nel+(1:Sel))=found;
                val(nel+(1:Sel))=N(:,ip,1);
                dval(nel+(1:Sel))=-N(:,ip,2).*tp;
                dvalo(nel+(1:Sel))=N(:,ip,2);
                ddval(nel+(1:Sel))=-N(:,ip,3).*tp;
                nel=nel+Sel;
            end
            
        end
        indp((nel+1):end)=[];
        indn((nel+1):end)=[];
        val((nel+1):end)=[];
        dval((nel+1):end)=[];
        ddval((nel+1):end)=[];
        
        phit=sparse(indp,indn,val,prod(sizeim),Nny);
        phis=sparse(indp,indn,dval,prod(sizeim),Nny);
        
        phias=sparse(indp,indn,(dsdx)*ddval,prod(sizeim),Nny);
        phits=sparse(prod(sizeim),Nny);
        phiss=sparse(prod(sizeim),Nny);
        
        if beamtype
            phit=[phit,sparse(prod(sizeim),Nny)];
            phis=[sparse(prod(sizeim),Nny),phis];
            phias=[sparse(prod(sizeim),Nny),phias];
            phits=sparse(prod(sizeim),2*Nny);
            phiss=sparse(indp,indn,dvalo*dsdx*0.5,prod(sizeim),Nny);
            phiss=[phiss,-phiss];
        end
        
        
        
        if model.exx
            p1=p-1;
            indp=zeros((p1+1)*Smesh,1);
            indn=zeros((p1+1)*Smesh,1);
            val=zeros((p1+1)*Smesh,1);
            dval=zeros((p1+1)*Smesh,1);
            Nny=length(si)+p1-1;
            nel=0;
            for iy=1:model.nb_element
                
                if iy==model.nb_element
                    found=find((S(:)>=si(iy))&(S(:)<=si(iy+1)));
                else
                    found=find((S(:)>=si(iy))&(S(:)<si(iy+1)));
                end
                sp=S(found);
                Sel=length(sp);
                [N]=NURBSBasisFunc(iy+p1,p1,sp,u1,2);
                for ip=1:(p1+1)
                    indn(nel+(1:Sel))=iy+ip-1;
                    indp(nel+(1:Sel))=found;
                    val(nel+(1:Sel))=N(:,ip,1);
                    dval(nel+(1:Sel))=N(:,ip,2);
                    nel=nel+Sel;
                end
                
            end
            indp((nel+1):end)=[];
            indn((nel+1):end)=[];
            val((nel+1):end)=[];
            dval((nel+1):end)=[];
            
            phis=[phis,sparse(indp,indn,val,prod(sizeim),Nny)];
            phit=[phit,sparse(prod(sizeim),Nny)];
            
            phias=[phias,sparse(indp,indn,dsdx*dval,prod(sizeim),Nny)];
            phits=[phits,sparse(prod(sizeim),Nny)];
            phiss=[phiss,sparse(prod(sizeim),Nny)];
        else
            phis=[phis,sparse(1:prod(sizeim),1,1,prod(sizeim),1)];
            phit=[phit,sparse(prod(sizeim),1)];
            
            phias=[phias,sparse(prod(sizeim),1)];
            phits=[phits,sparse(prod(sizeim),1)];
            phiss=[phiss,sparse(prod(sizeim),1)];
            
        end
        
        phixs=t(1)*phis+n(1)*phit;
        phiys=t(2)*phis+n(2)*phit;
        
end
delete([fullfile('VTK',filreso),'-0*.vtk']);

ExportVTKField(filreso,images,xo,yo,elt,conn,{[phixs*U;phiys*U],phias*U,phits*U,phiss*U},{'Displacement','Axial-strain','Transverse-strain','Shear-strain'});

rint=false;
ng=0;
rflag=0;
xo=xo-roi(1)+1;
yo=yo-roi(3)+1;
Up=U;
U=[phixs*Up;phiys*Up];
Eax=phias*Up;
Esh=phiss*Up;
save(filres,'U','Up','s','Eax','Esh','fleche','curv','rot','naxis','Nnodes','Nelems','xo','yo','param','model','nmod','conn','elt','rint','ng','rflag','-v7.3');
switch model.basis
    case 'beam'
        save(filres,'strain','-append');
end
%ExportImageToVTK(filres)
%%
end