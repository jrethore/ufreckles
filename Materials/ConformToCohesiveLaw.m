function [U]=ConformToCohesiveLaw(Uo,uno,tno,upo,tpo,nmodel)
load(fullfile('TMP','params'));
pix2m=param.pixel_size;
load(fullfile('TMP',sprintf('%d_params',nmodel)),'param');
nsub=length(param.sub_indices);
load(fullfile('TMP',sprintf('%d_dpsin_0',nmodel)),'dpsincz','indkcz');
load(fullfile('TMP',sprintf('%d_dphin_0',nmodel)),'dphincz');
dpsincz=dpsincz(:,indkcz);
dphincz=dphincz(:,indkcz)*pix2m;
check=1;
Ko=Uo(indkcz);
dtnodx=gradient(tno);
dunodx=gradient(uno);
dtnoduno=dtnodx./dunodx;
if check
    figure1 = figure('XVisual',...
        '0x23 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
        'PaperSize',[20.98 29.68]);
    axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16,...
        'FontName','Times');
    box('on');
    hold('all');
    plot(uno(:),tno(:),'b-','DisplayName','Prescribed','Parent',axes1,'LineWidth',2);
    title('Cohesive law','FontSize',20,'FontName','Times');
    xlabel('[U]','FontSize',20,'FontName','Times');
    ylabel('T','FontSize',20,'FontName','Times');
end

res=1.;
itera=0;
K=Ko;
alpha=1;
while (res>1.e-6)&&(itera<10/alpha)
    itera=itera+1;
    un=dphincz*K;
    tn=dpsincz*K;
    up=upo*norm(K)/norm(Ko);
    tp=tpo*norm(K)/norm(Ko);
    [dtdu,tu]=ComputeCohesiveTractionAndTangentStiffness(uno,tno,dtnoduno,un+up,tn+tp);
    if check
        if itera==1
            plot(un+up,tn+tp,'k-','DisplayName','Initial','Parent',axes1,'LineWidth',2,'MarkerSize',5);
        else
        h1=  plot(un+up,tn+tp,'r-','DisplayName','Corrected','Parent',axes1,'LineWidth',2,'MarkerSize',5);
        hl=legend(axes1,'show','Location','NorthEast');
        pause(0.5)
        delete(h1);
        delete(hl);
        end
    end
    R=dpsincz*K+tp-tu;
    L=dpsincz-diag(dtdu)*dphincz;

    M=L'*L;
    F=-L'*R;

    if itera==1
        normRo=norm(R);
        normKo=norm(K);
    end

    dK=M\F;

    K=K+alpha*dK;
    res=norm(dK)/normKo;
%     if check
%         disp(sprintf('|dK|/Ko=%f',res));
%         disp(sprintf('|F|/Fo=%f',norm(R)/normRo));
% 
%     end
end
U=Uo;
U(indkcz)=K;
disp(sprintf('After %d iterations: |dK|/Ko=%f and |F|/Fo=%f',itera,res,norm(R)/normRo));
if check
    h1=  plot(un+up,tn+tp,'r-','DisplayName','Corrected','Parent',axes1,'LineWidth',2,'MarkerSize',5);
    hl=legend(axes1,'show','Location','NorthEast');
end

end




