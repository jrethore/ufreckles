function buffer=writeAnsysMat(matmod,model)
buffer=[];
        buffer=[buffer,sprintf('MPTEMP,,,,,,,,\n')];
        buffer=[buffer,sprintf('MPTEMP,1,0\n')];
switch matmod
    case 'elastic_homogeneous_isotropic'
        young=model.young;
        nu=model.nu;
        buffer=[buffer,sprintf('MPDATA,EX,1,, %12.5e\n',young)];
        buffer=[buffer,sprintf('MPDATA,PRXY,1,, %12.5e\n',nu)];
    case 'elastic_plastic_homogeneous_isotropic'
        young=model.young;
        nu=model.nu;
        buffer=[buffer,sprintf('MPDATA,EX,1,, %12.5e\n',young)];
        buffer=[buffer,sprintf('MPDATA,PRXY,1,, %12.5e\n',nu)];
        buffer=[buffer,sprintf('TB,BISO,1,1\nTBTEMP,  0.00000000    ,   1\n')];
       buffer=[buffer,sprintf('TBDAT,1,%12.5e,%12.5e\n',model.Sy,model.Etan)];
       if isfield(model,'C1')&&isfield(model,'C2')&&isfield(model,'C3')
        buffer=[buffer,sprintf('TB,CHAB,1,1,1,3\nTBTEMP,  0.00000000    ,   1\n')];
       buffer=[buffer,sprintf('TBDAT,1,%12.5e,%12.5e,%12.5e\n',model.C1,model.C2,model.C3)];
       end
end



end