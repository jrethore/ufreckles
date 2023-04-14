function buffer=writeMooseMat(matmod,model)
buffer=[];
switch matmod
    case 'elastic_homogeneous_isotropic'
        young=model.young;
        nu=model.nu;
        buffer=[buffer,sprintf('matlaw = 1\n')];
        buffer=[buffer,sprintf('yngmod = %12.5e\n',young)];
        buffer=[buffer,sprintf('poisso = %12.5e\n',nu)];
        buffer=[buffer,sprintf('densit = %12.5e\n',1)];
end


end