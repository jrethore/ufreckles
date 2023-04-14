function buffer=writeAbaqusMat(matmod,model)
buffer=[];
switch matmod
    case 'aline'
        buffer=[buffer,sprintf('*DENSITY\n%d,\n',model.rho)];
        buffer=[buffer,sprintf('*INCLUDE, INPUT=aline.dat\n')];
        
        fwid = fopen('aline.dat','w');
        count = fprintf(fwid,'*DEPVAR\n36,\n');
        count = fprintf(fwid,'*USER MATERIAL, CONSTANTS=4\n');
        count = fprintf(fwid,'%d, %d, %d, %d\n',model.young,model.nu,model.Yf,model.Er);
        fclose(fwid);
    case 'elastic_homogeneous_isotropic'
        young=model.young;
        nu=model.nu;
        buffer=[buffer,sprintf('*ELASTIC,TYPE=ISOTROPIC\n')];
        buffer=[buffer,sprintf('%12.5e, %12.5e\n',young,nu)];
    case 'hyperelastic_homogeneous_isotropic'
        switch model.type
            case 'mooney-rivlin'
                buffer=[buffer,sprintf('*HYPERELASTIC,MOONEY-RIVLIN\n')];
                 buffer=[buffer,sprintf('%12.5e, %12.5e, %12.5e\n',model.C10,model.C01,model.D1)];
            case 'arruda-boyce'
                buffer=[buffer,sprintf('*HYPERELASTIC,ARRUDA-BOYCE\n')];
                 buffer=[buffer,sprintf('%12.5e, %12.5e, %12.5e\n',model.mu,model.lm,model.D)];
        end
    case 'elastic_plastic_homogeneous_isotropic'

        buffer=[buffer,sprintf('*ELASTIC,TYPE=ISOTROPIC\n')];
        buffer=[buffer,sprintf('%12.5e, %12.5e\n',model.young,model.nu)];
        switch model.hardening
            case 'perfect'
                buffer=[buffer,sprintf('*PLASTIC\n')];
                buffer=[buffer,sprintf('%12.5e\n',model.Sy)];
            case 'lin'
                buffer=[buffer,sprintf('*PLASTIC, HARDENING=KINEMATIC\n')];
%                buffer=[buffer,sprintf('*PLASTIC, HARDENING=ISOTROPIC\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e\n',model.Sy,0)];
                buffer=[buffer,sprintf('%12.5e, %12.5e\n',model.Sy+model.H,1)];
            case 'power'
                Ep=(0:1e-2:10^(1/3)).^3;
                Sy=model.Sy+(model.C)*Ep.^(model.delta);
                buffer=[buffer,sprintf('*PLASTIC\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e\n',[Sy;Ep])];
            case 'swift'
                Ep=(0:1e-2:10^(1/3)).^3;
                Sy=model.B*(model.C+Ep).^(model.delta);
                if isfield(model,'H')
                    a=1;
                    if isfield(model,'a')
                        a=model.a;
                    end
                    Sy=Sy+a*model.H*Ep;
                end
                 buffer=[buffer,sprintf('*PLASTIC\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e\n',[Sy;Ep])];                
            case {'isotropic','exp'}
                buffer=[buffer,sprintf('*PLASTIC,HARDENING=COMBINED,DATA TYPE=PARAMETERS\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e, %12.5e\n',model.Sy,0,0)];
                buffer=[buffer,sprintf('*CYCLIC HARDENING,PARAMETERS\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e, %12.5e\n',model.Sy,model.Q,model.b)];
            case 'kinematic'
                buffer=[buffer,sprintf('*PLASTIC,HARDENING=COMBINED,DATA TYPE=PARAMETERS\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e, %12.5e\n',model.Sy,model.C,model.gam)];
                buffer=[buffer,sprintf('*CYCLIC HARDENING,PARAMETERS\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e, %12.5e\n',model.Sy,0,0)];
            case 'combined'
                buffer=[buffer,sprintf('*PLASTIC,HARDENING=COMBINED,DATA TYPE=PARAMETERS\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e, %12.5e\n',model.Sy,model.C,model.gam)];
                buffer=[buffer,sprintf('*CYCLIC HARDENING,PARAMETERS\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e, %12.5e\n',model.Sy,model.Q,model.b)];
            case 'Johnson-Cook'
                 buffer=[buffer,sprintf('*PLASTIC,HARDENING=JOHNSON COOK\n')];
                buffer=[buffer,sprintf('%12.5e, %12.5e, %12.5e\n',model.Sy,model.H,model.n)];
        end
      if isfield(model,'ultimate_plastic_strain')
                  buffer=[buffer,sprintf('*DAMAGE INITIATION, CRITERION=%s\n',model.criterion)];
                 buffer=[buffer,sprintf('%12.5e, %12.5e\n',model.ultimate_plastic_strain')];
                  buffer=[buffer,sprintf('*DAMAGE EVOLUTION, TYPE=ENERGY, SOFTENING=%s\n',model.softening)];
                 buffer=[buffer,sprintf('%12.5e\n',model.fracture_energy)];

      end
end


end