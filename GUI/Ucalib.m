function filc=Ucalib(path0,fil0)

cd(path0);
[filc,pathres]=uiputfile({'*.cal','Ufreckles calibration file (*.cal)'},'Save as a new calibration data set');

fils=SortImageFiles(fil0);


nmod=3;

param.result_file=filc;
param.analysis='calibration';
param.deformed_image=fils;

                    target=inputdlg({'Grid step [m]','Steps along X','Steps along Y'},...
                        'Target parameters',1,{'0.001','9','9'});

model.grid_step=eval(target{1});
model.grid_size=[eval(target{2}),eval(target{3})];
LoadParameters(param);
LoadParameters(model,nmod);
%%
[T,residus]=CalibrationTarget(nmod);

 save(filc,'T','model','param','nmod','residus');

end