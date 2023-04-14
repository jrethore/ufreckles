clear all

[fil2,path1]=uigetfile('*.tif','Select a tiff image');
eval(['cd ' path1]);
        
        
       tmp=imread(fil2);
       answer=inputdlg({'X dim','Y dim' 'Z dim?',},...
            'Number of images',1,{num2str(size(tmp,1),'%3d'),num2str(size(tmp,2),'%3d'),'200'});
nbi=str2num(answer{3});
       jm3=zeros(size(tmp,1),size(tmp,2),nbi);
        for jj=1:nbi
        tmp=imread(fil2,0+jj);
        jm3(:,:,jj)=tmp;
        end
        imM=max((jm3(:)));
imm=min((jm3(:)));
immoy=(imM+imm)/2;
jm3=0.5+(jm3-immoy)/(imM-imm)*0.95; % uses 95% of the interval [0:1] centered on .5
jm3=(uint8(jm3*256));

        fil2=strrep(fil2,'.tif','.mat');
%        jm3=permute(jm3,[3 1 2]);
        save(fil2,'jm3');
