function [siz] = tiffdims(filename)                                                                               
                                                                                                                                     
    tmp = imread(filename);                                                                                                          
%    tmp=(permute(tmp,[2,1,3]));
%tmp=fliplr(tmp);
siz = size(tmp);                                                                                                       
                                                                                                                                     
    ii = 100;                                                                                                                         
    imin = 0;                                                                                                                        
    imax = -2;                                                                                                                       
                                                                                                                                     
    while abs(imax - imin) > 1                                                                                                       
        try                                                                                                                          
            imread(filename, ii);                                                                                                     
            imin = ii;                                                                                                                
        catch me                                                                                                                     
            imax = ii;                                                                                                                
        end                                                                                                                          
        if imax > 0                                                                                                                  
            ii = floor((imax + imin) / 2);                                                                                            
        else                                                                                                                         
            ii = 2 * ii;
        end
    end
    
    siz = [siz,ii];
    
end
