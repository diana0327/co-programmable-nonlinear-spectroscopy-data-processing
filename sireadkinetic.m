function [x,y]=sireadkinetic(filename)
    data = sifreadnk(filename);
    close all
    for i=1:data.kineticLength
        Y(:,i) = data.imageData(:,:,i)/2;
        x = data.axisWavelength;
        y=mean(Y,2);
    end
    [maxi,number]=max(max(Y(1900:2030,:)));
    y=Y(:,number);
end 
