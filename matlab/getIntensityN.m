
function [ intensity, overload ] = getIntensity( cobracenter, min_range, max_range, imageMatrix, resolution_angle, imgdir, nDir, numBuckets, ontime, pid )

figure (15)
caxis([0 20000])

title(strcat('Ontime: ', ontime, ' pId: ', pid));

intensity=zeros(1,round(2*pi/resolution_angle)+1);
ymin = uint32(imag(cobracenter)- max_range);
xmin = uint32(real(cobracenter) - max_range);

% reg.x = (-max_range:max_range) + round(real(cobracenter));
% reg.y = (-max_range:max_range) + round(imag(cobracenter));

bf = 3; %Boost FActor by which the resolution is increased. 
subMatrix = imageMatrix(ymin:uint32(floor(ymin)+max_range*2),xmin:uint32(floor(xmin)+max_range*2));
boostedMatrix = zeros((max_range*2)*bf,(max_range*2)*bf);

for rr = 1:length(subMatrix(1,:));
    for cc = 1:length(subMatrix(:,1));
    boostedMatrix(bf*(cc-1)+1:bf*cc,bf*(rr-1)+1:bf*rr) = double(subMatrix(cc,rr))./bf^2;
    end
end
 
boostedCenter = max_range*bf + 1i * max_range*bf;

numPixels = 0;
numOverloadedPixels = 0;
for rr = 1:length(boostedMatrix(1,:))
    for cc = 1:length(boostedMatrix(:,1))
        cur_vector = rr+1i*cc-boostedCenter;
        cur_dist=abs(cur_vector);
       % cur_angle=mod(angle(cur_vector),2*pi);
        cur_angle=angle(cur_vector);
        if cur_dist<=max_range*bf
           if cur_dist>=min_range*bf 
              %  bin=round(cur_angle/resolution_angle)+1;
              bin=round(cur_angle/resolution_angle)+numBuckets/2+1;
        
                %  Real Code (and not Calibration)
                % Filter out the spillage pixels
               % if(boostedMatrix(cc, rr) < 4*10^4/bf^2)
                intensity(bin)=intensity(bin)+boostedMatrix(cc, rr);
                J(cc, rr) = boostedMatrix (cc, rr);
                %else 
                 %   spilledPixelValue = 5;
                 %   J(cc, rr) = spilledPixelValue;
                %end
%               % Calibration Code
%               intensity(bin)=intensity(bin)+(mod(bin,2)*9000+2000)/bf^2;
%               J(cc, rr) = (mod(bin,2)*9000+2000)/bf^2;
%             intensity(bin)=intensity(bin)+10000/bf^2;
%             J(cc, rr) = 10000/bf^2;
                
                if(boostedMatrix(cc, rr)> 1000/bf^2)
                      numPixels = numPixels + 1;
                end
                if(boostedMatrix(cc, rr)> 40000/bf^2)
                    numOverloadedPixels = numOverloadedPixels + 1;
                end
           else
               J(cc,rr) = 1.5*10^3/bf^2;% To show the inner circle  
           end
       else
        J(cc,rr) = 10^3/bf^2; % To show the outer circle
        end
    end
end

overload = numOverloadedPixels/(numPixels/100);
intensity(1) = intensity(1)+ intensity(end);
intensity(end:end) = [];
J(uint16(imag(boostedCenter)-4:imag(boostedCenter)+4),uint16(real(boostedCenter)-4:real(boostedCenter)+4)) = 4*10^3/bf^2;
image(J,'CDataMapping','scaled')
 title(strcat('Ontime: ', num2str(ontime), ' µs. ', ' pId: ', num2str(pid)));

% figure(999) % Check image for consistency of boost method.
% image(subMatrix, 'CDataMapping', 'scaled');
%  saveas(gcf,horzcat(imgdir,'image',num2str(nDir),'.fig'));
%  saveas(gcf,horzcat(imgdir,'image',num2str(nDir),'.png'));

end

