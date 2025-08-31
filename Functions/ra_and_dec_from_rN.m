       function [ra, dec] = ra_and_dec_from_rN(r) 
        dec = atan2(r(:,3),sqrt(r(:,1).^2+r(:,2).^2))*180/pi;
            ra=atan2(r(:,2),r(:,1))*180/pi;
       end