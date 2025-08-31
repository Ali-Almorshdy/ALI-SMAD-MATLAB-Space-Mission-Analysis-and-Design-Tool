function [rmag,ra,dec,vmag,rav,decv] = C2SP(r,v)
                    % ––––––––––––––––––––––––––––––––––––––––––––––
                    
                    dec = atan2(r(3),sqrt(r(1)^2+r(2)^2))*180/pi;
                    ra=atan2(r(2),r(1))*180/pi;
                    decv=atan2(v(3),sqrt(v(1)^2+v(2)^2))*180/pi;
                    rav=atan2(v(2),v(1))*180/pi;
                    rmag=norm(r);
                    vmag=norm(v);
          end