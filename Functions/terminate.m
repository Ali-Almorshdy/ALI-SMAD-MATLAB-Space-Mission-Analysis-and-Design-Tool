 function [lookfor, stop, direction] = terminate(t,y)
                %
                % This function specifies the event at which ode45 terminates.
                % ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
                ae = 6378;             % [m]
                r1 = norm(y(1:3));
                alt = r1 - ae;
                lookfor = alt; % = 0 when altitude = 100 km
                stop = 1; % 1 means terminate at lookfor = 0; Otherwise 0
                direction = -1; % -1 means zero crossing is from above
            end %terminate