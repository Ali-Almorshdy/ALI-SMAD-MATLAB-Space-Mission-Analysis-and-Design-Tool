function [R0, V0] = COE2RV(a,e,TA,RA,incl,w)
            Pp=a*(1-e*e);
            mu = 3.986012E5;
            h=sqrt(Pp*mu);
            rp = (h^2/mu) * (1/(1 + e*cos(TA*pi/180))) * (cos(TA*pi/180)*[1;0;0] + sin(TA*pi/180)*[0;1;0]);
            vp = (mu/h) * (-sin(TA*pi/180)*[1;0;0] + (e + cos(TA*pi/180))*[0;1;0]);
            %% Rotation matrices
            R3_W = [ cos(RA*pi/180) sin(RA*pi/180) 0
                -sin(RA*pi/180) cos(RA*pi/180) 0
                0        0       1];
            R1_i = [1 0         0
                0 cos(incl*pi/180) sin(incl*pi/180)
                0 -sin(incl*pi/180) cos(incl*pi/180)];
            R3_w = [ cos(w*pi/180) sin(w*pi/180) 0
                -sin(w*pi/180) cos(w*pi/180) 0
                0       0      1];
            Q_pX = (R3_w*R1_i*R3_W)';
            r = Q_pX*rp;
            v = Q_pX*vp;
            r = r';
            v = v';
            V0=v;
            R0=r;
            X0=[R0;V0];
        end