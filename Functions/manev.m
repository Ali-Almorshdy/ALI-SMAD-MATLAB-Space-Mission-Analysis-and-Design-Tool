
    function manev(x,yy,z,s)
        
        S=1;
        mu = 3.986012E5;
                            dV=x;
                       dV1=yy;
                       dV2=z;
                        file='Files\end.txt';
%             file1='C:\Users\Ali\Desktop\pass.orbit';
            f=fileread(file);
        X=str2double( regexp(f,'(?<=X=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            Y=str2double( regexp(f,'(?<=Y=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            Z=str2double( regexp(f,'(?<=Z=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            VX=str2double( regexp(f,'(?<=V1=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            VY=str2double( regexp(f,'(?<=V2=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            VZ=str2double( regexp(f,'(?<=V3=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            t0=str2double( regexp(f,'(?<=t0=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            R00=[X Y Z]';
            V00=[VX VY VZ]';
            norm(V00)
            coe=ctok(R00,V00);
            e = coe(2)
                    RA = coe(3)*180/pi;
                    incl = coe(4)*180/pi;
                    w = coe(5)*180/pi;
                    TA0 = coe(6)*180/pi
                    a = coe(7);
                    
                     Pp=a*(1-e*e);
            mu = 3.986012E5;
            h=sqrt(Pp*mu);
            rp = (h^2/mu) * (1/(1 + e*cos(TA0*pi/180))) * (cos(TA0*pi/180)*[1;0;0] + sin(TA0*pi/180)*[0;1;0]);
            vp = (mu/h) * (-sin(TA0*pi/180)*[1;0;0] + (e + cos(TA0*pi/180))*[0;1;0])
%             (mu/h) * -sin(TA0*pi/180)*[1;0;0]
            nor=cross([vp(1) 0 0],[0 vp(2) 0])
            if vp(1)>=0
            vp(1)=vp(1)+dV1;
            else
             
             
             vp(1)=vp(1)-dV1;
            end
            k=vp(2);
             l=vp(1);
             j=dV;
            ince=dV+((2*j*(k^2 + l^2)^(1/2) + j^2 + k^2)^(1/2) - k - j);
            ince2=(- j - k - (2*j*(k^2 + l^2)^(1/2) + j^2 + k^2)^(1/2))+dV;
%             v3=dV2+( dV2 - (dV2*(dV2 + 2*(k^2 + l^2)^(1/2)))^(1/2))
%             v4=(dV2*( + 2*(k^2 + l^2)^(1/2)))^(1/2) - j
            if abs((ince)-(dV))<0.5
                vp(2)=vp(2)+ince;
%                 ince
            else
                vp(2)=vp(2)+ince2;
%                 ince2
            end
%        if vp(2)>=0
%            vp(2)=vp(2)+dV;
%        else
%            vp(2)=vp(2)-dV;
%        end
if vp(3)>=0
            vp(3)=vp(3)+dV2;
            else
             
             
             vp(3)=vp(3)-dV2;
            end
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
            V0=v';
            R0=r';
            X0=[R00;V0];
            norm(V0)
%                         R0=[R00(1) R00(2) R00(3)]';
%                         V0=[V00(1) V00(2) V00(3)]';
                      
%                         if V00(2)<dV
%                             V00(2)=V00(2)-dV;
%                         else
%                             V00(2)=V00(2)+dV;
%                         end
%                         V00(1)=V00(1)+dV1;
%                         V0(3)=V0(3)+dV2;
%                          if V00(1)<dV1
%                             V00(1)=V00(2)-dV1;
%                         else
%                             V00(1)=V00(2)+dV1;
%                          end
%                          if V00(3)<dV2
%                             V00(3)=V00(3)-dV2;
%                         else
%                             V00(3)=V00(3)+dV2;
%                         end
%                         X0=[R0;V0];
%                         EoM=@dfdf;

 coe2=ctok(R00,V0);
            e = coe2(2);
                    
                    incl = coe2(4)*180/pi
                    w = coe2(5)*180/pi;
                    TAo = coe2(6);
%                     a = coe2(7)
                        E = norm(V0)^2/2 - mu/norm(R00);
                        aa = -mu/(2*E);
                        n = sqrt(abs(mu/abs(aa^3)));
            T = abs(2*pi/n);
            n_periods=S;
            Eo = 2*atan(tan(TAo/2)*sqrt((1-e)/(1+e)));
            Mo = Eo - e*sin(Eo);
            to = Mo*(T/2/pi)
            tf =  real(to+n_periods*T*s);
            times = [to:10:tf];
                        % set(handles.text28, 'String', T);
                      options = odeset('reltol', 1.e-8, ...
                'abstol', 1.e-8, ...
                'initialstep', T/10000);
            
                        [time, State] = ode45(@(t,y) dfdf(t,y,t0), times, X0, options);
            
                        xv = real(State(:,1));
                        yv = real(State(:,2));
                        zv = real(State(:,3));
                        v1 = real(State(:,4));
                        v2 = real(State(:,5));
                        v3 = real(State(:,6));
                        delete(findobj('tag','S1')); delete(findobj('tag','S3')); delete(findobj('tag','S5'));
                         delete(findobj('tag','S2')); delete(findobj('tag','S4')); delete(findobj('tag','S6'));
              
                        [app.sat1,app.sat2,app.sat3,app.sat4,app.sat5,app.sat6] =AliCube(a/10,real(xv(end)),real(yv(end)),real(zv(end)));
                        plot3(xv,yv,zv,'.');
                        R000=[xv(end) yv(end) zv(end)];
                        V000=[v1(end) v2(end) v3(end)];
                        
                         coe=ctok(R000,V000);
            
                    incl = coe(4)*180/pi
                      e = coe2(2)
                    a = coe(7)
                    endfile='Files\end.txt';
            eend=fopen(endfile,'w');
            fprintf(eend,'X=%f\n',xv(end));
            fprintf(eend,'Y=%f\n',yv(end));
            fprintf(eend,'Z=%f\n',zv(end));
            fprintf(eend,'V1=%f\n',v1(end));
            fprintf(eend,'V2=%f\n',v2(end));
            fprintf(eend,'V3=%f\n',v3(end));
            tt0=(10/86400)*(length(xv))+t0;
            fprintf(eend,'t0=%f\n',tt0);
%         delete(findobj('tag','MV'),'string')
% delete(findobj('tag','MN'))
% delete(findobj('tag','MB'))
% delete(findobj('tag','MV'))
% delete(findobj('tag','bu1'))
%  delete(findobj('tag','V'))
% delete(findobj('tag','N'))
%  delete(findobj('tag','B'))
% delete(findobj('tag','S'))
%  delete(findobj('tag','MS'))
 
%  file1='G:result.txt';
% %  if(clc==1)
%          fr = fopen(file1,'w');
% %  else
%      fr = fopen(file1,'a');
% %  end
%          fprintf(fr,append('smaneuver',num2str(clc)));
%          fprintf(fr,'\n');
%          fprintf(fr,'1-State Vectot:\n');
%             fprintf(fr,'__________________________________________\n');
%             fprintf(fr,' 1-X= %f km \n',xv(end));
%             fprintf(fr,' 2-Y= %f km  \n',yv(end));
%             fprintf(fr,' 3-Z= %f km \n',zv(end));
%             fprintf(fr,' 4-VX= %f km/S \n',v1(end));
%             fprintf(fr,' 5-VY= %f km/S \n',v2(end));
%             fprintf(fr,' 6-VZ= %f km/S \n\n',v3(end));
%             
%             
%             coe=ctok([xv(end) yv(end) zv(end)]',[v1(end) v2(end) v3(end)]');
%                     e2 = coe(2);
%                     RA = coe(3)*180/pi;
%                     incl2 = coe(4)*180/pi;
%                     w = coe(5)*180/pi;
%                     TA22 = coe(6)*180/pi;
%                     a2 = coe(7);
%             EA=2*atan(sqrt((1-e2)/(1+e2))*tan((TA22*(pi/180))/2));
%             MA=EA-e2*sin(EA);
%             fprintf(fr,'1-Keplerian Elements:\n');
%             fprintf(fr,'__________________________________________\n');
%             fprintf(fr,' 1-Semi Majoraxis = %f km \n',a2);
%             fprintf(fr,' 2-eccentricity = %f  \n',e2);
%             fprintf(fr,' 3-Inclination = %f° \n',incl2);
%             fprintf(fr,' 4-Rightessention = %f° \n',RA);
%             fprintf(fr,' 5-Argument Of Perigee = %f° \n',w);
%             fprintf(fr,' 6-Mean Anomaly = %f° \n',MA*180/pi);
%             fprintf(fr,' 7-Eccentric Anomaly = %f° \n',EA*180/pi);
%             fprintf(fr,' 8-True Anomaly = %f° \n\n',TA22);
%             
%             
%              fprintf(fr,'2-Spherical State :\n');
%         [rmag,ra,dec,vmag,rav,decv] = C2SP([xv(end) yv(end) zv(end)],[v1(end) v2(end) v3(end)]);
% %             PP=a*(1-e*e);
%             fprintf(fr,'__________________________________________\n');
%             fprintf(fr,'RMAG = %f km \n',rmag);
%             fprintf(fr,'Right Ascension = %f deg \n',ra);
%             fprintf(fr,'Declination = %f km^2/S^2 \n',dec);
%             fprintf(fr,'VMAG = %f KM/S \n',vmag);
%             fprintf(fr,'Vel Right Ascension = %f deg \n',rav);
%             fprintf(fr,'Vel Declination  = %f deg \n\n',decv);
%             
%             
%             fprintf(fr,'2-Other Orbit Data:\n');
%             E = norm([v1(end) v2(end) v3(end)])^2/2 - mu/norm([xv(end) yv(end) zv(end)]);
%             aa = -mu/(2*E);
%             n = sqrt(abs(mu/abs(aa^3)));
%             T = abs(2*pi/n);
% %             PP=a*(1-e*e);
%             rp=a2*(1-e2);
%             hv = cross([xv(end) yv(end) zv(end)], [v1(end) v2(end) v3(end)]);
%             hmag = sqrt(dot(hv, hv));
%             vp=hmag/rp;
%             ra=a2*(1+e2);
%             va=hmag/ra;
%             p=sqrt(rp*ra);
%             mo=(2*pi)/T;
%             fprintf(fr,'__________________________________________\n');
%             fprintf(fr,'Semi minor Axis = %f km \n',p);
%             fprintf(fr,'specific angular momentum = %f KM^2/S \n',hmag);
%             fprintf(fr,'Orbit Energy = %f km^2/S^2 \n',E);
%             fprintf(fr,'Mean motion = %f RAD/S \n',mo);
%             fprintf(fr,'Periapsis = %f km \n',rp);
%             fprintf(fr,'VelPeriapsis  = %f km/S \n',vp);
%             fprintf(fr,'Apoapsis = %f km \n',ra);
%             fprintf(fr,'VelApoapsis = %f km/S \n',va);
%             fprintf(fr,'Orbit Period = %f S \n',T);
%             fprintf(fr,append('emaneuver',num2str(clc)));
%             fprintf(fr,'\n');

         
    end

      

