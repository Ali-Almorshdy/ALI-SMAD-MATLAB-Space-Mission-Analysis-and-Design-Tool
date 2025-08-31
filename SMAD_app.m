classdef SMAD_app < matlab.apps.AppBase
   
    % Properties that correspond to app components
    properties (Access = public)
        MainGUI                matlab.ui.Figure
        MissionMainButtonsBar  matlab.ui.container.Toolbar
        WriteScript            matlab.ui.container.toolbar.PushTool
        OpenScript             matlab.ui.container.toolbar.PushTool
        SaveChanges            matlab.ui.container.toolbar.PushTool
        Restart                matlab.ui.container.toolbar.PushTool
        RunButoon              matlab.ui.container.toolbar.PushTool
        AnimationButton        matlab.ui.container.toolbar.PushTool
        PauseButton            matlab.ui.container.toolbar.ToggleTool
        FasterAnimation        matlab.ui.container.toolbar.PushTool
        SlowerAnimation        matlab.ui.container.toolbar.PushTool
        OrbitMaximun           matlab.ui.container.toolbar.PushTool
        OrbitMinum             matlab.ui.container.toolbar.PushTool
        GroundMaximum          matlab.ui.container.toolbar.PushTool
        GroundMinum            matlab.ui.container.toolbar.PushTool
        PushHOLDTool4          matlab.ui.container.toolbar.PushTool
        GUIGrild               matlab.ui.container.GridLayout
        Panel                  matlab.ui.container.Panel
        PanelGrild             matlab.ui.container.GridLayout
        Image                  matlab.ui.control.Image
        GroundPanel            matlab.ui.container.Panel
        OrbitPanel             matlab.ui.container.Panel
        MissionButtoms         matlab.ui.container.TabGroup
        orbitTab               matlab.ui.container.Tab
        OrbitTree              matlab.ui.container.Tree
        SpaceCraftNode         matlab.ui.container.TreeNode
        BurnsNode              matlab.ui.container.TreeNode
        GroundTrackNode        matlab.ui.container.TreeNode
        GroundStationNode      matlab.ui.container.TreeNode
        PlotNode               matlab.ui.container.TreeNode
        STKPlotNode            matlab.ui.container.TreeNode
        EditNode               matlab.ui.container.TreeNode
        Node                   matlab.ui.container.TreeNode
        missionsTab            matlab.ui.container.Tab
        MissionTreee           matlab.ui.container.Tree
        OrbitNode              matlab.ui.container.TreeNode
        OutputTab              matlab.ui.container.Tab
        OutputTree             matlab.ui.container.Tree
        OrbitViewNode          matlab.ui.container.TreeNode
        GroundTrackNode_2      matlab.ui.container.TreeNode
        TargetsMenu            matlab.ui.container.ContextMenu
        addmaneuverMenu        matlab.ui.container.Menu
        atargetMenu            matlab.ui.container.Menu
        etargetMenu            matlab.ui.container.Menu
        EtargetMenu            matlab.ui.container.Menu
        itargetMenu            matlab.ui.container.Menu
    end



    properties (Access = public)

        count ; %% num of Manauvers;
        c ; %% for animation Restart
        RAS; %%Orbit RAAD
        DECL ; %% For Orbit Decl
        calc=1 ; %% animation speed
        x=1; %% Ground animation
        I;%% Grond Trace Image
        h; %% OrbitFigure
        f2; %% Ground Trace Figure
        ano=0; %%Orbit Or ground
        XV;
        YV;
        ZV;
        SMA;
        sat1;
        sat2;
        sat3;
        sat4;
        sat5;
        sat6;
        increase=10;
        earh;
        check=0;
        start;
        hgx;
        END;
        inc;
        tt;
        text1;
        ABS=0;
        number;
        constrain=0;
        Text;
        num=0;
        sun;
        ta;
        v1;
        v2;
        v3;
        man;
    end

    properties (Access = private)

    end

    methods (Access = public)



        function [ra, dec] = ra_and_dec_from_r(~,r)
            % ––––––––––––––––––––––––––––––––––––––––––––––
            dec = atan2(r(3),sqrt(r(1)^2+r(2)^2))*180/pi;
            ra=atan2(r(2),r(1))*180/pi;
        end

        function E = kepler_E(~,e, M)
            error = 1.e-8;
            %...Select a starting value for E:
            if M < pi
                E = M + e/2;
            else
                E = M - e/2;
            end

            ratio = 1;
            while abs(ratio) > error
                ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
                E=E-ratio;
            end
        end

        function [R0, V0] = KtoC(~,a,e,TA,RA,incl,w)
            Pp=a*(1-e*e);
            mu = 3.986012E5;
            hM=sqrt(Pp*mu);
            rp = (hM^2/mu) * (1/(1 + e*cos(TA*pi/180))) * (cos(TA*pi/180)*[1;0;0] + sin(TA*pi/180)*[0;1;0]);
            vp = (mu/hM) * (-sin(TA*pi/180)*[1;0;0] + (e + cos(TA*pi/180))*[0;1;0]);
            %             vp(3)=vp(3)+;

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
        end


        function gmstime = gmst(~,Mjd_UT1)

            Secs = 86400;                       % Seconds per day
            MJD_J2000 = 51544.5;

            Mjd_0 = floor(Mjd_UT1);
            UT1   = Secs*(Mjd_UT1-Mjd_0);       % [s]
            T_0   = (Mjd_0  -MJD_J2000)/36525;
            T     = (Mjd_UT1-MJD_J2000)/36525;

            gmst  = 24110.54841 + 8640184.812866.*T_0 + 1.002737909350795.*UT1...
                + (0.093104-6.2e-6.*T).*T.*T;  % [s]

            % [rad], 0..2pi
            gmstime = 2*pi*Frac2(gmst/Secs);
        end



        function [rmag,ra,dec,vmag,rav,decv] = C2SP(~,r,v)
            % ––––––––––––––––––––––––––––––––––––––––––––––

            dec = atan2(r(3),sqrt(r(1)^2+r(2)^2))*180/pi;
            ra=atan2(r(2),r(1))*180/pi;
            decv=atan2(v(3),sqrt(v(1)^2+v(2)^2))*180/pi;
            rav=atan2(v(2),v(1))*180/pi;
            rmag=norm(r);
            vmag=norm(v);
        end


    end

    methods (Access = public)



    end



    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
             
            app.check=0;
            app.count=1;
            app.c=0;
            app.ABS=0;
            app.num=0;
            app.constrain=0;
            set(app.AnimationButton,'Enable','off');
            set(app.OrbitMaximun,'Enable','off');
            set(app.OrbitMinum,'Enable','off');
            set(app.GroundMaximum,'Enable','off');
            set(app.GroundMinum,'Enable','off');
            set(app.FasterAnimation,'Enable','off');
            set(app.SlowerAnimation,'Enable','off');


        end

        % Callback function: RunButoon
        function RunButoonClicked(app, event)
            set(app.RunButoon,'Enable','off');
            close all;
            try
                delete(app.OrbitNode.Children)

            catch

            end
            app.ABS=0;
            ccc=app.check;
            app.ano=0;
            app.c=3;
            app.x=1;
            app.increase=20;
            app.count=1;
            file='Files\ali.orbit';
            f=fileread(file);
            if(ccc==1)
                file1='Files\pass.orbit';
                f1=fileread(file1);
                lat=str2double( regexp(f1,'(?<=lat=)([+-]?([0-9]*[.])?[0-9]+)','match'));
                long=str2double( regexp(f1,'(?<=long=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            else
                ccc=0;
            end
            a=str2double( regexp(f,'(?<=a=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            e=str2double( regexp(f,'(?<=e=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            incl=str2double( regexp(f,'(?<=i=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            w=str2double( regexp(f,'(?<=w=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            TA=str2double( regexp(f,'(?<=TA=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            RA=str2double( regexp(f,'(?<=RA=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            stop=str2double( regexp(f,'(?<=s=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            y=str2double( regexp(f,'(?<=y=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            m=str2double( regexp(f,'(?<=m=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            d=str2double( regexp(f,'(?<=d=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            hu=str2double( regexp(f,'(?<=h=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            mi=str2double( regexp(f,'(?<=mn=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            sec=str2double( regexp(f,'(?<=sec=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            t0=juliandate(y,m,d,hu,mi,sec);
            t000=juliandate(y,m,d,hu,mi,sec);
            t00=datetime(y,m,d,hu,mi,sec);
            set(app.AnimationButton,'Enable','on');
            set(app.OrbitMaximun,'Enable','on');
            set(app.OrbitMinum,'Enable','on');
            set(app.GroundMaximum,'Enable','on');
            set(app.GroundMinum,'Enable','on');
            set(app.FasterAnimation,'Enable','on');
            set(app.SlowerAnimation,'Enable','on');
            app.ta=t000;
            space_color = 'k';
            app.SMA=a;
            app.number=stop;
            app.h = figure('Color', space_color,'tag','orbits');

            l=app.OrbitPanel.Position;             k1=app.GroundPanel.Position;
            set(app.h, 'position', [l(1)+150 l(2)+30 l(3)+150 l(4)-4],'Toolbar','none', ...
                'menubar','none','name','OrbitView','NumberTitle','off');
            xx=axes;
            ha = axes('units','normalized', ...
                'position',[0 0 1 1]);
            uistack(ha,'bottom');

            Ii=imread('sky.jpg');
            imagesc(Ii);
            colormap gray
            set(ha,'handlevisibility','off', ...
                'visible','off');
            cla;
            mu = 3.986012E5;
            [R0,V0]=KtoC(app,a,e,TA,RA,incl,w);
            X0=[R0;V0];
            E = norm(V0)^2/2 - mu/norm(R0);
            aa = -mu/(2*E);
            n = sqrt(abs(mu/abs(a^3)));
            T = abs(2*pi/n);
            n_periods=stop;
            TAo = TA*pi/180;
            Eo = 2*atan(tan(TAo/2)*sqrt((1-e)/(1+e)));
            Mo = Eo - e*sin(Eo);
            to = Mo*(T/2/pi);
            tf = to + n_periods*T;
            times = to:10:tf;
            options = odeset('reltol', 1.e-8, ...
                'abstol', 1.e-8, ...
                'initialstep', T/10000, ...
                'events', @(t,y) terminate(t,y));
            [~, State] = ode89(@(t,y) dfdf(t,y,t000),times, X0, options);
            xv = State(:,1);
            yv = State(:,2);
            zv = State(:,3);
            xv1 = State(:,4);
            yv2 = State(:,5);
            zv3 = State(:,6);
            rhat=R0/norm(R0);
            vtmp =  V0/ mu;
            hv = cross(R0, V0);
            hmag = sqrt(dot(hv, hv));
            ecc = cross(vtmp, hv);
            ecc = ecc - rhat;
            Ee=norm(ecc);
            Pp=aa*(1-Ee*Ee);
            rmag=Pp/(1+Ee*cos(TA*pi/180));
            rp=aa*(1-Ee);
            vp=hmag/rp;
            ra=aa*(1+Ee);
            va=hmag/ra;
            space_color = 'k';
            npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
            alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
            %GMST0 = []; % Don't set up rotatable globe (ECEF)
            GMST0 = 4.89496121282306; % Set up a rotatable globe at J2000.0

            % Earth texture image
            % Anything imread() will handle, but needs to be a 2:1 unprojected globe
            % image.

            image_file = 'EARTH.jpg';

            % Mean spherical earth

            erad    = 6371.0087714; % equatorial radius (meters)
            prad    = 6371.0087714; % polar radius (meters)
            erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

            %% Create figure


            hold on;
            axis equal;
            view(30,0);
            %% Create wireframe globe

            % Create a 3D meshgrid of the sphere points using the ellipsoid function
            [xe, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);

            app.earh = surf(xe, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

            if ~isempty(GMST0)
                app.hgx = hgtransform;
                set(app.earh,'Parent',app.hgx);
            end

            %% Texturemap the globe

            % Load Earth image for texture map
            cdata = imread(image_file);

            % Set image as color data (cdata) property, and set face color to indicate
            % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.

            set(app.earh, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
            rotate3d on
            rot3d
            view(0,0);
            axis off;
            hold on
            eq_alpha = 0.28; % transparancy
            node=cross([0 0 1],hv);
            nnode=mArrow3(-node/2,node/2,'color','#808000','stemWidth',100);
            plot3(xx,xv,yv,zv,'.r');
            hold on

            [app.sat1,app.sat2,app.sat3,app.sat4,app.sat5,app.sat6] =AliCube(a/10,xv(end),yv(end),zv(end));
            axis_data = get(gca);
            xmin = axis_data.XLim(1);
            xmax = axis_data.XLim(2);
            ymin = axis_data.YLim(1);
            ymax = axis_data.YLim(2);
            zmin = axis_data.ZLim(1);
            zmax = axis_data.ZLim(2);
            xe = [xmax xmin;xmax xmin]; ye = [ymax ymax;ymin ymin]; ze = [0 0; 0 0];
            mesh(xe,ye,ze,'FaceAlpha',eq_alpha,'FaceColor',[0.753,0.753,0.753]);
            %
            h2 = mArrow3([0 0 0],[0 ymax 0],'color','blue','stemWidth',130);
            h1 = mArrow3([0 0 0],[xmax 0 0],'color','red','stemWidth',130);
            h3 = mArrow3([0 0 0],[0 0 zmax],'color','green','stemWidth',130);

            rhat=[xv(end) yv(end) zv(end)]/norm([xv(end) yv(end) zv(end)]);
            vtmp =  [xv1(end) yv2(end) zv3(end)]/ mu;
            hv = cross([xv(end) yv(end) zv(end)], [xv1(end) yv2(end) zv3(end)]);
            %             hmag = sqrt(dot(hv, hv));
            ecc = cross(vtmp, hv);
            ecc = ecc - rhat;
            %             er=norm(ecc)
            % ec= mArrow3([0 0 0],xmax*ecc,'color','white','stemWidth',100);
            %            h4 = mArrow3([0 0 0],hv/2.5,'color','#FFA500','stemWidth',100);
            [~, ~, u,~] = solar_position((10/86400)*(length(xv))+t000);
            app.sun=mArrow3([0 0 0],u*xmax,'color','yellow','stemWidth',100);
            raa=a*(1-e);
            app.XV=xv;
            app.YV=yv;
            app.ZV=zv;
            app.v1=xv1;
            app.v2=yv2;
            app.v3=zv3;
            view(90,0);
            %             drawnow;
            %                  r=norm([xv(end) yv(end) zv(end)]);
            %
            %  speedText=uicontrol(gcf,'Style','text','String',['RMAG:',num2str(r),' Km ',' VMAG=:',num2str(sqrt(mu*((2/r-(1/a))))),' Km/s'],'Units','normalized','Position',[0.007 0.01 0.4 0.1],'BackgroundColor','none','ForegroundColor','red','FontSize',10);
            %               m=1;
            %              while m<length(xv)
            %
            %                     degPerIteration=360/(24*60*60)*200*1.25;
            %                     rotate(app.earh,[0 0 90],degPerIteration);
            %                     m=m+200;
            %                     drawnow
            %              end
            %                     drawnow limitrate;
            %    hold(xx,'off')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Orbit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Rsults%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rfile='Files\Rsult.orbit';
            fr=fopen(Rfile,'w');
            fprintf(fr,'******************************************************Orbital Rsults*************************************************************** \n\n');

            fprintf(fr,'2-State Vectot:\n');
            fprintf(fr,'__________________________________________\n');
            fprintf(fr,' 1-X= %f km \n',xv(end));
            fprintf(fr,' 2-Y= %f km  \n',yv(end));
            fprintf(fr,' 3-Z= %f km \n',zv(end));
            fprintf(fr,' 4-VX= %f km/S \n',xv1(end));
            fprintf(fr,' 5-VY= %f km/S \n',yv2(end));
            fprintf(fr,' 6-VZ= %f km/S \n\n',zv3(end));
            endfile='Files\end.txt';
            eend=fopen(endfile,'w');
            fprintf(eend,'X=%f\n',xv(end));
            fprintf(eend,'Y=%f\n',yv(end));
            fprintf(eend,'Z=%f\n',zv(end));
            fprintf(eend,'V1=%f\n',xv1(end));
            fprintf(eend,'V2=%f\n',yv2(end));
            fprintf(eend,'V3=%f\n',zv3(end));
            tt0=(10/86400)*(length(xv))+t0;
            fprintf(eend,'t0=%f\n',tt0);
            coe=ctok([xv(end) yv(end) zv(end)]',[xv1(end) yv2(end) zv3(end)]');
            e2 = coe(2);
            RA = coe(3)*180/pi;
            incl2 = coe(4)*180/pi;
            w = coe(5)*180/pi;
            TA22 = coe(6)*180/pi;
            a2 = coe(7);
            EA=2*atan(sqrt((1-e2)/(1+e2))*tan((TA22*(pi/180))/2));
            MA=EA-e2*sin(EA);
            fprintf(fr,'1-Keplerian Elements:\n');
            fprintf(fr,'__________________________________________\n');
            fprintf(fr,' 1-Semi Majoraxis = %f km \n',a2);
            fprintf(fr,' 2-eccentricity = %f  \n',e2);
            fprintf(fr,' 3-Inclination = %f° \n',incl2);
            fprintf(fr,' 4-Rightessention = %f° \n',RA);
            fprintf(fr,' 5-Argument Of Perigee = %f° \n',w);
            fprintf(fr,' 6-Mean Anomaly = %f° \n',MA*180/pi);
            fprintf(fr,' 7-Eccentric Anomaly = %f° \n',EA*180/pi);
            fprintf(fr,' 8-True Anomaly = %f° \n\n',TA22);

            %    fprintf(fr,'__________________________________________ \n\n');

            %    fprintf(fr,'__________________________________________\n\n');



            %%%%%%%%%%%%%%%%%%%%%%%%%%%Maneuver Missions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             dV=str2double( regexp(f,'(?<=V=)([+-]?([0-9]*[.])?[0-9]+)','match'))
            %             dV1=str2double( regexp(f,'(?<=N=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            %             dV2=str2double( regexp(f,'(?<=P=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            %             S=str2double( regexp(f,'(?<=S=)([+-]?([0-9]*[.])?[0-9]+)','match'))
            %
            %
            %
            %                       R00=[xv(end) yv(end) zv(end)]';
            %                       V00=[xv1(end) yv2(end) zv3(end) ]';
            %                        dV=3;
            %                        dV1=0;
            %                        dV2=0;
            %                         R0=[R00(1) R00(2) R00(3)]';
            %                         V0=[V00(1) V00(2) V00(3)]';
            %
            %                         if V00(2)<dV
            %                             V00(2)=V00(2)-dV;
            %                         else
            %                             V00(2)=V00(2)+dV;
            %                         end
            %                         V00(1)=V00(1)+dV1;
            %                         V00(3)=V00(3)+dV2;
            %                         X0=[R00;V00]
            %                         EoM=@dfdf;
            %                         E = norm(V00)^2/2 - mu/norm(R00);
            %                         aa = -mu/(2*E);
            %                         n = sqrt(abs(mu/abs(aa^3)));
            %                         T = abs(2*pi/n);
            %                         % set(handles.text28, 'String', T);
            %                         options = odeset('MaxStep', 5);
            %                         [time, State] = ode45(@(t,y) dfdf(t,y,t0(end)), [0 T*1], X0, options);
            %
            %                         xv = State(:,1);
            %                         yv = State(:,2);
            %                         zv = State(:,3);
            %                         v1 = State(:,1);
            %                         v2 = State(:,2);
            %                         v3 = State(:,3);
            %                         plot3(xv,yv,zv,'.');
            %                         R0=[xv(end) yv(end) R00(end)]';
            %                         V0=[v(end) v2(end) v3(end)]';
            %                         dV=0;
            %                        dV1=0;
            %                        dV2=0.;
            %                         R0=[R00(1) R00(2) R00(3)]';
            %                         V0=[V00(1) V00(2) V00(3)]';
            %
            %                         if V00(2)<dV
            %                             V00(2)=V00(2)-dV;
            %                         else
            %                             V00(2)=V00(2)+dV;
            %                         end
            %                         V00(1)=V00(1)+dV1;
            %                         V00(3)=V00(3)+dV2;
            %                         X0=[R00;V00]
            %                         EoM=@dfdf;
            %                         E = norm(V00)^2/2 - mu/norm(R00);
            %                         aa = -mu/(2*E);
            %                         n = sqrt(abs(mu/abs(aa^3)));
            %                         T = abs(2*pi/n);
            %                         % set(handles.text28, 'String', T);
            %                         options = odeset('MaxStep', 5);
            %                         [time, State] = ode45(EoM, [0 T*1], X0, options);
            %
            %                         xv = State(:,1);
            %                         yv = State(:,2);
            %                         zv = State(:,3);
            %                         v1 = State(:,1);
            %                         v2 = State(:,2);
            %                         v3 = State(:,3);
            %                         plot3(xv,yv,zv,'.');
            %                         R0=[xv(end) yv(end) R00(end)]';
            %                         V0=[v(end) v2(2) v3(end)]';
            %

            %%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%ground traces%%%%%%%%%%%%%%

            app.f2=figure;
            set(app.f2, 'position', [k1(1)+150 k1(2)+50 k1(3)+150 k1(4)-50],'Toolbar','none','menubar','none','name','GroundTrack','NumberTitle','off','Color','black');
            %
            Re = 6378;
            %
            ra = zeros(1,length(xv));
            dec =zeros(1,length(xv));
            oo=[];
            %             theta = 0;
            app.tt=[];

            t0=t0:10/86400:(10/86400)*(length(xv)-1)+t0;
            %             t0=linspace(t0,10/86400,(10/86400)*length(xv)+t0);
            theta=gmst(app,t0-2400000.5);
            %            the= size(theta)
            %            xxvv= size(xv)
            %             sts=size(State(:,1:3))
            % syms u
            %            ref = zeros(length(times),3);

            r_rel=zeros(length(xv),3);
            wgs84 = wgs84Ellipsoid('kilometers');
            R5=State(:,1:3)';
            %             rell=exp(theta*1i);
            rell=exp(theta*1i).*([1i 1 0]*R5);
            RRR=[imag(rell);real(rell);R5(3,:)];
            [ra, dec] = ra_and_dec_from_rN(RRR');
            if ccc==1
                for i = 1:length(xv)
                    %                 t = times(i);
                    %                 M = 2*pi/T*t;
                    %                 % E = double(vpasolve(u-e*sin(u)==M));
                    %                 E=kepler_E(app,e, M);
                    %                 TA = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
                    %                 r =hh^2/mu/(1 + e*cos(TA))*[cos(TA) sin(TA) 0]';
                    %                 W = Wo + Wdot*t;
                    %                 wp = wpo + wpdot*t;
                    %                 R11 = [ cos(W) sin(W) 0
                    %                     -sin(W) cos(W) 0
                    %                     0 0 1];
                    %                 R2 = [1 0 0
                    %                     0 cos(incl) sin(incl)
                    %                     0 -sin(incl) cos(incl)];
                    %                 R3 = [ cos(wp) sin(wp) 0
                    %                     -sin(wp) cos(wp) 0
                    %                     0 0 1];
                    %                 QxX = (R3*R2*R11)';
                    %                 R = QxX*r;
                    %                   R=
                    %                 t0 = t0+10/86400;


                    %                  kk= datetime(t0(i),'convertfrom','juliandate');





                    %                  if(i==length(times))
                    %                       app.END=theta;
                    %                  end

                    %                theta(i)=gmst(app,t0(o)-2400000.5);
                    %                 Q = R_z(theta(i));
                    %                 r_rel(i,1:3)=Q*State(i,1:3)';


                    if ccc==1


                        try

                            Re=geocradius(lat)/1000;
                            thi=long+theta(i)*180/pi;
                            RR=Re*[(cosd(thi)*cosd(lat));(sind(thi)*cosd(lat));sind(lat)];
                            value=acos(dot(RR,State(i,1:3))/(norm(RR)*norm(State(i,1:3))))*180/pi;
                            delt=3/pi*acosd(Re/norm(State(i,1:3)));

                            %                    [~, ~, ~,r_S] = solar_position(t0(i));
                            %                    nu = los(State(i,1:3), r_S);
                            %               [~,elev,~] = ecef2aer( r_rel(i,1), r_rel(i,2), r_rel(i,3),lat,long,0,wgs84);
                            %               alt=norm(State(i,1:3))-6378;
                            %               k=6378/(6378+alt);
                            %               EE=acosd(k*cosd(elev))-elev


                            %               value
                            if(value<=delt)

                                %                         if(elev<=90)
                                kk=t0(i)-t000;
                                timef=t00+seconds(kk*24*60*60); an=timef; oo=[an;oo];
                                %                         end
                            end
                        catch M=0;

                        end

                    end

                    %                  if(norm(r_rel(i,:))<=6371)
                    %
                    %                     ii=i;
                    %
                    %
                    %                     break
                    %                 end

                    %                 [alpha, delta] = ra_and_dec_from_r(app,r_rel);
                    %                 ra = [ra; alpha];
                    %                 dec = [dec; delta];
                    %  if(ra(i)>180)
                    %      ra(i)=ra(i)-360;
                    %      end

                    %                 [ra(i), dec(i)] = ra_and_dec_from_r(app,r_rel(i,:));


                end
            end

            %               ra(end)
            app.start=theta(1);
            app.tt=datetime(t0,'convertfrom','juliandate');
            if ccc==1




                %                 event_no = 1;
                %                 n_event = 1;
                %                 n = 0;
                %                 timed = datevec(oo(1));
                %                 for i = 1:length(oo)
                %                     if abs(etime(datevec(oo(i)) , timed)) >=20
                %                         event_no = event_no + 1;
                %                         n_event= n_event + 1;
                %                         n = 0;
                %                     end
                %                     n = n + 1;
                %
                %                     T1{event_no}(n) = oo(i);
                %                     timed = datevec(oo(i));
                %                 end
                %
                %                 ali=0;
                try
                    T1=sep(oo);
                catch al=0;
                end

            end
            if ccc==1


                try


                    oo(1);
                    file11='Files\pass.orbit';
                    fileID11 = fopen(file11,'w');
                    fprintf(fileID11,'********************pass time******************************\n');
                    fprintf(fileID11,'Start(UTC)\t\t\tEnd(UTC)\t\t\tDuration (s)\n');
                    child = uitreenode(app.OutputTree,'Text','passTime',"Icon",'report.png');
                    for i=length(T1):-1:1
                        start1=datevec(T1{i}(end));endd=datevec(T1{i}(1));

                        duration=abs(etime(start1 , endd));


                        fprintf(fileID11,'%s\t%s\t%f\n',T1{i}(end),T1{i}(1),duration);

                    end

                catch ali=[];

                end

            end


            %find_ra_and_dec
            % comet(ra,dec)
            %
            % Breaks the ground track up into separate curves which start
            % and terminate at right ascensions in the range [0,360 deg].
            % –––––––––––––––––––––––––––
            %             tol = 100;
            %             curve_no = 1;
            %             n_curves = 1;
            %             k = 0;
            %             ra_prev = ra(1);
            %             for i = 1:length(ra)
            %                 if abs(ra(i) - ra_prev) > tol
            %                     curve_no = curve_no + 1;
            %                     n_curves = n_curves + 1;
            %                     k = 0;
            %                 end
            %                 k = k + 1;
            %                 RA1{curve_no}(k) = ra(i);
            %                 Dec{curve_no}(k) = dec(i);
            %                 ra_prev = ra(i);
            %             end
            %             [n_curves,RA1, Dec]=seperate(ra,dec);
            %form_separate_curves

            %
            o=imread('EARTH.jpg');
            xlim([-180 180]);
            ylim([-90 90]);
            hold on
            image([-180 180],[90 -90],o)

            xlabel('longitude (degrees)');
            ylabel('Latitude (degrees)');
            axis equal
            grid on
            for i = 1:11
                plot([-180+i*30,-180+i*30],[-90,90],'Color','#C0C0C0',...
                    'LineWidth',0.5,'LineStyle','-');
            end
            for i = 1:5
                plot([-180,180],[-90+i*30,-90+i*30],'color','#C0C0C0',...
                    'LineWidth',0.5,'LineStyle','-');
            end
            %             for i = 1:n_curves
            %
            %
            %                 plot(RA1{i}, Dec{i},'-r','LineWidth',2)
            %                 drawnow
            %             end
            plot(ra,dec,'.r')
            axis ([-180 180 -90 90]);
            xticks(-180:30:180);
            yticks(-90:30:90);
            % text( ra(1), dec(1), 'o Start')
            % text(ra(end), dec(end), 'o Finish')
            % set(ha,'position',[ra(end) dec(end) 0.1 0.1])
            % S=imread('sat2.png');
            [S, ~, ImageAlpha] = imread('sat2.png');
            S2 = imread('groundA.png');
            app.I = imagesc('XData', [ra(end)-10 ra(end)+10], 'YData',...
                [dec(end)-10 dec(end)+10],'CData', S,'AlphaData',ImageAlpha);
            if ccc==1
                ground= imagesc('XData', [long-10 long], 'YData',...
                    [lat-10 lat],'CData', S2);


            end

            app.text1=text(-160, -65,['epoch : ',datestr(app.tt(end))],'Color','red','FontSize',12,'fontweight','bold');
            app.Text=uicontrol(app.h,'Style','text','String',['epoch : ',datestr(app.tt(end))],'Units','normalized','Position',[0.007 0.005 0.5 0.06],'BackgroundColor','none','ForegroundColor','red','FontSize',12,'fontweight','bold');
            try
                tex=text(0, -65,['END LIFE : ',datestr(datetime(t0(ii),'convertfrom','juliandate'))],'Color','red','FontSize',10,'fontweight','bold');
            catch ali=0;
            end

            %  k=text(ra(end), dec(end),'')
            app.END=theta(end);
            theta(end)
            set(app.hgx,'Matrix', makehgtform('zrotate',app.END));
            set(gcf,'Toolbar','none','menubar','none','name','GroundTrack','NumberTitle','off','Color','#C0C0C0');
            app.RAS=ra;
            app.DECL=dec;
            % %             t = tiledlayout(app.OrbitPanel,2,1);
            %             fig=uifigure(app.MainGUI)
            % WinOnTop( app.h);
            % WinOnTop( app.f2);
            setFigureAlwaysOnTop(app.h)
            setFigureAlwaysOnTop(app.f2)
            set(app.RunButoon,'Enable','on');
            app.check=0;
            %             label=uilabel(app.f2);

            fprintf(fr,'2-Spherical State :\n');

            [rmag,ra,dec,vmag,rav,decv] = C2SP(app,[xv(end) yv(end) zv(end)],[xv1(end) yv2(end) zv3(end)]);
            %             PP=a*(1-e*e);
            fprintf(fr,'__________________________________________\n');
            fprintf(fr,'RMAG = %f km \n',rmag);
            fprintf(fr,'Right Ascension = %f deg \n',ra);
            fprintf(fr,'Declination = %f km^2/S^2 \n',dec);
            fprintf(fr,'VMAG = %f KM/S \n',vmag);
            fprintf(fr,'Vel Right Ascension = %f deg \n',rav);
            fprintf(fr,'Vel Declination  = %f deg \n\n',decv);

            fprintf(fr,'2-Other Orbit Data:\n');
            E = norm([xv1(end) yv2(end) zv3(end)])^2/2 - mu/norm([xv(end) yv(end) zv(end)]);
            aa = -mu/(2*E);
            n = sqrt(abs(mu/abs(aa^3)));
            T = abs(2*pi/n);
            %             PP=a*(1-e*e);
            rp=a2*(1-e2);
            hv = cross([xv(end) yv(end) zv(end)], [xv1(end) yv2(end) zv3(end)]);
            hmag = sqrt(dot(hv, hv));
            vp=hmag/rp;
            ra=a2*(1+e2);
            va=hmag/ra;
            p=sqrt(rp*ra);
            mo=(2*pi)/T;
            fprintf(fr,'__________________________________________\n');
            fprintf(fr,'Semi minor Axis = %f km \n',p);
            fprintf(fr,'specific angular momentum = %f KM^2/S \n',hmag);
            fprintf(fr,'Orbit Energy = %f km^2/S^2 \n',E);
            fprintf(fr,'Mean motion = %f RAD/S \n',mo);
            fprintf(fr,'Periapsis = %f km \n',rp);
            fprintf(fr,'VelPeriapsis  = %f km/S \n',vp);
            fprintf(fr,'Apoapsis = %f km \n',ra);
            fprintf(fr,'VelApoapsis = %f km/S \n',va);
            fprintf(fr,'Orbit Period = %f S \n',T);



            % app.c=2

            %  while(1)
            %
            %  if (app.c)==1
            %
            %
            %  for i = 3 : length(ra)
            % %         title(date(i));
            %         delete(I);
            %          delete(k)
            %         I = imagesc('XData', [ra(i)-5 ra(i)+5],...
            %             'YData', [dec(i) - 5 dec(i) + 5], 'CData', S,'AlphaData',ImageAlpha);
            % %        k=text(ra(11), dec(108), num2str(ra(i)))
            % %        k=text(ra(-100), dec(-190),['Latitude:' ,num2str(ra(i)),' Longitude:',num2str(dec(i))]);
            %               k=title(['Latitude:' ,num2str(ra(i)),' Longitude:',num2str(dec(i))]);
            %
            %
            % %         title(gca,['Latitude:' ,ra(i),' Longitude:',dec(i)]);
            %
            %         % other stuff (plots)
            %         drawnow
            %         if app.c==2
            %      break
            %  end
            %  end
            %  end
            %  if app.c==2
            %      break
            %  end
            % end


            % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
        end

        % Close request function: MainGUI
        function MainGUICloseRequest(app, event)
            % %              if app.x==0
            % %                 value=0
            % %             end
            %
            %
            %             switch value
            %                 case 'Space Craft'
            %                     value=1
            %                     appl
            %                     t=app.x;
            %                     t.TextAreaLabel.Text
            %                     app.ListBox.Value='';
            % %                     value
            % %
            % app.appl.Visible = 'on';
            % fig = app.MainGUI
            % msgbox('Close document?','Confirm Close',...
            %                         'Icon','warning');
            % selection = uiconfirm( 'Close document?','Confirm Close',...
            %                         'Icon','warning');
            % switch selection
            %     case'ok'
            %         close all
            %         delete(app)

            % msgbox('15');
            close all
            delete(app)
            close all


        end

        % Size changed function: MainGUI
        function MainGUISizeChanged(app, event)
            %             app.MainGUI.Position(3:4) = max(app.MainGUI.Position(3:4), [500 500]),min(app.MainGUI.Position(3:4), [100 100]);
            %             movegui(app.MainGUI)
        end

        % Callback function
        function ListBoxValueChanged(app, event)
            value = app.ListBox.Value;
            %              if app.x==0
            %                 value=0
            %             end


            switch value
                case 'Space Craft'

                    app22;
                    %                     t=app.x;
                    %                     t.TextAreaLabel.Text
                    app.ListBox.Value='';
                    %                     value
                    %                     app.appl.Visible = 'on';


            end
        end

        % Callback function: Restart
        function RestartClicked(app, event)
            closereq
            close all
            SMAD_app;
        end

        % Callback function: PushHOLDTool4
        function PushHOLDTool4Clicked(app, event)
            try

                clf(app.h)
                close(app.f2)
            catch M5=0;

            end
            xx=axes(app.h);
            % hold(app.h,'on');
            file='Files\end.txt';
            %             file1='C:\Users\Ali\Desktop\pass.orbit';
            f=fileread(file);
            X=str2double( regexp(f,'(?<=X=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            Y=str2double( regexp(f,'(?<=Y=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            Z=str2double( regexp(f,'(?<=Z=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            VX=str2double( regexp(f,'(?<=V1=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            VY=str2double( regexp(f,'(?<=V2=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            VZ=str2double( regexp(f,'(?<=V3=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            t000=str2double( regexp(f,'(?<=t0=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            mu = 3.986012E5;
            R0=[X Y Z]';
            V0=[VX VY VZ]';
            X0=[R0;V0];

            coe=ctok(R0,V0);
            e = coe(2);

            TAo = coe(6);
            E = norm(V0)^2/2 - mu/norm(R0);
            aa = -mu/(2*E);
            n = sqrt(abs(mu/abs(aa^3)));
            T = abs(2*pi/n);
            n_periods=1;
            %             TAo = TA*pi/180;
            Eo = 2*atan(tan(TAo/2)*sqrt((1-e)/(1+e)));
            Mo = Eo - e*sin(Eo);
            to = Mo*(T/2/pi);
            tf = to + n_periods*T;
            times = [to:10:tf];
            % set(handles.text28, 'String', T);
            options = odeset('reltol', 1.e-8, ...
                'abstol', 1.e-8, ...
                'initialstep', T/10000, ...
                'events', @(t,y) terminate(t,y));

            %          for i=1:30
            %                values = [1:10:5000];
            [time, State] = ode45(@(t,y) dfdf(t,y,t000),times, X0, options);

            %          end
            %   [time1, State1] = ode45(EoM, [0 T], R0, options);

            xv = State(:,1)
            yv = State(:,2);
            zv = State(:,3);
            xv1 = State(:,4);
            yv2 = State(:,5);
            zv3 = State(:,6);
            %             vv=[xv yv zv];
            %             vvv=[];
            %             for i =1:numel(time)
            %
            %                 vvv(i)=norm(vv(i:i,1:3));
            %             end

            % %  vvv=norm(vv)

            R = 6371;
            rhat=R0/norm(R0);
            vtmp =  V0/ mu;
            hv = cross(R0, V0);
            hmag = sqrt(dot(hv, hv));
            ecc = cross(vtmp, hv);
            ecc = ecc - rhat;
            Ee=norm(ecc);
            Pp=aa*(1-Ee*Ee);
            rmag=Pp/(1+Ee*cos(TAo));
            v=hmag/rmag;
            rp=aa*(1-Ee);
            vp=hmag/rp;
            ra=aa*(1+Ee);
            va=hmag/ra;
            space_color = 'k';
            npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
            alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
            %GMST0 = []; % Don't set up rotatable globe (ECEF)
            GMST0 = 4.89496121282306; % Set up a rotatable globe at J2000.0

            % Earth texture image
            % Anything imread() will handle, but needs to be a 2:1 unprojected globe
            % image.

            image_file = 'EARTH.jpg';

            % Mean spherical earth

            erad    = 6371.0087714; % equatorial radius (meters)
            prad    = 6371.0087714; % polar radius (meters)
            erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

            %% Create figure


            hold on;

            % Turn off the normal axes

            % set(xx, 'NextPlot','add', 'Visible','off');

            axis 'equal';
            % axis auto;

            % Set initial view

            view(30,0);

            %             axis vis3d;

            %% Create wireframe globe

            % Create a 3D meshgrid of the sphere points using the ellipsoid function

            [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);

            app.earh = surf(xx,x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

            if ~isempty(GMST0)
                app.hgx = hgtransform;
                %                 set(app.hgx,'Matrix', makehgtform('zrotate',GMST0));
                set(app.earh,'Parent',app.hgx);
            end

            %% Texturemap the globe

            % Load Earth image for texture map

            cdata = imread(image_file);

            % Set image as color data (cdata) property, and set face color to indicate
            % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.

            set(app.earh, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
            rotate3d on
            rot3d
            view(0,0);
            axis off;
            hold on
            %             xe = [1.1371e+04 -1.1371e+04;1.1371e+04 -1.1371e+04]; ye = [1.1371e+04 1.1371e+04;-1.1371e+04 -1.1371e+04]; ze = [0 0; 0 0];
            eq_alpha = 0.28; % transparancy
            %             mesh(xe,ye,ze,'FaceAlpha',eq_alpha,'FaceColor',[0.753,0.753,0.753]);

            %             il=acos(hv(3)/norm(hv))*180/pi;
            node=cross([0 0 1],hv);
            %             ww=acos(dot(node,ecc)/(norm(node)*norm(ecc)))*180/pi;
            %             ras=acos(node(1)/norm(node))*180/pi;


            %  h1 = plot([0 6371+5000],[0 0],'-r');
            %    h2 = plot([0 0],[0 6371+5000],'-g');
            %    h3 = plot3([0 0],[0 0],[0 6371+5000],'color',[0 0 .8]);
            %    set([h1 h2 h3],'linewidth',2);
            %
            %    plot3(6371+20000,0,0,'red','Linewidth',1);
            %   plot3(0,6371+20000,0,'greenY','Linewidth',1);
            % plot3(0,0,6371+20000,'blueZ','Linewidth',1);
            nnode=mArrow3(-node/2,node/2,'color','#808000','stemWidth',100);
            plot3(xx,xv,yv,zv,'.r');
            hold on
            a=aa;
            [app.sat1,app.sat2,app.sat3,app.sat4,app.sat5,app.sat6] =AliCube(a/10,xv(end),yv(end),zv(end));
            axis_data = get(gca);
            xmin = axis_data.XLim(1);
            xmax = axis_data.XLim(2);
            ymin = axis_data.YLim(1);
            ymax = axis_data.YLim(2);
            zmin = axis_data.ZLim(1);
            zmax = axis_data.ZLim(2);
            %                  xlim([xmin xmax]);
            %                  ylim([ymin ymax]);
            %                  zlim([zmin zmax]);
            xe = [xmax xmin;xmax xmin]; ye = [ymax ymax;ymin ymin]; ze = [0 0; 0 0];
            mesh(xe,ye,ze,'FaceAlpha',eq_alpha,'FaceColor',[0.753,0.753,0.753]);
            %
            h2 = mArrow3([0 0 0],[0 ymax 0],'color','blue','stemWidth',130);
            h1 = mArrow3([0 0 0],[xmax 0 0],'color','red','stemWidth',130);
            h3 = mArrow3([0 0 0],[0 0 zmax],'color','green','stemWidth',130);

            rhat=[xv(end) yv(end) zv(end)]/norm([xv(end) yv(end) zv(end)]);
            vtmp =  [xv1(end) yv2(end) zv3(end)]/ mu;
            hv = cross([xv(end) yv(end) zv(end)], [xv1(end) yv2(end) zv3(end)]);
            %             hmag = sqrt(dot(hv, hv));
            ecc = cross(vtmp, hv);
            ecc = ecc - rhat;
            %             er=norm(ecc)
            ec= mArrow3([0 0 0],xmax*ecc,'color','white','stemWidth',100);
            %            h4 = mArrow3([0 0 0],hv/2.5,'color','#FFA500','stemWidth',100);
            [~, ~, u,~] = solar_position((10/86400)*(length(xv))+t000);
            app.sun=mArrow3([0 0 0],u*xmax,'color','yellow','stemWidth',100);


            app.XV=xv;
            app.YV=yv;
            app.ZV=zv;
            app.v1=xv1;
            app.v2=yv2;
            app.v3=zv3;
            view(90,0);
            %             WinOnTop(figure(3))
            %             WinOnTop(figure(4))
            %       delete(app)
            %       close all
            %       SMAD_app

            %       app.UIFigure

            app.f2=figure;
            k1=app.GroundPanel.Position;
            set(app.f2, 'position', [k1(1)+150 k1(2)+50 k1(3)+150 k1(4)-50],'Toolbar','none','menubar','none','name','GroundTrack','NumberTitle','off','Color','black');
            ra = zeros(1,length(xv));
            dec =zeros(1,length(xv));
            oo=[];
            %             theta = 0;
            app.tt=[];

            t0=[t000:10/86400:(10/86400)*(length(xv))+t000];
            theta=gmst(app,t0-2400000.5);
            % syms u
            %            ref = zeros(length(times),3);

            r_rel=zeros(length(xv),3);
            %           wgs84 = wgs84Ellipsoid('kilometers');
            for i = 1:length(xv)
                %                 t = times(i);
                %                 M = 2*pi/T*t;
                %                 % E = double(vpasolve(u-e*sin(u)==M));
                %                 E=kepler_E(app,e, M);
                %                 TA = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
                %                 r =hh^2/mu/(1 + e*cos(TA))*[cos(TA) sin(TA) 0]';
                %                 W = Wo + Wdot*t;
                %                 wp = wpo + wpdot*t;
                %                 R11 = [ cos(W) sin(W) 0
                %                     -sin(W) cos(W) 0
                %                     0 0 1];
                %                 R2 = [1 0 0
                %                     0 cos(incl) sin(incl)
                %                     0 -sin(incl) cos(incl)];
                %                 R3 = [ cos(wp) sin(wp) 0
                %                     -sin(wp) cos(wp) 0
                %                     0 0 1];
                %                 QxX = (R3*R2*R11)';
                %                 R = QxX*r;
                %                   R=
                %                 t0 = t0+10/86400;


                %                  kk= datetime(t0(i),'convertfrom','juliandate');





                %                  if(i==length(times))
                %                       app.END=theta;
                %                  end

                %                theta(i)=gmst(app,t0(o)-2400000.5);
                Q = R_z(theta(i));
                r_rel(i,1:3)=Q*State(i,1:3)';







                %                  if(norm(r_rel(i,:))<=6371)
                %
                %                     ii=i;
                %
                %
                %                     break
                %                 end

                %                 [alpha, delta] = ra_and_dec_from_r(app,r_rel);
                %                 ra = [ra; alpha];
                %                 dec = [dec; delta];
                %  if(ra(i)>180)
                %      ra(i)=ra(i)-360;
                %      end

                [ra(i), dec(i)] = ra_and_dec_from_r(app,r_rel(i,:));


            end

            %               ra(end)
            app.start=theta(1);
            app.tt=datetime(t0,'convertfrom','juliandate');




            %find_ra_and_dec
            % comet(ra,dec)
            %
            % Breaks the ground track up into separate curves which start
            % and terminate at right ascensions in the range [0,360 deg].
            % –––––––––––––––––––––––––––
            %             tol = 100;
            %             curve_no = 1;
            %             n_curves = 1;
            %             k = 0;
            %             ra_prev = ra(1);
            %             for i = 1:length(ra)
            %                 if abs(ra(i) - ra_prev) > tol
            %                     curve_no = curve_no + 1;
            %                     n_curves = n_curves + 1;
            %                     k = 0;
            %                 end
            %                 k = k + 1;
            %                 RA1{curve_no}(k) = ra(i);
            %                 Dec{curve_no}(k) = dec(i);
            %                 ra_prev = ra(i);
            %             end
            [n_curves,RA1, Dec]=seperate(ra,dec);
            %form_separate_curves

            %
            o=imread('EARTH.jpg');
            xlim([-180 180]);
            ylim([-90 90]);
            hold on
            image([-180 180],[90 -90],o)

            xlabel('longitude (degrees)');
            ylabel('Latitude (degrees)');
            axis equal
            grid on
            for i = 1:11
                plot([-180+i*30,-180+i*30],[-90,90],'Color','#C0C0C0',...
                    'LineWidth',0.5,'LineStyle','-');
            end
            for i = 1:5
                plot([-180,180],[-90+i*30,-90+i*30],'color','#C0C0C0',...
                    'LineWidth',0.5,'LineStyle','-');
            end
            for i = 1:n_curves


                plot(RA1{i}, Dec{i},'-r','LineWidth',2)
                drawnow
            end
            axis ([-180 180 -90 90]);
            xticks([-180:30:180]);
            yticks([-90:30:90]);
            % text( ra(1), dec(1), 'o Start')
            % text(ra(end), dec(end), 'o Finish')
            % set(ha,'position',[ra(end) dec(end) 0.1 0.1])
            % S=imread('sat2.png');
            [S, ~, ImageAlpha] = imread('sat2.png');
            S2 = imread('groundA.png');
            app.I = imagesc('XData', [ra(end)-10 ra(end)+10], 'YData',...
                [dec(end)-10 dec(end)+10],'CData', S,'AlphaData',ImageAlpha);

            app.text1=text(-160, -65,['epoch : ',datestr(app.tt(end))],'Color','red','FontSize',12,'fontweight','bold');
            app.Text=uicontrol(app.h,'Style','text','String',['epoch : ',datestr(app.tt(end))],'Units','normalized','Position',[0.007 0.005 0.5 0.06],'BackgroundColor','none','ForegroundColor','red','FontSize',12,'fontweight','bold');
            try
                tex=text(0, -65,['END LIFE : ',datestr(datetime(t0(ii),'convertfrom','juliandate'))],'Color','red','FontSize',10,'fontweight','bold');
            catch ali=0;
            end

            %  k=text(ra(end), dec(end),'')
            app.END=theta(end);
            set(app.hgx,'Matrix', makehgtform('zrotate',app.END));
            set(gcf,'Toolbar','none','menubar','none','name','GroundTrack','NumberTitle','off','Color','#C0C0C0');
            app.RAS=ra;
            app.DECL=dec;
            % %             t = tiledlayout(app.OrbitPanel,2,1);
            %             fig=uifigure(app.MainGUI)
            setFigureAlwaysOnTop( app.h);
            setFigureAlwaysOnTop( app.f2);

            file1='Files\ali.orbit';
            fileID1 = fopen(file1,'w');
            t=datetime(t0(1),'convertfrom','juliandate');
            tim=datevec(t);
            fprintf(fileID1,'y=%f\n',tim(1));
            fprintf(fileID1,'m=%f\n',tim(2));
            fprintf(fileID1,'d=%f\n',tim(3));
            fprintf(fileID1,'h=%f\n',tim(4));
            fprintf(fileID1,'mn=%f\n',tim(5));
            fprintf(fileID1,'sec=%f\n',tim(6));

            coe=ctok([xv(end) yv(end) zv(end)]',[xv1(end) yv2(end) zv3(end)]');
            e2 = coe(2);
            RA = coe(3)*180/pi;
            incl2 = coe(4)*180/pi;
            w = coe(5)*180/pi;
            TA22 = coe(6)*180/pi;
            a2 = coe(7);


            fprintf(fileID1,'a=%f\n',a2);
            fprintf(fileID1,'e=%f\n',e2);
            fprintf(fileID1,'i=%f\n',incl2);
            fprintf(fileID1,'TA=%f\n',TA22);
            fprintf(fileID1,'RA=%f\n',RA);
            fprintf(fileID1,'w=%f\n',w);
            fprintf(fileID1,'s=%f\n',1);

            %            C2SP(app,State)

        end

        % Callback function: AnimationButton
        function AnimationButtonClicked(app, event)
            %        f=fileread('C:\Users\Ali\Desktop\ali.orbit');
            %         f=fileread('f:orbit.orbit');
            %
            %        a=str2double( regexp(f,'(?<=a=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            %        e=str2double( regexp(f,'(?<=e=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            %        i=str2double( regexp(f,'(?<=i=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            %        w=str2double( regexp(f,'(?<=w=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            %        TA=str2double( regexp(f,'(?<=TA=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            %        RA=str2double( regexp(f,'(?<=RA=)([+-]?([0-9]*[.])?[0-9]+)','match'));
            %        app.Label.Text=num2str(a);
            %        app.Label2.Text=num2str(e);
            %        app.Label3.Text=num2str(i);
            %        app.Label4.Text=num2str(w);
            %        app.Label5.Text=num2str(TA);
            %        app.Label6.Text=num2str(RA);


            app.c=1;
            ali = app.x;
            set(app.AnimationButton,'Enable','off')
            set(app.PauseButton,'Enable','on')
            ra1=app.RAS;
            dec1=app.DECL;
            xv1=app.XV;
            yv1=app.YV;
            zv1=app.ZV;
            [S, ~, ImageAlpha] = imread('sat2.png');
            %   if(app.calc==1)
            %  I = imagesc('XData', [ra1(1)-5 ra1(1)+5], 'YData',...
            %     [dec1(1)-5 dec1(1)+5],'CData', S,'AlphaData',ImageAlpha);
            if(app.ano==0)
                delete(app.I);
                delete(app.text1);
                app.text1=text(-160, -65,['epoch : ',datestr(app.tt(1))],'Color','red','FontSize',12,'fontweight','bold');
                app.I = imagesc('XData', [ra1(1)-5 ra1(1)+5], 'YData',...
                    [dec1(1)-5 dec1(1)+5],'CData', S,'AlphaData',ImageAlpha);
            end
            if(app.ano==1)
                taa=[app.ta:10/86400:(10/86400)*(length(xv1))+app.ta];
                axis_data = get(gca);
                xmin = axis_data.XLim(1);
                xmax = axis_data.XLim(2);
                ymin = axis_data.YLim(1);
                ymax = axis_data.YLim(2);
                zmin = axis_data.ZLim(1);
                zmax = axis_data.ZLim(2);
                xlim([xmin xmax]);
                ylim([ymin ymax]);
                zlim([zmin zmax]);
                delete(app.sun);
                delete(app.sat1); delete(app.sat2); delete(app.sat3);
                delete(app.sat4); delete(app.sat5); delete(app.sat6) ;
                [app.sat1,app.sat2,app.sat3,app.sat4,app.sat5,app.sat6]=  AliCube(app.SMA/10,xv1(1),yv1(1),zv1(1));
                [~, ~, u,~] = solar_position(taa(1));
                app.sun=mArrow3([0 0 0],u*xmax,'color','yellow','stemWidth',100);
                %                 degPerIteration=(app.start*(180/pi))/(24*60*60);
                %                     rotate(app.earh,[0 90],-10000*degPerIteration,[0 0 0]);
                %           set(app.hgx,'Matrix', makehgtform('zrotate',app.start));
                %                 set(app.hgx,'Matrix',eye(4))




                % set(app.hgx,'Matrix', makehgtform('zrotate',app.start));
                % xe = [xv1(1)+5000 xv1(1)-5000;xv1(1)+5000 xv1(1)-5000]; ye = [yv1(1)+5000 yv1(1)+5000;yv1(1)-5000 yv1(1)-5000]; ze = [zv1(1)+5000 zv1(1)+5000;zv1(1)-5000 zv1(1)-5000];
                %            [M, ~, ImageAlpha] = imread('AliSat.png');
                %           app.sat= mesh(xe,ye,ze,'FaceColor', 'texturemap', 'CData', M,'FaceAlpha','texturemap','AlphaData',ImageAlpha,'EdgeColor', 'none');
            end

            %  k=text(ra1(1), dec1(1),'')


            I='0';
            if(app.ano==0)

                while app.x<length(ra1)&&(app.c==1)





                    %         title(date(i));
                    delete(app.I);
                    delete(app.text1);
                    app.text1=text(-160, -65,['epoch : ',datestr(app.tt(app.x))],'Color','red','FontSize',12,'fontweight','bold');
                    %          delete(k)
                    app.I = imagesc('XData', [ra1(app.x)-10 ra1(app.x)+10],...
                        'YData', [dec1(app.x) - 10 dec1(app.x) + 10], 'CData', S,'AlphaData',ImageAlpha);

                    %        k=text(ra(11), dec(108), num2str(ra(i)))
                    %        k=text(ra(-100), dec(-190),['Latitude:' ,num2str(ra(i)),' Longitude:',num2str(dec(i))]);
                    k=title(['Latitude:' ,num2str(ra1(app.x)),' Longitude:',num2str(dec1(app.x))]);


                    %         title(gca,['Latitude:' ,ra(i),' Longitude:',dec(i)]);

                    % other stuff (plots)


                    app.x=app.x+app.increase/10;

                    drawnow
                    if(app.x>=length(ra1))
                        delete(app.I);
                        set(app.AnimationButton,'Enable','on');
                        app.x=1;
                        app.I = imagesc('XData', [ra1(end)-10 ra1(end)+10], 'YData',...
                            [dec1(end)-10 dec1(end)+10],'CData', S,'AlphaData',ImageAlpha);
                        delete(app.text1);
                        app.text1=text(-160, -65,['epoch : ',datestr(app.tt(end))],'Color','red','FontSize',12,'fontweight','bold');
                        break;
                    end


                end
            end
            if(app.ano==1)
                GM1=398600.44189;
                a=app.SMA;

                absvalue=(sqrt(a^3/42164^3)*2*pi)*app.number;

                set(app.hgx,'Matrix', makehgtform('zrotate',app.start-absvalue*app.ABS));
                %                 if(app.num==0)
                %
                %                 end
                %                app.num=1;
                %

                while app.x<length(xv1)&&(app.c==1)



                    app.Text=uicontrol(app.h,'Style','text','String',['epoch : ',datestr(app.tt(app.x))],'Units','normalized','Position',[0.007 0.005 0.5 0.06],'BackgroundColor','none','ForegroundColor','red','FontSize',12,'fontweight','bold');

                    degPerIteration=360/(24*60*6)*app.increase;
                    rotate(app.earh,[0 90],degPerIteration,[0 0 0]);
                    delete(app.sun);
                    delete(app.sat1); delete(app.sat2); delete(app.sat3);
                    delete(app.sat4); delete(app.sat5); delete(app.sat6);
                    try
                        r=norm([xv1(app.x) yv1(app.x) zv1(app.x)]);
                        k=title(['RMAG:',num2str(r),' Km ',' Velocity:',num2str(sqrt(GM1*((2/r-(1/a))))),' Km/s'],'Color','white');
                    catch
                        continue

                    end

                    %                                 xe = [xv1(app.x)+5000 xv1(app.x)-5000;xv1(app.x)+5000 xv1(app.x)-5000]; ye = [yv1(app.x)+5000 yv1(app.x)+5000;yv1(app.x)-5000 yv1(app.x)-5000]; ze = [zv1(app.x)+5000 zv1(app.x)+5000;zv1(app.x)-5000 zv1(app.x)-5000];
                    %                     %
                    %                     %
                    %                     %           app.sat= mesh(xe,ye,ze,'FaceColor', 'texturemap', 'CData', M,'FaceAlpha','texturemap','AlphaData',ImageAlpha,'EdgeColor', 'none');
                    [app.sat1,app.sat2,app.sat3,app.sat4,app.sat5,app.sat6]=  AliCube(app.SMA/10,xv1(app.x),yv1(app.x),zv1(app.x));
                    [~, ~, u,~] = solar_position(taa(app.x));
                    app.sun=mArrow3([0 0 0],u*xmax,'color','yellow','stemWidth',100);
                    %
                    %
                    %                     %
                    %                     % pause(0.00000000000000000000000000000000000001)
                    %
                    %
                    app.x=app.x+app.increase;
                    %
                    %
                    if(app.x>length(xv1))
                        app.ABS=app.ABS+1;
                        app.num=0;
                        app.inc=degPerIteration;
                        axis_data = get(gca);
                        delete(app.sat1); delete(app.sat2); delete(app.sat3);
                        delete(app.sat4); delete(app.sat5); delete(app.sat6);
                        set(app.AnimationButton,'Enable','on')
                        app.x=1;
                        [app.sat1,app.sat2,app.sat3,app.sat4,app.sat5,app.sat6]=  AliCube(app.SMA/10,xv1(end),yv1(end),zv1(end));
                        delete(app.sun);
                        [~, ~, u,~] = solar_position(taa(end));
                        app.sun=mArrow3([0 0 0],u*xmax,'color','yellow','stemWidth',100);
                        %                    set(app.hgx,'Matrix', makehgtform('zrotate',app.END));
                        % zmax = axis_data.ZLim(2);
                        %
                        %             h3 = plot3([0 0],[0 0],[0 20000],'color','red');
                        break;
                    end
                    %
                    % %                     axis equal


                    drawnow limitrate;
                end
            end

        end

        % Callback function: SaveChanges
        function SaveChangesClicked(app, event)
            % set(app.hgx,'Matrix', makehgtform('zrotate',app.start-2));
            % app.start
            % GMST0 = 4.89496121282306;
            %                 set(app.hgx,'Matrix', makehgtform('zrotate',GMST0));
        end

        % Menu selected function: addmaneuverMenu
        function addmaneuverMenuSelected(app, event)
            %                 parent=app.MissionTreee.SelectedNodes;
            app.man = uitreenode(app.OrbitNode,'Text',append('maneuver',num2str(app.count)),'Icon','G:deltav.png','tag','man');

            %             sat6=[app.sat1,app.sat2,app.sat3,app.sat4,app.sat5,app.sat6];
            maneuvr(app.count);
            app.count=app.count+1;

            % app.OrbitNode=uitreenode()
        end

        % Selection changed function: MissionTreee
        function MissionTreeeSelectionChanged(app, event)
            selectedNodes = app.MissionTreee.SelectedNodes;
            switch selectedNodes.Text
                case 'Orbit'
                    1+1;
                otherwise
                    file="Files\result.txt";
                    s=selectedNodes.Text(9:length(selectedNodes.Text));
                    %               S=str2num(s);
                    %             file1='C:\Users\Ali\Desktop\pass.orbit';
                    f=fileread(file);
                    t=extractBetween(f,append('smaneuver',s),append('emaneuver',s));
                    uifig = uifigure('WindowStyle','alwaysontop','Name',append('maneuver',s),'Icon','logo.png');
                    try


                        uitextarea(uifig,'position',[1,1,560,420],'Value',t);
                    catch asal=0;

                    end

                    %  R00=[app.XV(end) app.YV(end) app.ZV(end) ]';
                    %  V00=[app.v1(end) app.v2(end) app.v3(end)]';

                    % maneuvr();

            end


        end

        % Selection changed function: OrbitTree
        function OrbitTreeSelectionChanged(app, event)
            selectedNodes = app.OrbitTree.SelectedNodes;
            %             selectedNodes.Value
            switch selectedNodes.Text
                case 'Space Craft'
                    app22;
                case 'Ground Station'
                    perdict;
                    app.check=1;
                case 'Plot'
                    xyplot;
                case 'STK Plot'
                    STKplot
                    %                     app.OrbitTree.SelectedNode=app.Node
                    %                     app.OrbitTreeSelectionChanged(app.Node)

            end

        end

        % Selection changed function: OutputTree
        function OutputTreeSelectionChanged(app, event)
            selectedNodes = app.OutputTree.SelectedNodes;
            switch selectedNodes.Text
                case 'OrbitView'

                    file22
                case 'GroundTrack'
                case 'passTime'
                    pass

            end
        end

        % Context menu opening function: TargetsMenu
        function TargetsMenuOpening(app, event)

        end

        % Callback function
        function TextAreaValueChanged(app, event)
            value = app.TextArea.Value;
            z
        end

        % Callback function: PauseButton
        function PauseButtonClicked(app, event)
            app.calc=app.calc+1;
            %              app.num=2;
            app.c=2;
            set(app.PauseButton,'Enable','off')
            set(app.AnimationButton,'Enable','on')
        end

        % On callback: PauseButton
        function PauseButtonOn(app, event)

        end

        % Callback function: OrbitMaximun
        function OrbitMaximunClicked(app, event)
            app.x=1;
            app.ano=1;
            l=app.Panel.Position;
            set(app.h, 'position', [l(1) l(2)+30 l(3) l(4)-20]);
        end

        % Callback function: OrbitMinum
        function OrbitMinumClicked(app, event)
            l=app.OrbitPanel.Position;
            set(app.h, 'position', [l(1)+150 l(2)+30 l(3)+150 l(4)-4]);
        end

        % Callback function: GroundMaximum
        function GroundMaximumClicked(app, event)
            app.x=1;
            app.ano=0;
            l=app.Panel.Position;
            try


                set(app.f2, 'position', [l(1) l(2)+30 l(3) l(4)-20]);

            catch 1+1;

            end
        end

        % Callback function: GroundMinum
        function GroundMinumClicked(app, event)
            l=app.GroundPanel.Position;
            try


                set(app.f2, 'position', [l(1)+150 l(2)+50 l(3)+150 l(4)-50]);

            catch 1+1;

            end
        end

        % Callback function: FasterAnimation
        function FasterAnimationClicked(app, event)
            app.increase=app.increase+10;
            if(app.increase==200)
                set(app.FasterAnimation,'Enable','off')
            end
            if(app.increase>10)
                set(app.SlowerAnimation,'Enable','on')
            end

        end

        % Callback function: SlowerAnimation
        function SlowerAnimationClicked(app, event)

            if(app.increase==10)
                set(app.SlowerAnimation,'Enable','off');
            end
            if(app.increase<200)
                set(app.FasterAnimation,'Enable','on');
            end
            if(app.increase>10)

                app.increase=app.increase-10;
            end
        end

        % Callback function
        function Panel_2SizeChanged(app, event)
            position = app.OrbitPanel.Position;

        end

        % Callback function
        function PanelButtonDown(app, event)

        end

        % Callback function
        function PanelSizeChanged(app, event)
            position = app.Panel.Position;

        end

        % Menu selected function: atargetMenu
        function atargetMenuSelected(app, event)
            atarget ();
        end

        % Menu selected function: EtargetMenu
        function EtargetMenuSelected(app, event)
            Etarget ();
        end

        % Menu selected function: itargetMenu
        function itargetMenuSelected(app, event)
            itarget ();
        end

        % Menu selected function: etargetMenu
        function etargetMenuSelected(app, event)
            ectarget();
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            warning("off")
          addpath("Functions")
             addpath("Icons")
             addpath("mlapp")
             addpath("Graphics")
             addpath("Files")
            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create MainGUI and hide until all components are created
            app.MainGUI = uifigure('Visible', 'off');
            app.MainGUI.IntegerHandle = 'on';
            app.MainGUI.Color = [0.502 0.502 0.502];
            app.MainGUI.Position = [100 100 778 593];
            app.MainGUI.Name = 'ALI SMAD';
            app.MainGUI.Icon =  'logo.png';
            app.MainGUI.CloseRequestFcn = createCallbackFcn(app, @MainGUICloseRequest, true);
            app.MainGUI.SizeChangedFcn = createCallbackFcn(app, @MainGUISizeChanged, true);

            % Create MissionMainButtonsBar
            app.MissionMainButtonsBar = uitoolbar(app.MainGUI);

            % Create WriteScript
            app.WriteScript = uipushtool(app.MissionMainButtonsBar);
            app.WriteScript.Tooltip = {'add script'};
            app.WriteScript.Icon = 'NewScript.png';

            % Create OpenScript
            app.OpenScript = uipushtool(app.MissionMainButtonsBar);
            app.OpenScript.Tooltip = {'Open Mission'};
            app.OpenScript.Icon = 'OpenScript.png';

            % Create SaveChanges
            app.SaveChanges = uipushtool(app.MissionMainButtonsBar);
            app.SaveChanges.Tooltip = {'Save Mission'};
            app.SaveChanges.ClickedCallback = createCallbackFcn(app, @SaveChangesClicked, true);
            app.SaveChanges.Icon = 'SaveMissionOld.png';

            % Create Restart
            app.Restart = uipushtool(app.MissionMainButtonsBar);
            app.Restart.Tooltip = {'restart'};
            app.Restart.ClickedCallback = createCallbackFcn(app, @RestartClicked, true);
            app.Restart.Icon = 'Flat_restart_icon.svg.png';

            % Create RunButoon
            app.RunButoon = uipushtool(app.MissionMainButtonsBar);
            app.RunButoon.Tooltip = {'Run'};
            app.RunButoon.ClickedCallback = createCallbackFcn(app, @RunButoonClicked, true);
            app.RunButoon.Icon = 'RunCalculation.png';

            % Create AnimationButton
            app.AnimationButton = uipushtool(app.MissionMainButtonsBar);
            app.AnimationButton.Tooltip = {'Animation'};
            app.AnimationButton.ClickedCallback = createCallbackFcn(app, @AnimationButtonClicked, true);
            app.AnimationButton.Icon = 'RunAnimation.png';
            app.AnimationButton.Separator = 'on';

            % Create PauseButton
            app.PauseButton = uitoggletool(app.MissionMainButtonsBar);
            app.PauseButton.Tooltip = {'Pause'};
            app.PauseButton.ClickedCallback = createCallbackFcn(app, @PauseButtonClicked, true);
            app.PauseButton.Icon = 'PauseMission.png';
            app.PauseButton.OnCallback = createCallbackFcn(app, @PauseButtonOn, true);

            % Create FasterAnimation
            app.FasterAnimation = uipushtool(app.MissionMainButtonsBar);
            app.FasterAnimation.Tooltip = {'Faster Animation'};
            app.FasterAnimation.ClickedCallback = createCallbackFcn(app, @FasterAnimationClicked, true);
            app.FasterAnimation.Icon = 'FasterAnimation.png';

            % Create SlowerAnimation
            app.SlowerAnimation = uipushtool(app.MissionMainButtonsBar);
            app.SlowerAnimation.Tooltip = {'Slower Animation'};
            app.SlowerAnimation.ClickedCallback = createCallbackFcn(app, @SlowerAnimationClicked, true);
            app.SlowerAnimation.Icon = 'SlowerAnimation.png';

            % Create OrbitMaximun
            app.OrbitMaximun = uipushtool(app.MissionMainButtonsBar);
            app.OrbitMaximun.Tooltip = {'maximize  Orbit View'};
            app.OrbitMaximun.ClickedCallback = createCallbackFcn(app, @OrbitMaximunClicked, true);
            app.OrbitMaximun.Icon = 'download.jpg';
            app.OrbitMaximun.Separator = 'on';

            % Create OrbitMinum
            app.OrbitMinum = uipushtool(app.MissionMainButtonsBar);
            app.OrbitMinum.Tooltip = {'Mnimize Orbit View'};
            app.OrbitMinum.ClickedCallback = createCallbackFcn(app, @OrbitMinumClicked, true);
            app.OrbitMinum.Icon = 'downloa2d.png';

            % Create GroundMaximum
            app.GroundMaximum = uipushtool(app.MissionMainButtonsBar);
            app.GroundMaximum.Tooltip = {'Maximize Groun Track'};
            app.GroundMaximum.ClickedCallback = createCallbackFcn(app, @GroundMaximumClicked, true);
            app.GroundMaximum.Icon = 'images.png';

            % Create GroundMinum
            app.GroundMinum = uipushtool(app.MissionMainButtonsBar);
            app.GroundMinum.Tooltip = {'Minimize Ground Track'};
            app.GroundMinum.ClickedCallback = createCallbackFcn(app, @GroundMinumClicked, true);
            app.GroundMinum.Icon = 'imag2es.jpg';

            % Create PushHOLDTool4
            app.PushHOLDTool4 = uipushtool(app.MissionMainButtonsBar);
            app.PushHOLDTool4.Tooltip = {'Hold Last Orbit'};
            app.PushHOLDTool4.ClickedCallback = createCallbackFcn(app, @PushHOLDTool4Clicked, true);
            app.PushHOLDTool4.Icon = 'hold.png';
            app.PushHOLDTool4.Separator = 'on';

            % Create GUIGrild
            app.GUIGrild = uigridlayout(app.MainGUI);
            app.GUIGrild.ColumnWidth = {127, '2x'};
            app.GUIGrild.RowHeight = {'1x'};
            app.GUIGrild.Padding = [10 13 10 13];
            app.GUIGrild.BackgroundColor = [0 0 0];

            % Create MissionButtoms
            app.MissionButtoms = uitabgroup(app.GUIGrild);
            app.MissionButtoms.Layout.Row = 1;
            app.MissionButtoms.Layout.Column = 1;

            % Create orbitTab
            app.orbitTab = uitab(app.MissionButtoms);
            app.orbitTab.Title = 'orbit';

            % Create OrbitTree
            app.OrbitTree = uitree(app.orbitTab);
            app.OrbitTree.SelectionChangedFcn = createCallbackFcn(app, @OrbitTreeSelectionChanged, true);
            app.OrbitTree.Position = [2 2 122 542];

            % Create SpaceCraftNode
            app.SpaceCraftNode = uitreenode(app.OrbitTree);
            app.SpaceCraftNode.Icon = 'rt_Spacecraft.png';
            app.SpaceCraftNode.Text = 'Space Craft';

            % Create BurnsNode
            app.BurnsNode = uitreenode(app.OrbitTree);
            app.BurnsNode.Icon = 'rt_ImpulsiveBurn.png';
            app.BurnsNode.Text = 'Burns';

            % Create GroundTrackNode
            app.GroundTrackNode = uitreenode(app.OrbitTree);
            app.GroundTrackNode.Icon = 'rt_groundtrackplot.png';
            app.GroundTrackNode.Text = 'Ground Track';

            % Create GroundStationNode
            app.GroundStationNode = uitreenode(app.OrbitTree);
            app.GroundStationNode.Icon = 'rt_GroundStation.png';
            app.GroundStationNode.Text = 'Ground Station';

            % Create PlotNode
            app.PlotNode = uitreenode(app.OrbitTree);
            app.PlotNode.Icon = 'rt_XYPlot.png';
            app.PlotNode.Text = 'Plot';

            % Create STKPlotNode
            app.STKPlotNode = uitreenode(app.OrbitTree);
            app.STKPlotNode.Icon =  'STKG.png';
            app.STKPlotNode.Text = 'STK Plot';

            % Create EditNode
            app.EditNode = uitreenode(app.OrbitTree);
            app.EditNode.Text = 'Edit';

            % Create Node
            app.Node = uitreenode(app.OrbitTree);
            app.Node.Text = '';

            % Create missionsTab
            app.missionsTab = uitab(app.MissionButtoms);
            app.missionsTab.Title = 'missions';

            % Create MissionTreee
            app.MissionTreee = uitree(app.missionsTab);
            app.MissionTreee.SelectionChangedFcn = createCallbackFcn(app, @MissionTreeeSelectionChanged, true);
            app.MissionTreee.Editable = 'on';
            app.MissionTreee.Position = [0 4 124 540];

            % Create OrbitNode
            app.OrbitNode = uitreenode(app.MissionTreee);
            app.OrbitNode.Icon = 'propagateevent.png';
            app.OrbitNode.Text = 'Orbit';

            % Create OutputTab
            app.OutputTab = uitab(app.MissionButtoms);
            app.OutputTab.Title = 'Output';

            % Create OutputTree
            app.OutputTree = uitree(app.OutputTab);
            app.OutputTree.SelectionChangedFcn = createCallbackFcn(app, @OutputTreeSelectionChanged, true);
            app.OutputTree.Position = [1 2 122 542];

            % Create OrbitViewNode
            app.OrbitViewNode = uitreenode(app.OutputTree);
            app.OrbitViewNode.Icon = 'openglplot.png';
            app.OrbitViewNode.Text = 'OrbitView';

            % Create GroundTrackNode_2
            app.GroundTrackNode_2 = uitreenode(app.OutputTree);
            app.GroundTrackNode_2.Icon = 'rt_groundtrackplot.png';
            app.GroundTrackNode_2.Text = 'GroundTrack';

            % Create Panel
            app.Panel = uipanel(app.GUIGrild);
            app.Panel.ForegroundColor = [0.8 0.8 0.8];
            app.Panel.BackgroundColor = [0.502 0.502 0.502];
            app.Panel.Layout.Row = 1;
            app.Panel.Layout.Column = 2;

            % Create PanelGrild
            app.PanelGrild = uigridlayout(app.Panel);
            app.PanelGrild.ColumnWidth = {406, '1x', 100};
            app.PanelGrild.RowHeight = {'100x', '100x', 100};
            app.PanelGrild.RowSpacing = 2.25;
            app.PanelGrild.Padding = [10 2.25 10 2.25];
            app.PanelGrild.BackgroundColor = [0.502 0.502 0.502];

            % Create OrbitPanel
            app.OrbitPanel = uipanel(app.PanelGrild);
            app.OrbitPanel.Visible = 'off';
            app.OrbitPanel.Layout.Row = 1;
            app.OrbitPanel.Layout.Column = 1;

            % Create GroundPanel
            app.GroundPanel = uipanel(app.PanelGrild);
            app.GroundPanel.Visible = 'off';
            app.GroundPanel.Layout.Row = [2 3];
            app.GroundPanel.Layout.Column = 1;

            % Create Image
            app.Image = uiimage(app.PanelGrild);
            app.Image.Layout.Row = 3;
            app.Image.Layout.Column = 3;
            app.Image.ImageSource =  'logo.png';

            % Create TargetsMenu
            app.TargetsMenu = uicontextmenu(app.MainGUI);
            app.TargetsMenu.ContextMenuOpeningFcn = createCallbackFcn(app, @TargetsMenuOpening, true);

            % Create addmaneuverMenu
            app.addmaneuverMenu = uimenu(app.TargetsMenu);
            app.addmaneuverMenu.MenuSelectedFcn = createCallbackFcn(app, @addmaneuverMenuSelected, true);
            app.addmaneuverMenu.Text = 'add maneuver';

            % Create atargetMenu
            app.atargetMenu = uimenu(app.TargetsMenu);
            app.atargetMenu.MenuSelectedFcn = createCallbackFcn(app, @atargetMenuSelected, true);
            app.atargetMenu.Text = '''a'' target';

            % Create etargetMenu
            app.etargetMenu = uimenu(app.TargetsMenu);
            app.etargetMenu.MenuSelectedFcn = createCallbackFcn(app, @etargetMenuSelected, true);
            app.etargetMenu.Text = '''e'' target';

            % Create EtargetMenu
            app.EtargetMenu = uimenu(app.TargetsMenu);
            app.EtargetMenu.MenuSelectedFcn = createCallbackFcn(app, @EtargetMenuSelected, true);
            app.EtargetMenu.Text = '''E'' target';

            % Create itargetMenu
            app.itargetMenu = uimenu(app.TargetsMenu);
            app.itargetMenu.MenuSelectedFcn = createCallbackFcn(app, @itargetMenuSelected, true);
            app.itargetMenu.Text = '''i'' target';
            
            % Assign app.TargetsMenu
            app.MissionTreee.ContextMenu = app.TargetsMenu;
            app.OrbitNode.ContextMenu = app.TargetsMenu;

            % Show the figure after all components are created
            app.MainGUI.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SMAD_app

            % Create UIFigure and components
            createComponents(app)
            addpath("Functions")
            addpath("Icons")
            addpath("mlapp")
            addpath("Graphics")
            addpath("Files")

            % Register the app with App Designer
            registerApp(app, app.MainGUI)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.MainGUI)
        end
    end
end