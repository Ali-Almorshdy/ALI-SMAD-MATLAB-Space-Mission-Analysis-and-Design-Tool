function atarget ()
     
        values = {'0.2574','0','1.9550'};
        strings2 = {'Acheive'};
        tags = {'a','MN','MB','MS'};
        callbacks={'manav(''manauver'')'};
       
%         if isempty(findobj('tag','orbits'))
%             www = figure('tag','orbits');
%         else
%             www = findobj('tag','orbits');
%             figure(www);
%             clf
%         end
%         set(www,'position',[25   75   600   500],'color','b');
        
        
        % uicontrol('Units','pixels','Position',[100 400 20 20],...
        %          'String','V','style','text','foregroundcolor','w','fontsize',9,'backgroundcolor','none');
%         l=uicontrol('Units','pixels','Position',[55 (20*(1-1)+16) 50 20],...
%             'tag',tags{1},'style','edit','string',values{1},'backgroundcolor',[1 1 1]);
        %                 uicontrol('Units','pixels','Position',[170 400 20 20],...
        %          'String','N','style','text','foregroundcolor','w','fontsize',9,'backgroundcolor','none');
        %                uicontrol('Units','pixels','Position',[190 400 50 20],...
        %          'tag','MN','style','edit','string','0','backgroundcolor',[1 1 1]);
        %                uicontrol('Units','pixels','Position',[240 400 20 20],...
        %          'String','B','style','text','foregroundcolor','w','fontsize',9,'backgroundcolor','none');
        %                uicontrol('Units','pixels','Position',[260 400 50 20],...
        %          'tag','MB','style','edit','string','0','backgroundcolor',[1 1 1]);
%         btn=uicontrol('Units','pixels','Position',[10 (25*(1-1)+150) 70 22],...
%             'string',strings2{1},'callback',@manavv,'backgroundcolor','red','tag','bu1');
        
        uicontrol('Units','pixels','Position',[100 400 20 20],...
         'String','a','style','text','foregroundcolor','w','fontsize',9,'backgroundcolor','r','tag','V');
               uicontrol('Units','pixels','Position',[120 400 50 20],...
         'tag',tags{1},'style','edit','string','0','backgroundcolor',[1 1 1]);
               
          s=uicontrol('Units','pixels','Position',[170 400 20 20],...
         'String','s','style','text','foregroundcolor','w','fontsize',9,'backgroundcolor','r','tag','S');
     
     uicontrol('Units','pixels','Position',[190 400 50 20],...
         'tag','MS','style','edit','string','1','backgroundcolor',[1 1 1]);
    
             btn=uicontrol('Units','pixels','Position',[240 400 50 20],...
            'string',strings2{1},'callback',@manavv,'backgroundcolor','red','tag','bu1');
        
        
        
%         clear all
    function manavv(btn,~)
        x = str2num(get(findobj('tag','a'),'string'));
              file='end.txt';
              S=str2num(get(findobj('tag','MS'),'string'));
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
            
             
            coe=ctok(R00,V00);
            e = coe(2);
                    RA = coe(3)*180/pi;
                    incl = coe(4)*180/pi;
                    w = coe(5)*180/pi;
                    TA0 = coe(6)*180/pi;
                    a = coe(7)';
                    
                     Pp=a*(1-e*e);
            mu = 3.986012E5;
            h=sqrt(Pp*mu);
            rp = (h^2/mu) * (1/(1 + e*cos(TA0*pi/180))) * (cos(TA0*pi/180)*[1;0;0] + sin(TA0*pi/180)*[0;1;0]);
            vp = (mu/h) * (-sin(TA0*pi/180)*[1;0;0] + (e + cos(TA0*pi/180))*[0;1;0]);
            dV=0;
            syms dV
             
            V=sym(vp);
            v1=sym(vp);
%             xx=0;
%             syms xx
%           oo=sqrt(V(1)^2+V(2)^2)-sqrt(V(1)^2+(V(2)+(dV+xx))^2)==-dV;
          
          k=vp(2);
             l=vp(1);
             j=dV;
            ince=dV+((2*j*(k^2 + l^2)^(1/2) + j^2 + k^2)^(1/2) - k - j);
            ince2=(- j - k - (2*j*(k^2 + l^2)^(1/2) + j^2 + k^2)^(1/2))+dV;
%           ince=solve(oo,xx)
         V(2)=V(2)+ince;
         v1(2)=v1(2)+ince2;
         
         
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
            v = Q_pX*V;
             v2 = Q_pX*v1;
            r = r';
            v = v';
            v2=v2';
            V0=v';
            V02=v2';
            R0=r';
%             X0=[R00;V0];
            E = norm(V0)^2/2 - mu/norm(R00);
             value= -mu/(2*E)==x;
            a1= double(vpasolve(value,dV));
         manev(a1,0,0,S);
 solver=sprintf(' delta V in velDirection = %d',a1);
(msgbox(solver));
         
            delete(findobj('tag','bu1'))
            delete(findobj('tag','a'))
            delete(findobj('tag','V'))
            delete(findobj('tag','MS'))
            delete(findobj('tag','S'))
    end
end
      

