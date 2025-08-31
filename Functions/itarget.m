function itarget ()
     
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
         'String','i','style','text','foregroundcolor','w','fontsize',9,'backgroundcolor','r','tag','V');
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
            V=norm(V00)
             
            coe=ctok(R00,V00);
           
                    incl = coe(4)*180/pi
                 chang=abs(x-incl);
                 vn=V*sind(chang)
                 vv=V*(1-cosd(chang))
                    
                  
         manev(-vv,0,vn,S);
 solver=sprintf(' delta V in velDirection = %d and in normal direction =%d',-vv,vn);
(msgbox(solver));
         
            delete(findobj('tag','bu1'))
            delete(findobj('tag','a'))
            delete(findobj('tag','V'))
            delete(findobj('tag','MS'))
            delete(findobj('tag','S'))
    end
end
      

