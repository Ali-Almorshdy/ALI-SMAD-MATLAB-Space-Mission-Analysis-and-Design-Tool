function rot3d(huh)

if nargin<1
   set(gcf,'WindowButtonDownFcn','orbits(''down'')');
   set(gcf,'WindowButtonUpFcn','orbits(''up'')');
   set(gcf,'WindowButtonMotionFcn','');
else

   
   switch huh

   case 'zoom'
      rdata = get(gca,'userdata');
      oldpt = rdata.oldpt;
      newpt = get(0,'PointerLocation');
      dy = (newpt(2) - oldpt(2))/abs(oldpt(2));
      camzoom(gca,1+dy)
      rdata.oldpt = newpt;
      set(gca,'userdata',rdata)
   end
end
   end
