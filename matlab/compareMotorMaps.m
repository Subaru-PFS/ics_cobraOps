for x = [3 5 7 9 11 13 15 17 19 21 23 25 27]
   figure(1)
   [oy ox] = getMMap(xml, x,'FAST',2,'rev');
   [cy cx] = getMMap(currentXML, x,'FAST',2,'rev');
   plot(oy(3:end)); 
   hold on;
   plot(cy(3:end),'r');
   keyboard;
   hold off
end