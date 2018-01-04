function output = stow_cobras(targetlist,bench)
% figure out the stow angles for cobras in a rail
% want elbows to point at neighbors and fibers to be inboard.
% replace the targetlist with stowed triangles.

   % define the angles for the elbow to point in, minimizing chance of contact
   thtLft = mean(angle(diff(bench.center(1:2:end)))) + pi/2;
   thtRgt = thtLft - pi;

   phi = 1 - pi;
   
   tht(1:2:57) = thtLft;
   tht(2:2:56) = thtRgt;
   tht = tht(:);
   
   stowxy = bench.center + bench.L1 .* exp(i*tht) + bench.L2 .* exp(i*(tht + phi));

   tgt1 = stowxy;
   tgt2 = stowxy;
   tgt3 = stowxy;
   
   tgt1(1:3:end) = targetlist(1:3:end);
   tgt2(2:3:end) = targetlist(2:3:end);
   tgt3(3:3:end) = targetlist(3:3:end);

   
   output = packstruct(tgt1,tgt2,tgt3);
   output.stow = stowxy;