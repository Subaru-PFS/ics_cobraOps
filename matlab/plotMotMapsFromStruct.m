function output = plotMotMapsFromStruct(mapStructure,varargin)

J1FwR = mMapStr2NumArr(mapStructure.Joint1_fwd_regions.Text);
J1RvR = mMapStr2NumArr(mapStructure.Joint1_rev_regions.Text);
J2FwR = mMapStr2NumArr(mapStructure.Joint2_fwd_regions.Text);
J2RvR = mMapStr2NumArr(mapStructure.Joint2_rev_regions.Text);

J1FwS = mMapStr2NumArr(mapStructure.Joint1_fwd_stepsizes.Text);
J1RvS = mMapStr2NumArr(mapStructure.Joint1_rev_stepsizes.Text);
J2FwS = mMapStr2NumArr(mapStructure.Joint2_fwd_stepsizes.Text);
J2RvS = mMapStr2NumArr(mapStructure.Joint2_rev_stepsizes.Text);

subplot(2,2,1);hold on
output.ph(1) = plot(J1FwR(3:2+J1FwR(1)),J1FwS(3:2+J1FwR(1)),varargin{:});
title('Theta Forward')
ylabel('deg/step')
xlabel('stage angle')

subplot(2,2,2);hold on
output.ph(2) = plot(J1RvR(3:2+J1RvR(1)),J1RvS(3:2+J1RvR(1)),varargin{:});
title('Theta Reverse')
ylabel('deg/step')
xlabel('stage angle')

subplot(2,2,3);hold on
output.ph(3) = plot(J2FwR(3:2+J2FwR(1)),J2FwS(3:2+J2FwR(1)),varargin{:});
title('Phi Forward')
ylabel('deg/step')
xlabel('stage angle')

subplot(2,2,4);hold on
output.ph(4) = plot(J2RvR(3:2+J2RvR(1)),J2RvS(3:2+J2RvR(1)),varargin{:});
title('Phi Reverse')
ylabel('deg/step')
xlabel('stage angle')

end