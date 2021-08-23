function [correctionfit1,correctionfit2,correctionfit3,correctionfit4] = expfit4(linefluorescenceseries1,timeline1,timelinefull1,linefluorescenceseries2,timeline2,timelinefull2,linefluorescenceseries3,timeline3,timelinefull3,linefluorescenceseries4,timeline4,timelinefull4)
global ff1 ff2 ff3 ff4
ff1=fit(timeline1,linefluorescenceseries1,'exp2');
correctionfit1=ff1.a.*exp(ff1.b*timelinefull1)+ff1.c.*exp(ff1.d*timelinefull1);
ff2=fit(timeline2,linefluorescenceseries2,'exp2');
correctionfit2=ff2.a.*exp(ff2.b*timelinefull2)+ff2.c.*exp(ff2.d*timelinefull2);
ff3=fit(timeline3,linefluorescenceseries3,'exp2');
correctionfit3=ff3.a.*exp(ff3.b*timelinefull3)+ff3.c.*exp(ff3.d*timelinefull3);
ff4=fit(timeline4,linefluorescenceseries4,'exp2');
correctionfit4=ff4.a.*exp(ff4.b*timelinefull4)+ff4.c.*exp(ff4.d*timelinefull4);
end
