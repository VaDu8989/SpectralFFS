function [linefluorescenceseries_corrected] =depletioncorrectionfunc(timeline,linefluorescenceseries,f)
% Depletion correction of time trace using double-exponential and correction formula from Ries, Chiantia, Schwille, BJ, 2006;
% Exponential
correctionfit=f.a.*exp(f.b*timeline)+f.c.*exp(f.d*timeline);
size(correctionfit)
size(linefluorescenceseries)
linefluorescenceseries_corrected=linefluorescenceseries./sqrt(correctionfit./correctionfit(1))+correctionfit(1)*(1-sqrt(correctionfit./correctionfit(1)));

