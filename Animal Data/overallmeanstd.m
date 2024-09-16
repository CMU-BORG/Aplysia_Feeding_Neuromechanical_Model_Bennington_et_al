function [alln, allmean, allstd] = overallmeanstd(grpn, grpmean, grpstd)
% [alln, allmean, allstd] = overallmeanstd(grpn, grpmean, grpstd)
% calculates the overall number of observations (alln), the mean (allmean)
% and standard deviation (allstd) of multiple groups with given number of 
% observations (grpn), mean values (grpmean) and standard deviations 
% (grpstd) for each group
% this is usefull if you don't have your observations stored anymore or
% you have a really high number of observations
% solution comes from ...
% http://www.burtonsys.com/climate/composite_standard_deviations.html
% 
% Example: copy this example (uncommented) to your command window and press enter
% grp1 = [50,43,45,56]
% grp2 = [65,67,50,81,83,90]
% grp3 = [43,77,66,84,94,54,67]
% grpn = [size(grp1,2);size(grp2,2);size(grp3,2)]
% grpmean = [mean(grp1);mean(grp2);mean(grp3)]
% grpstd = [std(grp1);std(grp2);std(grp3)]
% grps = [grp1,grp2,grp3];
% [alln, allmean, allstd] = overallmeanstd(grpn, grpmean, grpstd)
% n_ = size(grps,2)
% mean_ = mean(grps)
% std_ = std(grps)

% overall Number of observations
alln = sum(grpn);

% overall Mean
allmean = sum(grpmean.*grpn/alln);

% overall Error Sum of Squares
ESS = sum(grpstd.^2.*(grpn-1));

% overall Group Sum of Squares
GSS = sum((grpmean-allmean).^2.*grpn);

% overall Standard Deviation
allstd = sqrt((ESS+GSS)/(alln-1));
    
    
end