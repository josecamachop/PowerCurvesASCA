%% Plot of the bootstrap figures in the paper for B.
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 09/Feb/2022
%
% Copyright (C) 2022  University of Granada, Granada
% Copyright (C) 2022  Jose Camacho Paez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

clear

load ../random_e1 eD
eD2 = eD;
load random_e1 eD 

label = { 'Obs TreeFM k(C(A))', 'Obs BranchFM k(C(A))', 'Obs minFM k(C(A))', 'Obs TreeFM k(A) eu(AB) F', 'Res BranchRM','Res minRM','Res minRM eu(AB)','Res minRM F','Obs TreeFM','Obs TreeFM 2PCs','Obs BranchFM','Obs TreeFM F'};

f=figure;
Alpha0 = [squeeze(eD([2:3 13 15 6:9],1,:));squeeze(eD2(1:2,1,2,:));squeeze(eD([11 16],1,:))]';
j=1;
for i=1:1000 % bootstrap
    A0(j,:) = mean(Alpha0(ceil(100*rand(100,1)),:));
    j = j+1;
end
    
aboxplot(A0,'OutlierMarker','+','OutlierMarkerSize',6), ylabel('Power','FontSize',24), %ax = axis; axis([ax(1:2) 0 600])
a=get(f,'Children');
set(a(1),'FontSize',16);
set(a(1),'XTickLabel',label);
set(a(1),'XTickLabelRotation',45);
saveas(gcf,'../Fig/B_TNboots');
saveas(gcf,'../Fig/B_TNboots.eps','epsc');