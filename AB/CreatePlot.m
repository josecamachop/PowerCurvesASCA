%% Plot of figures in the paper for AB.
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

load ../random_e1 e
e2 = e;
load random_e1 e alpha

f=figure; hold on, box on  
plot(alpha,e(1,:),'r-o');   plot(alpha,e(2,:),'r-o','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(3,:),'r--o'); plot(alpha,e(6,:),'r--o','MarkerFaceColor',[0.5,0.5,0.5]); 
plot(alpha,e2(1,:,1),'g-d');  plot(alpha,e2(2,:,1),'g-d','MarkerFaceColor',[0.5,0.5,0.5]);  plot(alpha,e(4,:),'g--d','MarkerSize',10); 
plot(alpha,e(7,:),'g--d','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(8,:),'g:d','MarkerSize',10); plot(alpha,e(9,:),'g:d','MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('Effect Size'),ylabel('Power'),title('Time x Treatment Interaction')
a = axis; a(1:2) = [0 0.5]; axis(a);

legend('Res TreeRM','Res BranchRM', 'Res TreeRM F', 'Res TreeRM III', ... % residuals RM
    'Obs TreeFM','Obs TreeFM 2PCs','Obs BranchFM', ...
    'Obs TreeFM F', 'Obs TreeFM III', 'Obs TreeFM F III', 'Location','southeast')  % unrestricted raw

a=get(f,'Children');
set(a(1),'Fontsize',12);
set(a(2),'Fontsize',12);
saveas(gcf,'../Fig/AB_TN');
saveas(gcf,'../Fig/AB_TN.eps','epsc');


load ../random_e2 e
e2 = e;
load random_e2 e

f=figure; hold on, box on   
plot(alpha,e(1,:),'r-o');   plot(alpha,e(2,:),'r-o','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(3,:),'r--o'); plot(alpha,e(6,:),'r--o','MarkerFaceColor',[0.5,0.5,0.5]); 
plot(alpha,e2(1,:,1),'g-d');  plot(alpha,e2(2,:,1),'g-d','MarkerFaceColor',[0.5,0.5,0.5]);  plot(alpha,e(4,:),'g--d','MarkerSize',10); 
plot(alpha,e(7,:),'g--d','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(8,:),'g:d','MarkerSize',10); plot(alpha,e(9,:),'g:d','MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('Effect Size'),ylabel('Power'),%title('Time x Treatment Interaction')
a = axis; a(1:2) = [0 0.5]; axis(a);

%legend('Res TreeRM','Res BranchRM', 'Res TreeRM F', 'Res TreeRM III', ... % residuals RM
%    'Obs TreeFM','Obs TreeFM 2PCs','Obs BranchFM', ...
%    'Obs TreeFM F', 'Obs TreeFM III', 'Obs TreeFM F III', 'Location','southeast')  % unrestricted raw

a=get(f,'Children');
%set(a(1),'Fontsize',12);
%set(a(2),'Fontsize',16);
set(a,'Fontsize',18);
saveas(gcf,'../Fig/AB_LB');
saveas(gcf,'../Fig/AB_LB.eps','epsc');

load ../random_e3 e
e2 = e;
load random_e3 e

f=figure; hold on, box on   
plot(alpha,e(1,:),'r-o');   plot(alpha,e(2,:),'r-o','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(3,:),'r--o'); plot(alpha,e(6,:),'r--o','MarkerFaceColor',[0.5,0.5,0.5]); 
plot(alpha,e2(1,:,1),'g-d');  plot(alpha,e2(2,:,1),'g-d','MarkerFaceColor',[0.5,0.5,0.5]);  plot(alpha,e(4,:),'g--d','MarkerSize',10); 
plot(alpha,e(7,:),'g--d','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(8,:),'g:d','MarkerSize',10); plot(alpha,e(9,:),'g:d','MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('Effect Size'),ylabel('Power'),%title('Time x Treatment Interaction')
a = axis; a(1:2) = [0 0.5]; axis(a);

%legend('Res TreeRM','Res BranchRM', 'Res TreeRM F', 'Res TreeRM III', ... % residuals RM
%    'Obs TreeFM','Obs TreeFM 2PCs','Obs BranchFM', ...
%    'Obs TreeFM F', 'Obs TreeFM III', 'Obs TreeFM F III', 'Location','southeast')  % unrestricted raw

a=get(f,'Children');
%set(a(1),'Fontsize',12);
%set(a(2),'Fontsize',16);
set(a,'Fontsize',18);
saveas(gcf,'../Fig/AB_TNb');
saveas(gcf,'../Fig/AB_TNb.eps','epsc');