%% Plot of figures in the paper for A.
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

f=figure; plot(alpha,e(1,:),'b-s'); hold on, plot(alpha,e(2,:),'b-s','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(3,:),'b--s'); plot(alpha,e(11,:),'b--s','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(12,:),'b:s','MarkerSize',10); 
plot(alpha,e(5,:),'r-o');   plot(alpha,e(6,:),'r-o','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(7,:),'r--o'); 
plot(alpha,e2(1,:,3),'g-d');  plot(alpha,e2(2,:,3),'g-d','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(14,:),'g--d','MarkerFaceColor',[0.5,0.5,0.5]); 
xlabel('Effect Size'),ylabel('Power'),title('Treatment Factor')
a = axis; a(1:2) = [0 0.5]; axis(a);

legend('Obs TreeFM k(B) eu(AB)', 'Obs TreeFM k(B) eu(C(A))', 'Obs BranchFM k(B) eu(C(A))', 'Obs minFM k(B) eu(C(A)))', 'Obs TreeFM k(B) eu(C(A)) F',... % pseudo-exact tests
    'Res TreeRM','Res BranchRM','Res minRM', ... % % residuals RM
    'Obs TreeFM','Obs TreeFM 2PCs', 'Obs TreeFM F', 'Location','southeast')  % unrestricted raw

a=get(f,'Children');
set(a(1),'Fontsize',10);
set(a(2),'Fontsize',12);
saveas(gcf,'../Fig/A_TN');
saveas(gcf,'../Fig/A_TN.eps','epsc');


load ../random_e2 e
e2 = e;
load random_e2 e

f=figure; plot(alpha,e(1,:),'b-s'); hold on, plot(alpha,e(2,:),'b-s','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(3,:),'b--s'); plot(alpha,e(11,:),'b--s','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(12,:),'b:s','MarkerSize',10); 
plot(alpha,e(5,:),'r-o');   plot(alpha,e(6,:),'r-o','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(7,:),'r--o'); 
plot(alpha,e2(1,:,3),'g-d');  plot(alpha,e2(2,:,3),'g-d','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(14,:),'g--d','MarkerFaceColor',[0.5,0.5,0.5]); 
xlabel('Effect Size'),ylabel('Power')
a = axis; a(1:2) = [0 0.5]; axis(a);

% legend('Obs TreeFM k(B) eu(AB)', 'Obs TreeFM k(B) eu(C(A))', 'Obs BranchFM k(B) eu(C(A))', 'Obs minFM k(B) eu(C(A)))', 'Obs TreeFM k(B) eu(C(A)) F',... % pseudo-exact tests
%    'Res TreeRM','Res BranchRM','Res minRM', ... % % residuals RM
%    'Obs TreeFM','Obs TreeFM 2PCs', 'Obs TreeFM F', 'Location','southeast')  % unrestricted raw

a=get(f,'Children');
%set(a(1),'Fontsize',12);
%set(a(2),'Fontsize',16);
set(a,'Fontsize',18);
saveas(gcf,'../Fig/A_LB');
saveas(gcf,'../Fig/A_LB.eps','epsc');


load ../random_e3 e
e2 = e;
load random_e3 e

f=figure; plot(alpha,e(1,:),'b-s'); hold on, plot(alpha,e(2,:),'b-s','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(3,:),'b--s'); plot(alpha,e(11,:),'b--s','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(12,:),'b:s','MarkerSize',10); 
plot(alpha,e(5,:),'r-o');   plot(alpha,e(6,:),'r-o','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(7,:),'r--o'); 
plot(alpha,e2(1,:,3),'g-d');  plot(alpha,e2(2,:,3),'g-d','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(14,:),'g--d','MarkerFaceColor',[0.5,0.5,0.5]); 
xlabel('Effect Size'),ylabel('Power')
a = axis; a(1:2) = [0 0.5]; axis(a);

% legend('Obs TreeFM k(B) eu(AB)', 'Obs TreeFM k(B) eu(C(A))', 'Obs BranchFM k(B) eu(C(A))', 'Obs minFM k(B) eu(C(A)))', 'Obs TreeFM k(B) eu(C(A)) F',... % pseudo-exact tests
%    'Res TreeRM','Res BranchRM','Res minRM', ... % % residuals RM
%    'Obs TreeFM','Obs TreeFM 2PCs', 'Obs TreeFM F', 'Location','southeast')  % unrestricted raw

a=get(f,'Children');
%set(a(1),'Fontsize',12);
%set(a(2),'Fontsize',16);
set(a,'Fontsize',18);
saveas(gcf,'../Fig/A_TNb');
saveas(gcf,'../Fig/A_TNb.eps','epsc');

