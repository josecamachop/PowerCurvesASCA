%% Computiation of Power Curves for AB
% 
% Experimental Design: X = m + A + B + AB + C(A) 
% where: 
%   - A: Treatment 
%   - B: Time 
%   - C(A): Subject, nested in Treatment
%   - AB: Interaction A & B
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


%% Inicialization

clear
close all
clc

tic

data = importdata('../TN.xlsx');
%data = importdata('../LB.xlsx');
%data = importdata('../TN_balanced.xlsx');

X = data.data;
obs_l = data.textdata(2:end,1);
class = data.textdata(2:end,2);
time = X(2:end,1);
for i=1:length(X(1,:))
    var_l{i} = num2str(X(1,i));
end
var_l(1) = [];
pac_l = obs_l;
for i=1:length(pac_l)
    pac_l{i} = pac_l{i}(1:5);
end

% vc: factor class
uc = unique(class);
vc = zeros(length(class),1);
for i=1:length(uc)
    vc(find(ismember(class,uc(i)))) = i;
end

% vt: factor time
vt = time;
ut = unique(vt);

% vp: factor person
up = unique(pac_l);
vp = zeros(length(pac_l),1);
for i=1:length(up)
    vp(find(ismember(pac_l,up(i)))) = i;
end

% Design matrix
F = [vc vt vp];
interactions = [1 2];
n_factors = 3;
n_interactions = 1;

% Create Design Matrix
n = 1;
D = ones(size(X,1)-1,1);

for f = 1 : n_factors
    uF = unique(F(:,f));
    for i = 1:length(uF)-1
        D(find(F(:,f)==uF(i)),n+i) = 1;
    end
    parglmo.factors{f}.Dvars = n+(1:length(uF)-1);
    D(find(F(:,f)==uF(end)),parglmo.factors{f}.Dvars) = -1;
    n = n + length(uF) - 1;
end

for i = 1 : n_interactions
    for j = parglmo.factors{interactions(i,1)}.Dvars
        for k = parglmo.factors{interactions(i,2)}.Dvars
            D(:,end+1) = D(:,j).* D(:,k);
        end
    end
    parglmo.interactions{i}.Dvars = n+1:size(D,2);
    n = size(D,2);
end

perm_tot = 100;

%% Repite permutations

alpha = 0:0.05:0.5; % This controls the compromise of true significance vs random
eD = zeros(9,length(alpha),perm_tot);

for i2=1:perm_tot
    
    disp(i2)
    
    rng(i2);
    
    % Build data
    
    Xpac = randn(length(up),length(var_l));
    Xpac = Xpac/norm(Xpac);
    Xtime = randn(length(ut),length(var_l));
    Xtime = Xtime/norm(Xtime);
    for i = 1:length(obs_l)
        if strcmp(class{i},'R')
            Xstruct(i,:) = Xpac(vp(i),:) + Xtime(vt(i),:);
        else
            Xstruct(i,:) = Xpac(vp(i),:);
        end
    end
    
    Xnoise = randn(length(obs_l),length(var_l));
    Xnoise = Xnoise/norm(Xnoise);

    for a = 1:length(alpha)
        
        % Make a blend
         
        Xm = alpha(a)*Xstruct + (1-alpha(a))*Xnoise; 
   
        eD(1,a,i2) = ResTreeRM(Xm,D,parglmo);
        eD(2,a,i2) = ResBranchRM(Xm,D,parglmo);
        eD(3,a,i2) = ResTreeRMF(Xm,D,parglmo);
        eD(4,a,i2) = ObsBranchFM(Xm,D,parglmo);
        eD(5,a,i2) = ObsBranchFM2PCs(Xm,D,parglmo);
        eD(6,a,i2) = ResTreeRMIII(Xm,D,parglmo);
        eD(7,a,i2) = ObsTreeFMF(Xm,D,parglmo);
        eD(8,a,i2) = ObsTreeFMIII(Xm,D,parglmo);
        eD(9,a,i2) = ObsTreeFMFIII(Xm,D,parglmo);        
    end
end

toc

e = mean(eD,3);
eST = std(eD,0,3);
save random_e1 e eST eD alpha
%save random_e2 e eST eD alpha
%save random_e3 e eST eD alpha

%% Compare 
load ../random_e1 e
e2 = e;
load random_e1 e

f=figure; hold on, box on  
plot(alpha,e(1,:),'r-o');   plot(alpha,e(2,:),'r-o','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(3,:),'r--o'); plot(alpha,e(6,:),'r--o','MarkerFaceColor',[0.5,0.5,0.5]); 
plot(alpha,e2(1,:,1),'g-d');  plot(alpha,e2(2,:,1),'g-d','MarkerFaceColor',[0.5,0.5,0.5]);  plot(alpha,e(4,:),'g--d','MarkerSize',10); 
plot(alpha,e(7,:),'g--d','MarkerFaceColor',[0.5,0.5,0.5]); plot(alpha,e(8,:),'g:d','MarkerSize',10); plot(alpha,e(9,:),'g:d','MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('Effect Size'),ylabel('Power'),title('Time x Treatment Interaction')
a = axis; a(1:2) = [0 0.6]; axis(a);

legend('Res TreeRM','Res BranchRM', 'Res TreeRM F', 'Res TreeRM III', ... % residuals RM
    'Obs TreeFM','Obs TreeFM 2PCs','Obs BranchFM', ...
    'Obs TreeFM F', 'Obs TreeFM III', 'Obs TreeFM F III', 'Location','southeast')  % unrestricted raw
