%% Computiation of Power Curves for unconstrained permutations on raw observations
% 
% Experimental Design: X = m + A + B + AB + C(A) 
% where: 
%   - A: Treatment 
%   - B: Time 
%   - C(A): Subject, nested in Treatment
%   - AB: Interaction A & B
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 05/Jun/2023
%
% Copyright (C) 2023  University of Granada, Granada
% Copyright (C) 2023  Jose Camacho Paez
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

%data = importdata('TN.xlsx');
%data = importdata('LB.xlsx');
data = importdata('TN_balanced.xlsx');

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

perm_tot = 100;

%% Repite permutations

alpha = 0:0.05:0.5; % This controls the compromise of true significance vs random
eD = zeros(5,length(alpha),4,perm_tot);

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
        
        % raw data perm
        [T, paranovao] = parglm(Xm,[vc vt vp],[1 2],0,[],0);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(1,a,o,i2) = 1;
            end
        end
        
        % raw data perm 2PCs
        paranovao = parglm2(Xm,[vc vt vp],[1 2],0);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(2,a,o,i2) = 1;
            end
        end
        
        % raw data perm ETIII
        paranovao = parglmIII(Xm,[vc vt vp],[1 2],0);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(3,a,o,i2) = 1;
            end
        end
        
        % raw data perm ETIII
        paranovao = parglmIII(Xm,[vc vt],[1 2],0);
        reo = [3 2 1];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(4,a,o,i2) = 1;
            end
        end
     
        % raw data perm F
        [T, paranovao] = parglm(Xm,[vc vt vp],[1 2],0,[],1);
        reo = [4 2 1 3];
        for o = 1:length(reo)
            if paranovao.p(reo(o))<0.05
                eD(5,a,o,i2) = 1;
            end
        end
        
    end
end

toc

e = mean(eD,4);
eST = std(eD,0,4);
%save random_e1 e eST eD alpha
%save random_e2 e eST eD alpha
save random_e3 e eST eD alpha

%% Compare Interaction

i=1;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); 
xlabel('Alpha'),ylabel('Power'),title('Interaction Time x Class')
legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F')

%% Compare time

i=2;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m');  
xlabel('Alpha'),ylabel('Power'),title('Factor Time')
legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F')

%% Compare class

i=3;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m');  
xlabel('Alpha'),ylabel('Power'),title('Factor Class')
legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F')

%% Compare individual

i=4;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); 
xlabel('Alpha'),ylabel('Power'),title('Factor Individual')
legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F')

