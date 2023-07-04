%% Plot of figures: exp^3 distribution 
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 05/Jul/2023
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

%% TN.xlsx

load random_e1 e alpha

% Compare Interaction

i=1;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');  
xlabel('Alpha'),ylabel('Power'),title('Interaction Time x Class')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

% Compare time

i=2;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');   
xlabel('Alpha'),ylabel('Power'),title('Factor Time')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

% Compare class

i=3;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');   
xlabel('Alpha'),ylabel('Power'),title('Factor Class')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

% Compare individual

i=4;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');  
xlabel('Alpha'),ylabel('Power'),title('Factor Individual')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

%% LB.xlsx

load random_e2 e alpha

% Compare Interaction

i=1;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');  
xlabel('Alpha'),ylabel('Power'),title('Interaction Time x Class')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

% Compare time

i=2;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');   
xlabel('Alpha'),ylabel('Power'),title('Factor Time')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

% Compare class

i=3;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');   
xlabel('Alpha'),ylabel('Power'),title('Factor Class')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

% Compare individual

i=4;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');  
xlabel('Alpha'),ylabel('Power'),title('Factor Individual')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

%% TN_balanced

load random_e3 e alpha

% Compare Interaction

i=1;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');  
xlabel('Alpha'),ylabel('Power'),title('Interaction Time x Class')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

% Compare time

i=2;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');   
xlabel('Alpha'),ylabel('Power'),title('Factor Time')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

% Compare class

i=3;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');   
xlabel('Alpha'),ylabel('Power'),title('Factor Class')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')

% Compare individual

i=4;
figure, plot(alpha,e(1,:,i),'b'); hold on, plot(alpha,e(2,:,i),'r'); plot(alpha,e(3,:,i),'g.-'); plot(alpha,e(4,:,i),'k--'); plot(alpha,e(5,:,i),'m'); plot(alpha,e(6,:,i),'c');  
xlabel('Alpha'),ylabel('Power'),title('Factor Individual')
 legend('Obs TreeFM','Obs TreeFM 2PCs','Obs TreeFM III','Obs TreeFM A+B+AB','Obs TreeFM F','BoxCox Obs TreeFM F')