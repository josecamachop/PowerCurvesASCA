%%  A permutation test for the computiation of power curves at 0.05 
% 
% Experimental Design: X = m + A + B + AB + C(A) 
% where: 
%   - A: Treatment 
%   - B: Time 
%   - C(A): Subject, nested in Treatment
%   - AB: Interaction A & B
%
% Permute: Residuals from Reduced Model
% Tree Reduced Model: X = m + B + C(A) + AB
% Constrain: None
% Ex Units: Unit samples
% Statistic: SSE
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

function e = ResTreeRM(Xm,D,parglmo)
        
        n_perm = 1000;
        e = 0;
        
        % GLM model calibration with LS, only fixed factors
        
        D1 = D;
        D1(:,parglmo.factors{1}.Dvars) = [];
        B = pinv(D1'*D1)*D1'*Xm;
        X_residuals = Xm - D1*B;
        
        D2 = D(:,parglmo.factors{1}.Dvars);
        B = pinv(D2'*D2)*D2'*X_residuals;
        SSQ_class(1) = sum(sum((D2*B).^2));
        
        % Permutations
        for j = 1 : n_perm
            
            perms = randperm(size(X_residuals,1)); % permuted data (permute whole data matrix)
            
            B = pinv(D2'*D2)*D2'*X_residuals(perms, :);
            
            SSQ_class(1 + j) = sum(sum((D2*B).^2));
            
        end        % permutations
        
        p_value = (length(find(SSQ_class(2:n_perm + 1) >= SSQ_class(1))) + 1)/(n_perm);
        
        disp('p-value factor class')
        disp(p_value)
        
        if (p_value<0.05), e = 1; end
        
