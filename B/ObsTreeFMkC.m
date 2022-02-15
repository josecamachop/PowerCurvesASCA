%% B permutation test for the computiation of power curves at 0.05 
% 
% Experimental Design: X = m + A + B + AB + C(A) 
% where: 
%   - A: Treatment 
%   - B: Time 
%   - C(A): Subject, nested in Treatment
%   - AB: Interaction A & B
%
% Permute: Raw samples Xm (given by the exact test)
% Tree Full Model: X = m + A + B + C(A) + AB (changes due to unballanceness)
% Constrain: We permute within C(A)
% Ex Units: unit samples
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

function e = ObsTreeFMkC(Xm,D,parglmo,vp)
        
        n_perm = 1000;
        e = 0;
      
        % GLM model calibration with LS, only fixed factors
        D1 = D;
        B = pinv(D1'*D1)*D1'*Xm;
        SSQ_time(1) = sum(sum((D1(:,parglmo.factors{2}.Dvars)*B(parglmo.factors{2}.Dvars,:)).^2));
        
        % Permutations
        up = unique(vp);
        for j = 1 : n_perm
            
            perms = zeros(size(Xm,1),1);
            for p = 1:length(up)
                ind = find(vp == up(p));
                perms(ind) = ind(randperm(length(ind)));
            end
            
            B = pinv(D1'*D1)*D1'*Xm(perms,:);
            
            SSQ_time(1 + j) = sum(sum((D1(:,parglmo.factors{2}.Dvars)*B(parglmo.factors{2}.Dvars,:)).^2));
            
        end    % permutations
        
        
        p_value = (length(find(SSQ_time(2:n_perm + 1) >= SSQ_time(1))) + 1)/(n_perm);
        
        disp('p-value factor time')
        disp(p_value)
        
        if (p_value<0.05), e = 1; end

        
