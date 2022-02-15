%% C(A) permutation test for the computiation of power curves at 0.05 
% 
% Experimental Design: X = m + A + B + AB + C(A) 
% where: 
%   - A: Treatment 
%   - B: Time 
%   - C(A): Subject, nested in Treatment
%   - AB: Interaction A & B
%
% Permute: Raw samples Xm (given by the exact test)
% Model: X = m + A + B + C(A) + AB (changes due to unballanceness)
% Constrain: We permute within A and B (basically within AB)
% Ex Units: unit samples (given by the exact test)
% Statistic: F-ratio
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

function e = ObsTreeFMkABF(Xm,D,parglmo,vt,vc)
        
        n_perm = 1000;
        e = 0;
        
        dfA = size(unique(D(:,parglmo.factors{1}.Dvars),'Rows'),1) -1;
        dfB = size(unique(D(:,parglmo.factors{2}.Dvars),'Rows'),1) -1;
        dfAB = (size(unique(D(:,parglmo.factors{1}.Dvars),'Rows'),1)*size(unique(D(:,parglmo.factors{2}.Dvars),'Rows'),1) -1);
        dfC = size(unique(D(:,parglmo.factors{1}.Dvars),'Rows'),1)*(size(unique(D(:,parglmo.factors{3}.Dvars),'Rows'),1) -1);
       
        % GLM model calibration with LS, only fixed factors
        D1 = D;
        B = pinv(D1'*D1)*D1'*Xm;
        Xr = Xm - D1*B;
        
        Tdf = size(Xm,1)-1;    
        Rdf = Tdf - dfA - dfB - dfAB - dfC;
        
        SSQ_ind(1) = (sum(sum((D1(:,parglmo.factors{3}.Dvars)*B(parglmo.factors{3}.Dvars,:)).^2))/dfC) / (sum(sum((Xr).^2))/Rdf);
        
        % Permutations
        uc = unique(vc);
        ut = unique(vt);
        for j = 1 : n_perm     
            
            perms = zeros(size(Xm,1),1);
            for c = 1:length(uc)
                for t = 1:length(ut)
                    ind = find(vc == uc(c) & vt == ut(t));
                    perms(ind) = ind(randperm(length(ind)));
                end
            end
            
            B = pinv(D1'*D1)*D1'*Xm(perms,:);
            Xr = Xm(perms,:) - D1*B;
            
            SSQ_ind(1 + j) = (sum(sum((D1(:,parglmo.factors{3}.Dvars)*B(parglmo.factors{3}.Dvars,:)).^2))/dfC) / (sum(sum((Xr).^2))/Rdf);
            
        end        % permutations
        
        
        p_value = (length(find(SSQ_ind(2:n_perm + 1) >= SSQ_ind(1))) + 1)/(n_perm);
        
        disp('p-value factor individual')
        disp(p_value)
        
        if (p_value<0.05), e = 1; end

        
