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
% Tree Full Model: X = m + A + B + AB (changes due to unballanceness)
% Constrain: We permute within A
% Ex Units: AB cells (given by the exact test)
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

function e = ObsTreeFMkCF(Xm,D,parglmo,vt,vc)
        
        n_perm = 1000;
        e = 0;
        
        dfA = size(unique(D(:,parglmo.factors{1}.Dvars),'Rows'),1) -1;
        dfB = size(unique(D(:,parglmo.factors{2}.Dvars),'Rows'),1) -1;
        dfAB = (size(unique(D(:,parglmo.factors{1}.Dvars),'Rows'),1)*size(unique(D(:,parglmo.factors{2}.Dvars),'Rows'),1) -1);
        dfC = size(unique(D(:,parglmo.factors{1}.Dvars),'Rows'),1)*(size(unique(D(:,parglmo.factors{3}.Dvars),'Rows'),1) -1);
        
        % GLM model calibration with LS, only fixed factors
        D1 = D;
        B = pinv(D1'*D1)*D1'*Xm;
        
        SSQ_time(1) = (sum(sum((D1(:,parglmo.factors{2}.Dvars)*B(parglmo.factors{2}.Dvars,:)).^2))/dfB) / (sum(sum((D1(:,parglmo.interactions{1}.Dvars)*B(parglmo.interactions{1}.Dvars,:)).^2))/dfAB);
        
        % Permutations: This strategy of permitation (the Ex Units) does
        % not make sense if we include C(A) in the model
        % but this just doesn't work even for a branch or min FM
        D1p = D1;
        uc = unique(vc);
        ut = unique(vt);
        D1u = zeros(length(uc),length(ut),size(D1,2));
        eu = 1;
        for c = 1:length(uc)
            for t = 1:length(ut)
                ind = find(vc == uc(c) & vt == ut(t),1); % C(A) won't make sense
                D1u(c,t,:) = D1(ind,:);
                eu = eu + 1;
            end
        end
        
        for j = 1 : n_perm     
            
            for c = 1:length(uc)
                ind = find(vc == uc(c));
                pereu = randperm(length(ut));
                for t = 1:length(ut)
                    ind2 = find(vt(ind) == ut(t));
                    D1p(ind(ind2),:) = ones(length(ind2),1)*squeeze(D1u(c,pereu(t),:))';
                end
            end
            
            B = pinv(D1p'*D1p)*D1p'*Xm;
            
            SSQ_time(1 + j) = (sum(sum((D1(:,parglmo.factors{2}.Dvars)*B(parglmo.factors{2}.Dvars,:)).^2))/dfB) / (sum(sum((D1(:,parglmo.interactions{1}.Dvars)*B(parglmo.interactions{1}.Dvars,:)).^2))/dfAB);
            
        end        % permutations
        
        
        p_value = (length(find(SSQ_time(2:n_perm + 1) >= SSQ_time(1))) + 1)/(n_perm);
        
        disp('p-value factor time')
        disp(p_value)
        
        if (p_value<0.05), e = 1; end

        
