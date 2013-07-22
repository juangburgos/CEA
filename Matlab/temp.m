% % QR Decomposition 1 (50%)
% % % Initialize variables
% % n = size(R,1);
% % c = zeros(n,1);
% % d = zeros(n,1);
% % 
% % % Perform Householder QR decomposition
% % for k = 1:(n-1)
% %     scale = 0.0;
% %     for i = k:n
% %         scale = max(scale, abs(R(i,k)));
% %     end
% %     if (scale == 0.0)
% %         c(k) = 0.0;
% %         d(k) = 0.0;
% %     else
% %         for i = k:n
% %             R(i,k) = R(i,k) / scale;
% %         end
% %         sum = 0.0;
% %         for i = k:n
% %             sum = sum + R(i,k) * R(i,k);
% %         end
% %         sigma = sqrt(sum) * sign(R(k,k));
% %         R(k,k) = R(k,k) + sigma;
% %         c(k) = sigma * R(k,k);
% %         d(k) = -1.0 * scale * sigma;
% %         for j = (k+1):n
% %             sum = 0.0;
% %             for i = k:n
% %                 sum = sum + R(i,k) * R(i,j);
% %             end
% %             tau = sum / c(k);
% %             for i = k:n
% %                 R(i,j) = R(i,j) - tau * R(i,k);
% %             end
% %         end
% %     end
% % end
% % d(n) = R(n,n);
% % 
% % % Construct Q and erase temporary data in R
% % % TODO: Could inegrate with loop in next step
% % % TODO: Rewrite in C compatibility ****************************************
% % Q = eye(n);
% % for j = (n-1):-1:1
% %     Q(j:n,:) = Q(j:n,:) - ( (1.0 / c(j)) * R(j:n,j) ) * ( R(j:n,j)' * Q(j:n,:) );
% %     R((j+1):n,j) = 0;
% %     R(j,j) = d(j);
% % end