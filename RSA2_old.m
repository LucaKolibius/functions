% + misses
% +other +hipp

% %% cortex
% % RSA for hits + respLocked
% cd ../cortex
% if enoughC==1
%     glmY=[];
%     glmX=zeros(size(hitsIdx,1),2);
%     for i=1:size(hitsIdx,1)*size(hitsIdx,1)
%         glmY(i) = RSAhits_respLocked_c(i);
%         if RSA_mask_hits(i)==1 % main diagnonal
%             glmX(i,1)=1;
%             glmX(i,2)=1;
%         elseif RSA_mask_hits(i)==2
%             glmX(i,2)=1;
%         end
%     end
%
%     % [RSAcoeff_hitsResp, RSAdev_hitsResp, RSAstats_hitsResp]=glmfit(glmX,glmY,'normal');
%     glm_hitsResp_c =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%     % GLM for hits + cueLocked
%     glmY=[];
%     glmX=zeros(size(hitsIdx,1),2);
%     for i=1:size(hitsIdx,1)*size(hitsIdx,1)
%         glmY(i)=RSAhits_cueLocked_c(i);
%         if RSA_mask_hits(i)==1 % main diagnonal
%             glmX(i,1)=1;
%             glmX(i,2)=1;
%         elseif RSA_mask_hits(i)==2
%             glmX(i,2)=1;
%         end
%     end
%
%     % [RSAcoeff_hitsCue, RSAdev_hitsCue, RSAstats_hitsCue] = glmfit(glmX,glmY,'normal');
%     glm_hitsCue_c =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%     % miss cuelocked cortex
%     if enoughMiss == 1
%         glmY=[];
%         glmX=zeros(size(missIdx,1),2);
%         for i=1:size(missIdx,1)*size(missIdx,1)
%             glmY(i)=RSAmiss_cueLocked_c(i);
%             if RSA_mask_miss(i)==1 % main diagnonal
%                 glmX(i,1)=1;
%                 glmX(i,2)=1;
%             elseif RSA_mask_miss(i)==2
%                 glmX(i,2)=1;
%             end
%         end
%
%         % [RSAcoeff_missCue, RSAdev_missCue, RSAstats_missCue] = glmfit(glmX,glmY,'normal');
%         glm_missCue_c =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%         % miss resplocked cortex
%         glmY=[];
%         glmX=zeros(size(missIdx,1),2);
%         for i=1:size(missIdx,1)*size(missIdx,1)
%             glmY(i)=RSAmiss_respLocked_c(i);
%             if RSA_mask_miss(i)==1 % main diagnonal
%                 glmX(i,1)=1;
%                 glmX(i,2)=1;
%             elseif RSA_mask_miss(i)==2
%                 glmX(i,2)=1;
%             end
%         end
%
%         % [RSAcoeff_missResp, RSAdev_missResp, RSAstats_missResp] = glmfit(glmX,glmY,'normal');
%         glm_missResp_c =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%         % save it to GLM_output
%         save(['glmfit_c_', subjID, '.mat'], 'glm_hitsCue_c', 'glm_hitsResp_c', 'glm_missCue_c', 'glm_missResp_c');
%
%     elseif enoughC== 1 && enoughMiss == 0 %
%         save(['glmfit_c_noMISS', subjID, '.mat'], 'glm_hitsCue_c', 'glm_hitsResp_c');
%     end
% end
%
% %% hippocampus
% cd ../hippocampus/
% if enoughH == 1
%     % GLM for hits + respLocked
%     glmY=[];
%     glmX=zeros(size(hitsIdx,1),2);
%     for i=1:size(hitsIdx,1)*size(hitsIdx,1)
%         glmY(i) = RSAhits_respLocked_h(i);
%         if RSA_mask_hits(i)==1 % main diagnonal
%             glmX(i,1)=1;
%             glmX(i,2)=1;
%         elseif RSA_mask_hits(i)==2
%             glmX(i,2)=1;
%         end
%     end
%
%     % [RSAcoeff_hitsResp, RSAdev_hitsResp, RSAstats_hitsResp]=glmfit(glmX,glmY,'normal');
%     glm_hitsResp_h =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%     % GLM for hits + cueLocked
%     glmY=[];
%     glmX=zeros(size(hitsIdx,1),2);
%     for i=1:size(hitsIdx,1)*size(hitsIdx,1)
%         glmY(i)=RSAhits_cueLocked_h(i);
%         if RSA_mask_hits(i)==1 % main diagnonal
%             glmX(i,1)=1;
%             glmX(i,2)=1;
%         elseif RSA_mask_hits(i)==2
%             glmX(i,2)=1;
%         end
%     end
%
%     % [RSAcoeff_hitsCue, RSAdev_hitsCue, RSAstats_hitsCue] = glmfit(glmX,glmY,'normal');
%     glm_hitsCue_h =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%     % Misses cue locked hipp
%     if enoughMiss == 1 % already within enoughH==1 if statement
%         glmY=[];
%         glmX=zeros(size(missIdx,1),2);
%         for i=1:size(missIdx,1)*size(missIdx,1)
%             glmY(i)=RSAmiss_cueLocked_h(i);
%             if RSA_mask_miss(i)==1 % main diagnonal
%                 glmX(i,1)=1;
%                 glmX(i,2)=1;
%             elseif RSA_mask_miss(i)==2
%                 glmX(i,2)=1;
%             end
%         end
%
%         % [RSAcoeff_missCue, RSAdev_missCue, RSAstats_missCue] = glmfit(glmX,glmY,'normal');
%         glm_missCue_h =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%         % Misses resp locked hipp
%         glmY=[];
%         glmX=zeros(size(missIdx,1),2);
%         for i=1:size(missIdx,1)*size(missIdx,1)
%             glmY(i)=RSAmiss_respLocked_h(i);
%             if RSA_mask_miss(i)==1 % main diagnonal
%                 glmX(i,1)=1;
%                 glmX(i,2)=1;
%             elseif RSA_mask_miss(i)==2
%                 glmX(i,2)=1;
%             end
%         end
%
%         % [RSAcoeff_missResp, RSAdev_missResp, RSAstats_missResp] = glmfit(glmX,glmY,'normal');
%         glm_missResp_h =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%         % save it to GLM_output
%         save(['glmfit_h_', subjID, '.mat'], 'glm_hitsCue_h', 'glm_hitsResp_h', 'glm_missCue_h', 'glm_missResp_h');
%
%     elseif enoughMiss == 0
%         save(['glmfit_h_noMiss', subjID, '.mat'], 'glm_hitsCue_h', 'glm_hitsResp_h');
%     end
% end
%
% %% other
% cd ../other/
% if enoughO == 1
%     % GLM for hits + respLocked
%     glmY=[];
%     glmX=zeros(size(hitsIdx,1),2);
%     for i=1:size(hitsIdx,1)*size(hitsIdx,1)
%         glmY(i) = RSAhits_respLocked_o(i);
%         if RSA_mask_hits(i)==1 % main diagnonal
%             glmX(i,1)=1;
%             glmX(i,2)=1;
%         elseif RSA_mask_hits(i)==2
%             glmX(i,2)=1;
%         end
%     end
%
%     % [RSAcoeff_hitsResp, RSAdev_hitsResp, RSAstats_hitsResp]=glmfit(glmX,glmY,'normal');
%     glm_hitsResp_o =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%     % GLM for hits + cueLocked
%     glmY=[];
%     glmX=zeros(size(hitsIdx,1),2);
%     for i=1:size(hitsIdx,1)*size(hitsIdx,1)
%         glmY(i)=RSAhits_cueLocked_o(i);
%         if RSA_mask_hits(i)==1 % main diagnonal
%             glmX(i,1)=1;
%             glmX(i,2)=1;
%         elseif RSA_mask_hits(i)==2
%             glmX(i,2)=1;
%         end
%     end
%
%     % [RSAcoeff_hitsCue, RSAdev_hitsCue, RSAstats_hitsCue] = glmfit(glmX,glmY,'normal');
%     glm_hitsCue_o =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%     % Misses cue locked other
%     if enoughMiss == 1 % already within enoughO==1 if statement
%         glmY=[];
%         glmX=zeros(size(missIdx,1),2);
%         for i=1:size(missIdx,1)*size(missIdx,1)
%             glmY(i)=RSAmiss_cueLocked_o(i);
%             if RSA_mask_miss(i)==1 % main diagnonal
%                 glmX(i,1)=1;
%                 glmX(i,2)=1;
%             elseif RSA_mask_miss(i)==2
%                 glmX(i,2)=1;
%             end
%         end
%
%         % [RSAcoeff_missCue, RSAdev_missCue, RSAstats_missCue] = glmfit(glmX,glmY,'normal');
%         glm_missCue_o =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%         % Misses resp locked other
%         glmY=[];
%         glmX=zeros(size(missIdx,1),2);
%         for i=1:size(missIdx,1)*size(missIdx,1)
%             glmY(i)=RSAmiss_respLocked_o(i);
%             if RSA_mask_miss(i)==1 % main diagnonal
%                 glmX(i,1)=1;
%                 glmX(i,2)=1;
%             elseif RSA_mask_miss(i)==2
%                 glmX(i,2)=1;
%             end
%         end
%
%         glm_missResp_o =  fitglm(glmX,glmY,'linear','Distribution','normal', 'Intercept', true)
%
%         % save it to GLM_output
%         save(['glmfit_o_', subjID, '.mat'], 'glm_hitsCue_o', 'glm_hitsResp_o', 'glm_missCue_o', 'glm_missResp_o');
%
%     elseif enoughMiss == 0
%         save(['glmfit_o_noMiss', subjID, '.mat'], 'glm_hitsCue_o', 'glm_hitsResp_o');
%     end
% end