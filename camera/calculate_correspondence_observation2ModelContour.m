function [indexCorrespondence_observation2ModelContour,observation_removedNotCorrespondendingPoint] = ...
    calculate_correspondence_observation2ModelContour(observation,observation_normal,modelContour,modelContour_normal,theta_threshold,dist_threshold)
%Calculate correspondence according to normal contour
% input: observation M*2
%        observation_normal
%        model_contour N*2 (N>M)
%        model_normal
%        theta_threshold: threshold for normal selection
%        dist_threshold: distance threshold, only the correspondence smaller than this value are output
% input data has been processed in the main scrips: calculate_downsampleObservation_normal_correspondence

%% deal with theta_threshold
% do not need to normalize the nv, theta only has one value from 0 to pi, therefore we do not need acos
% dot product should be larger than cos(theta_threshold)
theta_threshold_normalProduct = cos(theta_threshold);

num_observation = size(observation,1);
num_modelContour = size(modelContour,1);

indexCorrespondence_observation2ModelContour = zeros(num_observation,1);  % the index of model contour associated by downsampled observation. Normal is used


% put distance first
% Yanhao20190328: faster 5 times than before
% tic
for i=1:num_observation
    position_observationI = observation(i,:);
    normal_observationI = observation_normal(i,:);
    
    % distance within threshold
    tem_dist = repmat(position_observationI,num_modelContour,1)-modelContour;                         % same as repmat(A(a,:),m,1)-B;
    dist_observation2ModelContour = sqrt((tem_dist(:,1)).^2 + (tem_dist(:,2)).^2);    % only use the euclidean distance
    contourindex_distwithinThreshold = find(dist_observation2ModelContour < dist_threshold); % the index of contour point whose index is within threshold
    
    if ~isempty(contourindex_distwithinThreshold)
        
        % angle within threshold
        contournormal_distwithinThreshold = modelContour_normal(contourindex_distwithinThreshold,:);
        innerProduct_contournormal_distwithinThreshold = sum(normal_observationI.*contournormal_distwithinThreshold,2);   % do not need to normalize the nv
        
        tem_contourindex_thetawithinThreshold = find(innerProduct_contournormal_distwithinThreshold>theta_threshold_normalProduct); % tem index
        
        % nearest distance
        if ~isempty(tem_contourindex_thetawithinThreshold)
            contourindex_withinAllThreshold = contourindex_distwithinThreshold(tem_contourindex_thetawithinThreshold);
            contourpoint_withinAllThreshold = modelContour(contourindex_withinAllThreshold,:);
            
            if size(contourpoint_withinAllThreshold,1) ==1                   % only one point within all threshold
                association_index_contour = contourindex_withinAllThreshold;
            elseif size(contourpoint_withinAllThreshold,1) ==2                   % two point within all threshold
                tem_dist = dist_observation2ModelContour(contourindex_withinAllThreshold);  % distance between observation to these two points
                [~,tem_id] = min(tem_dist);
                association_index_contour = contourindex_withinAllThreshold(tem_id);    
            else
                tem_DT = delaunayTriangulation(contourpoint_withinAllThreshold);     % trangularize model boundary
                [association_temindex, ~] = nearestNeighbor(tem_DT,position_observationI(1),position_observationI(2));   % notice 这只是contour_point_tobechecked的index
                association_index_contour = contourindex_withinAllThreshold(association_temindex);                       % put it back to model contour index
            end
        else
            association_index_contour = -10000;
        end        
    else
        association_index_contour = -10000;
    end
    
    indexCorrespondence_observation2ModelContour(i) = association_index_contour;                               % write like this for readable purpose
end
% toc

%% remove the not corrected points
index_toBeRemovedFromObservation = find(indexCorrespondence_observation2ModelContour(:)==-10000);   % remove the un correct correspondence

indexCorrespondence_observation2ModelContour(index_toBeRemovedFromObservation) = [];
observation_removedNotCorrespondendingPoint = observation;
observation_removedNotCorrespondendingPoint(index_toBeRemovedFromObservation,:) = [];


% indexCorrespondence_observation2ModelContour(indexCorrespondence_observation2ModelContour(:)==-10000) = [];     % remove the un correct correspondence

if(0)
%% check result
observation_normal_removed = observation_normal;
observation_normal_removed(index_toBeRemovedFromObservation,:) = [];

figure
quiver(observation_removedNotCorrespondendingPoint(:,1),observation_removedNotCorrespondendingPoint(:,2),observation_normal_removed(:,1),observation_normal_removed(:,2),'g'); hold on
plot(observation(:,1),observation(:,2),'r.'); hold on;
plot(observation_removedNotCorrespondendingPoint(:,1),observation_removedNotCorrespondendingPoint(:,2),'r*'); hold on;
quiver(modelContour(:,1),modelContour(:,2),modelContour_normal(:,1),modelContour_normal(:,2),'y'); hold on
plot(modelContour(:,1),modelContour(:,2),'b.'); hold on;

for i=1:size(indexCorrespondence_observation2ModelContour,1)
    index_observation = i;
    index_model_contour = indexCorrespondence_observation2ModelContour(i);
    
    tem_observation = observation_removedNotCorrespondendingPoint(index_observation,:);
    tem_model = modelContour(index_model_contour,:);
    
    plot([tem_observation(1) tem_model(1)],[tem_observation(2) tem_model(2)],'k'); hold on
end
title(['Threshold: theta ', num2str(rad2deg(theta_threshold)), ', dist ', num2str(dist_threshold)]); hold on
axis equal

    
end

end