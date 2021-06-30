function [points_new] = update_pts_position(obj, points, num_nearestpts)
%update_model_vertex: update model vertices accouding to new EDNode rotation and translation






% deside which to be updated
% num_nearestpts = obj.num_nearestpts_skel;
% points = obj.modelVertices;                   % start from the initial instead of updated

[pts_weights_id, pts_weights ] = updateWeight_knn(points', obj.node_position', num_nearestpts); 
num_pts = size(points,2);
points_new = zeros(3,num_pts);

for i=1:num_pts
    pts_i  = points(:,i);
    weight_i  = pts_weights(i,:);
    weight_id = pts_weights_id(i,:);
    pts_tem = zeros(3,1);
    for j=1:num_nearestpts
        
        node_id = weight_id(j);
        
        node_R = obj.node_rotation(3*node_id-2:3*node_id,:);
        node_p = obj.node_position(:,node_id);
        node_t = obj.node_translation(:,node_id);
        mapped_pts = node_R * (pts_i - node_p) + node_p + node_t;   %%%Yanhao20181225: Eq. (1)
        pts_tem = pts_tem + weight_i(j) * mapped_pts;
    end
    points_new(:,i) = pts_tem;
end
end