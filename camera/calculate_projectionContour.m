function [contour, contour_index] = calculate_projectionContour(vertices, radius)
%calculate projection contour using projection vertices




shp = alphaShape(vertices(:,1),vertices(:,2), radius); % the value of obj.alphaShape 
[~,P] = boundaryFacets(shp);
boundary_index = zeros(size(P,1),1);
for ll=1:size(P,1)
    id_tem = find(vertices(:,1)==P(ll,1) & vertices(:,2)==P(ll,2));
    
    if (size(id_tem,1)==1)
        boundary_index(ll) = id_tem;   % only select the points that is not overlap
    else
        boundary_index(ll) = NaN;      % use a simple way to remove the redundant boundary
    end
end
boundary_index(isnan(boundary_index)) = [];  % remove nan

contour_index = boundary_index;
contour = vertices(contour_index,1:2);        % contour of model

end