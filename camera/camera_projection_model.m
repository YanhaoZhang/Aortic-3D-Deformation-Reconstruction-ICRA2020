function [proj] = camera_projection_model(v, R, t, s)

U = [1 0 0;0 1 0];
proj = U*s*R*(v+t);


end