function [Jacobian] = camera_projection_model_Jacobian(dv, R, t, s)

U = [1 0 0;0 1 0];
% proj = U*s*R*(v+t);

Jacobian = U*s*R*dv;

end