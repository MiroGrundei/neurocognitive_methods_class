function D = SBL_rm_scale(D, val)

% INPUT
% D = meeg object with at least one forward model
% val = index of D.inv field - inversion to be descaled.
% OUTPUT
% D = meeg object with scaling removed for D.inv{val} specifically

if D.inv{val}.forward.sensors.tra(1,1) > 1e15 % check if scaling is present
    temp = D.inv{val}.forward.sensors.tra;
    temp = temp * 1e-18;
    D.inv{val}.forward.sensors.tra = temp;
    if isfield(D.inv{val}.forward, 'scale')
        D.inv{val}.forward = rmfield(D.inv{val}.forward, 'scale');
    end
    disp('Scaling found and removed')
else
    disp('No scaling seems to be present. No modifications made.')
end
