function [mat] = setOrientation(quat)
    % setOrientation converts a scalar-first quaternion to a rotation matrix.
    % Input: quat = [w, x, y, z]
    % Output: 3x3 Rotation Matrix

    % 1. Normalize the quaternion to ensure a valid rotation matrix
    n = norm(quat);
    if n == 0
        % Handle zero-vector case (return identity or error)
        mat = eye(3); 
        return;
    end
    quat = quat / n;

    % 2. Extract components
    w = quat(1);
    x = quat(2);
    y = quat(3);
    z = quat(4);

    % 3. Calculate Matrix
    % Uses the user-provided sign convention (Inverted/Frame Transform)
    mat = [1-2*(y^2+z^2),  2*(x*y+z*w),    2*(x*z-y*w);
           2*(x*y-z*w),    1-2*(x^2+z^2),  2*(y*z+x*w);
           2*(x*z+y*w),    2*(y*z-x*w),    1-2*(x^2+y^2)];
end