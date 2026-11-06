function [theta_x, theta_y] = pixels2eyePosition(eye_to_monitor_distance, pixel_size, monitor_resolution)
% PIXELS2EYEPOSITION  Visual-field position (eccentricity) of every pixel on a flat screen.
%
%   [theta_x, theta_y] = pixels2eyePosition(eye_to_monitor_distance, pixel_size, monitor_resolution)
%
%   Inputs (same as pixels2eyeDegrees):
%     eye_to_monitor_distance  distance from the eye to the screen centre, in cm (e.g. 21.5)
%     pixel_size               side length of one pixel, in cm (e.g. 33/1080)
%     monitor_resolution       [width height] in pixels (e.g. [1920 1080])
%
%   Outputs: [height x width] maps, in degrees, signed and 0 at the screen centre
%   (theta_x negative to the left, theta_y negative towards the top row):
%     theta_x = atan(dx / D),  theta_y = atan(dy / D)
%   where dx, dy are the offsets (cm) of each PIXEL CENTRE from the screen centre and D is
%   the eye-to-screen distance. The eye is assumed to face the screen centre perpendicularly.
%
%   WHICH FUNCTION TO USE
%     pixels2eyePosition  WHERE a point is: axis ticks of receptive-field maps, positions of
%                         trajectories/cells. At the screen edge (16.5 cm at 21.5 cm) it gives
%                         atan(16.5/21.5) = 37.5 deg.
%     pixels2eyeDegrees   HOW BIG an object centred on the line of sight is: 2*atan(d/(2D))
%                         is the angle subtended by an object of size d at the screen centre
%                         (ball diameters, bar widths, spacing). Used as a position it
%                         overestimates the periphery (42.0 deg instead of 37.5 deg at the edge).
%   Both functions assume a flat screen; off-centre sizes subtend less than at the centre.

screen_center = monitor_resolution / 2;   % Centre of the screen in pixel units ([x y], measured from the outer edge)

[x, y] = meshgrid(1:monitor_resolution(1), 1:monitor_resolution(2));   % Pixel indices (1-based) for every pixel

delta_x = (x - 0.5 - screen_center(1)) * pixel_size;   % Horizontal offset (cm) of each pixel centre from the screen centre
delta_y = (y - 0.5 - screen_center(2)) * pixel_size;   % Vertical offset (cm) of each pixel centre from the screen centre

theta_x = atand(delta_x / eye_to_monitor_distance);    % Horizontal eccentricity in degrees: atan(dx / D)
theta_y = atand(delta_y / eye_to_monitor_distance);    % Vertical eccentricity in degrees: atan(dy / D)

end
