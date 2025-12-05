% Returns orientation matrix

function A = orientationMatrix2D(angle)
    A = [cos(angle), -sin(angle);
        sin(angle) cos(angle)];
end