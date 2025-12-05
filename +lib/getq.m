% Returns the vector q based on the number of bodies and x, y and phi
% symbolic expressions
function q = getq(x, y, phi, nb)
    % Define
    q = sym(zeros(nb*3, 1));
    itterations = 0;

    % Create q
    for i = 1:3:nb*3
        itterations = itterations + 1;
        q(i) = x(itterations);
        q(i+1) = y(itterations);
        q(i+2) = phi(itterations);
    end

end