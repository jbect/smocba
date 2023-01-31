function alpha = sMOCBA_iteration (obj_mean, obj_var)

% Equation numbers in the comments refer to Chapter 6 of Chen & Lee (2010).

m = size (obj_mean, 1);  % Input set size (number of designs)
q = size (obj_mean, 2);  % Output dimension  (number of objectives)
assert (isequal (size (obj_mean), [m, q]));
assert (isequal (size (obj_var),  [m, q]));

% `kk` will contain the indices computed using (6.41)
kk = nan (m, m);

% `trucmax` will contain the maximal value for the argmax in (6.41)
trucmax = nan (m, m);

for j = 1:m
    for i = 1:m

        % Temp vector to store what's inside the argmax of (6.41)
        truc = nan (1, q);

        % This `k` is the index `l` for the argmax of (6.41)
        for k = 1:q

            % `delta_{ijk}` is defined by (6.39)
            delta = obj_mean(j,k) - obj_mean(i,k);

            % The thing inside the argmax of (6.41)
            truc(k) = delta * abs (delta) / (obj_var(i,k) + obj_var(j,k));
        end

        % Compute the argmax in (6.41), and save the max too
        [trucmax(i,j), kk(i,j)] = max (truc);
    end
end

% `jj` contains the (argmin) indices defined by (6.42)
% `trucminmax` contains the associated minimal values
[trucminmax, jj] = min (trucmax' + diag (inf * ones(1, m)));

% `S_A` will contain the indicator vector for the set S_A in (6.43)
S_A = true (1, m);

% Construct `S_A`
for h = 1:m

    % Construct the indicator of the set Theta_h in (6.45)
    Theta = (jj == h);

    % Decide if h is in S_A, according to (6.43)
    if any (Theta)
        tmp_LHS = abs (trucminmax(h));
        tmp_RHS = min (abs (trucmax(Theta, h)));
        S_A(h) = tmp_LHS < tmp_RHS;
    end

    % Implementation note: when Theta_h is empty, we set S_A(h) to true.
    % This special case does not appear to be discussed by Chen & Lee,
    % but our choice is consistent with the fact that the minimum of the
    % empty set is usually taken to be +Inf.
end

% `alpha` will contain the proportion computed from (6.38) or (6.39)
alpha = nan(1, m);

% Construct the part of `alpha` that corresponds to `S_A`,
% omitting the normalization which is not needed at this stage
for h = 1:m
    if S_A(h)
        j = jj(h);
        k = kk(h, j);
        alpha(h) = obj_var(h,k) / ((obj_mean(h,k) - obj_mean(j,k))^2);
    end
end

% Construct the part of `alpha` that corresponds to `S_B`
for d = 1:m
    if ~ S_A(d)

        % Construct the indicator of the set Theta^*_h defined by (6.45)
        ThetaStar = S_A & (jj == d);

        % Compute the sum in (6.39)
        S = 0;
        for h = 1:m
            if ThetaStar(h)
                k = kk(h, d);
                S = S + obj_var(d,k) / obj_var(h,k) * alpha(h)^2;
            end
        end

        % And don't forget the square root
        alpha(d) = sqrt (S);
    end
end

% Finally, normalize to get the proportions
alpha = alpha / (sum (alpha));

end % function
