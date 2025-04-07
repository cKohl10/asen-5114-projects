function [Lg, C] = loopgain(G,s, params, K)

    for i = 1:length(params)
        C_set(i) = comp_select(params{i}, s);
    end

    C = K;
    for i = 1:length(C_set)
        C = C*C_set(i);
    end

    Lg = C*G;

end