function popchange = NEWu_opioid_postvuln_social(t,pops,p)

S1 = pops(p.k+1);
S2 = pops(p.k+2);
A  = pops(p.k+3);
R  = pops(p.k+4);

popchange = zeros(p.num_classes,1);

social_rate = p.beta * A / p.N;
fval = f_c_eval(p.c,p);

%% Prescription classes
popchange(1) = p.rho*S1 ...
             - (p.alpha*p.k + p.sigma_vec(1) + social_rate + p.mu)*pops(1);

if p.cutoff > 1
    for i = 2:p.cutoff
        popchange(i) = p.alpha*p.k*pops(i-1) ...
                     - (p.alpha*p.k + p.sigma_vec(i) + social_rate + p.mu)*pops(i);
    end
end

if p.cutoff < p.k
    for i = p.cutoff+1:p.k
        popchange(i) = 0;
    end
end

%% Flows
prescription_exit = p.alpha * p.k * pops(p.cutoff);

presc_only_to_A = sum(pops(1:p.cutoff)'.*p.sigma_vec(1:p.cutoff));
social_from_P   = social_rate * sum(pops(1:p.cutoff));
social_from_S1  = social_rate * S1;
s2_to_A         = fval * S2;

first_to_A      = presc_only_to_A + social_from_P + social_from_S1 + s2_to_A;
relapse_to_A    = p.delta * R;

%% S1
popchange(p.k+1) = -(p.rho + social_rate)*S1 ...
                 + p.gammaS2*S2 ...
                 - p.mu*S1 ...
                 + p.mu*(S1 + sum(pops(1:p.k)) + S2 + A + R) ...
                 + p.muOD*A;

%% S2
popchange(p.k+2) = prescription_exit ...
                 - s2_to_A ...
                 - p.gammaS2*S2 ...
                 - p.mu*S2;

%% A
popchange(p.k+3) = first_to_A ...
                 - p.eta*A ...
                 + relapse_to_A ...
                 - (p.mu + p.muOD)*A;

%% R
popchange(p.k+4) = p.eta*A - p.delta*R - p.mu*R;

%% Cumulative outputs
popchange(p.k+5)  = p.rho*S1;
popchange(p.k+6)  = first_to_A;
popchange(p.k+7)  = relapse_to_A;
popchange(p.k+8)  = presc_only_to_A;
popchange(p.k+9)  = social_from_S1 + social_from_P;
popchange(p.k+10) = s2_to_A;
popchange(p.k+11) = social_from_S1;

end

function val = f_c_eval(c,p)

    switch lower(p.f_form)

        case 'decreasing'
            g = p.cshift / (c + p.cshift);
            val = p.fbase + p.Adec * (g - p.gdec_bar);

        case 'increasing'
            g = c / (c + p.cshift);
            val = p.fbase + p.Ainc * (g - p.ginc_bar);

        case 'ushaped'
            c_days = 365*c;
            val = p.fbase + p.mU * ...
                ((c_days - p.u_mid_day)^2 - p.u_quad_bar);

        otherwise
            error('Unknown f_form.');
    end
end
