function projected_power_objective = objective_power(spec_Hz, spec_Ey, Hz, Ey)
    project_Ey = spec_Ey * (spec_Ey' * Ey);
    project_Hz = spec_Hz * (spec_Hz' * Hz);
    
    % Ey'*Hz + conj(Ey'*Hz) where Ey and Hz are projected into the mode Ey
    % and Hz so that we get how much power is flowing through the end of
    % the waveguide in the proper mode
    
    projected_power_objective = -real(project_Ey' * project_Hz);    
end

