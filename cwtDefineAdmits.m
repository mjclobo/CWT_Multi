%% Define admittance output structures

n = 1;
admittances.wl = dataIn;

for k=ref.constits.names
    admittances.constits.(k).amps = constits.(k).amps ./ ref.constits.(k).amps;
    admittances.constits.(k).abs_phase = constits.(k).phases - ref.constits.(k).phases;
    % admittances.constits.(k).abs_phase = ref.constits.(k).phases - constits.(k).phases;
    n=n+1;
end

n = 1;
for k=ref.species.names
    admittances.species.(k).amps = species.(k).amps ./ ref.species.(k).amps;
    admittances.species.(k).abs_phase = species.(k).phases - ref.species.(k).phases;
    n=n+1;
end

