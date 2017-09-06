function [ particleOutput ] = dynUpdateParticles( particle, fDYN, timeVector)
% Non-Linear MonteCarlo propagation

% particle (nxm), n state vector size, m number of particles
% fDYN matlab function handle @(t,x)
% timeVector = [init time, final time]
    sizeP = size(particle);

    for i=1:sizeP(2)
        particleOutput(:,i) = deval(ode45( @(t,x) tempCall(fDYN,t,x),timeVector, particle(:,i)),timeVector(2));
    end


end

