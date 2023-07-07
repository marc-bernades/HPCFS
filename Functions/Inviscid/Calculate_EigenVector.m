function K = Calculate_EigenVector(Flux,dFlux,sos,P)

u = dFlux(2,:)./dFlux(1,:);
v = dFlux(3,:)./dFlux(1,:);
w = dFlux(4,:)./dFlux(1,:);
E = dFlux(5,:)./dFlux(1,:);

for ii = 1:length(Flux(1,:))

    % First row
    K{ii}(1,:) = [1 0 0 0 1];

    % Second row
    K{ii}(2,1) = u(ii) - sos(ii);
    K{ii}(2,2) = u(ii);
    K{ii}(2,3) = 0;
    K{ii}(2,4) = 0;
    K{ii}(2,5) = u(ii) + sos(ii);

    % Third row
    K{ii}(2,1) = v(ii);
    K{ii}(2,2) = v(ii);
    K{ii}(2,3) = 1;
    K{ii}(2,4) = 0;
    K{ii}(2,5) = dFlux(3,ii)./dFlux(1,ii);

    % Fourth row
    K{ii}(2,1) = dFlux(4,ii)./dFlux(1,ii);
    K{ii}(2,2) = dFlux(4,ii)./dFlux(1,ii);
    K{ii}(2,3) = 0;
    K{ii}(2,4) = 1;
    K{ii}(2,5) = dFlux(4,ii)./dFlux(1,ii);

    % Fifth row
    K{ii}(2,1) = dFlux(5,ii)./dFlux(1,ii) + P(ii)/dFlux(1,ii) - dFlux(5,ii)./dFlux(1,ii)*sos(ii);
    K{ii}(2,2) = dFlux(5,ii)./dFlux(1,ii) + P(ii)/dFlux(1,ii) - dFlux(5,ii)./dFlux(1,ii)*sos(ii);
    K{ii}(2,3) = 0;
    K{ii}(2,4) = 0;
    K{ii}(2,5) = dFlux(2,ii)./dFlux(1,ii) + sos(ii);

end


end