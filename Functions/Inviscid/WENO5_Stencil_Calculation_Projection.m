function [Fp,Fn] = WENO5_Stencil_Calculation_Projection(i,j,k,rho,u,v,w,E,P,sos, dimension)
% Initialize
index_Stencil = cell([1,2]);
Flux          = cell([1,2]);
dFlux         = cell([1,2]);
alpha         = cell([1,2]);


% Positive Stencil (j-2, j-1, j ,j+1, j+2) > Left
% Negative Stencil (j-1, j ,j+1, j+2, j+3) > Right
for ii = 1:2

    % Check position X, Y, Z

    %% X-Direction
    if dimension == 1

        % Define function evaluate stencil(dimension,current position,
        % substencil)
        index_Stencil{ii} = WENO5_Periodic_SubStencil(length(u(1,:,1)),j,ii);
  
        % Flux
        Flux{ii} = [rho(i,index_Stencil{ii},k).*u(i,index_Stencil{ii},k);
            rho(i,index_Stencil{ii},k).*u(i,index_Stencil{ii},k).*u(i,index_Stencil{ii},k) + P(i,index_Stencil{ii},k);
            rho(i,index_Stencil{ii},k).*u(i,index_Stencil{ii},k).*v(i,index_Stencil{ii},k);
            rho(i,index_Stencil{ii},k).*u(i,index_Stencil{ii},k).*w(i,index_Stencil{ii},k);
            rho(i,index_Stencil{ii},k).*u(i,index_Stencil{ii},k).*E(i,index_Stencil{ii},k) + P(i,index_Stencil{ii},k).*(u(i,index_Stencil{ii},k))];

        % Conserved variables
        dFlux{ii} = [rho(i,index_Stencil{ii},k);
            rho(i,index_Stencil{ii},k).*u(i,index_Stencil{ii},k);
            rho(i,index_Stencil{ii},k).*v(i,index_Stencil{ii},k);
            rho(i,index_Stencil{ii},k).*w(i,index_Stencil{ii},k);
            rho(i,index_Stencil{ii},k).*E(i,index_Stencil{ii},k)];

        % Project Flux and dFlux to characteristic space
        K = Calculate_EigenVector(Flux{ii},dFlux{ii},sos(i,index_Stencil{ii},k));
        [Flux,dFlux] = Project_CharacteristicSpace(Flux,dFlux,K);

        % Eigen values on each stencil
        alpha{ii} = max([max(abs(u(i,index_Stencil{ii},k))),max(abs(u(i,index_Stencil{ii},k) + sos(i,index_Stencil{ii},k))),max(abs(u(i,index_Stencil{ii},k) - sos(i,index_Stencil{ii},k)))]);
%         alpha{ii} = [max(abs(u(i,index_Stencil{ii},k) + sos(i,index_Stencil{ii},k)));
%             max(abs(u(i,index_Stencil{ii},k)));
%             max(abs(u(i,index_Stencil{ii},k)));
%              max(abs(u(i,index_Stencil{ii},k)));
%             max(abs(u(i,index_Stencil{ii},k) - sos(i,index_Stencil{ii},k)))];


    %% Y - Direction
    elseif dimension == 2

        % Define function evaluate stencil(dimension,current position,
        % substencil)
        index_Stencil{ii} = WENO5_Periodic_SubStencil(length(u(:,1,1)),i,ii);

        % Flux
        Flux{ii} = [rho(index_Stencil{ii},j,k)'.*v(index_Stencil{ii},j,k)';
            rho(index_Stencil{ii},j,k)'.*v(index_Stencil{ii},j,k)'.*u(index_Stencil{ii},j,k)';
            rho(index_Stencil{ii},j,k)'.*v(index_Stencil{ii},j,k)'.*v(index_Stencil{ii},j,k)' + P(index_Stencil{ii},j,k)';
            rho(index_Stencil{ii},j,k)'.*v(index_Stencil{ii},j,k)'.*w(index_Stencil{ii},j,k)';
            rho(index_Stencil{ii},j,k)'.*v(index_Stencil{ii},j,k)'.*E(index_Stencil{ii},j,k)' + P(index_Stencil{ii},j,k)'.*(v(index_Stencil{ii},j,k)')];

        % Conserved variables
        dFlux{ii} = [rho(index_Stencil{ii},j,k)';
            rho(index_Stencil{ii},j,k)'.*u(index_Stencil{ii},j,k)';
            rho(index_Stencil{ii},j,k)'.*v(index_Stencil{ii},j,k)';
            rho(index_Stencil{ii},j,k)'.*w(index_Stencil{ii},j,k)';
            rho(index_Stencil{ii},j,k)'.*E(index_Stencil{ii},j,k)'];

        % Eigen values on each stencil
        alpha{ii} = max([max(abs(v(index_Stencil{ii},j,k))),max(abs(v(index_Stencil{ii},j,k) + sos(index_Stencil{ii},j,k))),max(abs(v(index_Stencil{ii},j,k) - sos(index_Stencil{ii},j,k)))]);


    %% Z - Direction
    elseif dimension == 3

       % Define function evaluate stencil(dimension,current position,
        % substencil)
        index_Stencil{ii} = WENO5_Periodic_SubStencil(length(u(1,1,:)),k,ii);

        % Flux
        Flux{ii} = [rho(i,j,index_Stencil{ii}).*w(i,j,index_Stencil{ii});
            rho(i,j,index_Stencil{ii}).*w(i,j,index_Stencil{ii}).*u(i,j,index_Stencil{ii});
            rho(i,j,index_Stencil{ii}).*w(i,j,index_Stencil{ii}).*v(i,j,index_Stencil{ii});
            rho(i,j,index_Stencil{ii}).*w(i,j,index_Stencil{ii}).*w(i,j,index_Stencil{ii}) + P(i,j,index_Stencil{ii});
            rho(i,j,index_Stencil{ii}).*w(i,j,index_Stencil{ii}).*E(i,j,index_Stencil{ii}) + P(i,j,index_Stencil{ii}).*(w(i,j,index_Stencil{ii}))];

        Flux{ii} = reshape(Flux{ii},[5 5]);

        % Conserved variables
        dFlux{ii} = [rho(i,j,index_Stencil{ii});
            rho(i,j,index_Stencil{ii}).*u(i,j,index_Stencil{ii});
            rho(i,j,index_Stencil{ii}).*v(i,j,index_Stencil{ii});
            rho(i,j,index_Stencil{ii}).*w(i,j,index_Stencil{ii});
            rho(i,j,index_Stencil{ii}).*E(i,j,index_Stencil{ii})];

        dFlux{ii} = reshape(dFlux{ii},[5 5]);

        % Eigen values on each stencil
        alpha{ii} = max([max(abs(w(i,j,index_Stencil{ii}))),max(abs(w(i,j,index_Stencil{ii}) + sos(i,j,index_Stencil{ii}))),max(abs(w(i,j,index_Stencil{ii}) - sos(i,j,index_Stencil{ii})))]);    

    end



    

end

% Max eigen value
alpha_T = max([alpha{1},alpha{2}]);

% Local Fp , Fn = Lax-Friedrichs
Fp = 0.5*(Flux{1} + alpha_T.*dFlux{1});  %Right
Fn = 0.5*(Flux{2} - alpha_T.*dFlux{2});  %Left


end