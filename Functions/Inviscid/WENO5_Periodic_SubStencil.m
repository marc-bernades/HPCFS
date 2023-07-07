function index_Stencil = WENO5_Periodic_SubStencil(dimension,position,nStencil)



% Define stencil
index_Stencil = [position-2,position-1,position,position+1,position+2] + nStencil - 1; % Stencil 2 shifted right

% Periodic BC

if position == 2
    if nStencil == 1
        index_Stencil = [dimension-2,dimension-1,position,position+1,position+2];
    end

    if nStencil == 2
        index_Stencil = [dimension-1,position,position+1,position+2,position+3];
    end

elseif position == 3
    if nStencil == 1
        index_Stencil = [dimension-1,position-1,position,position+1,position+2];
    end

elseif position == dimension-3
    if nStencil == 2
        index_Stencil = [position-1,position,position+1,position+2,2];
    end

elseif position == dimension-2
    if nStencil == 1
        index_Stencil = [position-2,position-1,position,position+1,2];
    end
    if nStencil == 2
        index_Stencil = [position-1,position,position+1,2,3];
    end
elseif position == dimension-1
    if nStencil == 1
        index_Stencil = [position-2,position-1,position,2,3];
    else
        index_Stencil = [position-1,position,2,3,4];
    end
end



end