function result = tensor_ops(M, operation, varargin)

% Turns tensor operations into sparse matrix-vector products (and viceversa)
%
% author: Antonio Stanziola
% date: 11th April 2020
% last update: 11th Feb 2021 (Antonio Stanziola)

    switch operation
        case 'matrix_to_vector'
            result = M(:);

        case 'matrix_to_elementwise'
            result = sparse(diag(M(:)));

        case 'matrix_to_LHS_product'
            result = sparse(kron(eye(size(M)),M)); 

        case 'matrix_to_RHS_product'
            result = sparse(kron(transpose(M),eye(size(M))));

        case 'vector_to_matrix'
            matrix_shape = varargin{1};
            result = reshape(M, matrix_shape(1),matrix_shape(2));

        otherwise
            error('Unknown operation')
    end
end
    