function mustBeEqualSize(a,b)
    % Test for equal size
    if ~isequal(size(a),size(b))
        eid = 'Size:notEqual';
        msg = 'Size of density_map, sos_map and source must be the same.';
        throwAsCaller(MException(eid,msg))
    end
end
