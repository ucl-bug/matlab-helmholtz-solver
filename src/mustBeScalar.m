function mustBeScalar(x)
    % Test for equal size
    if ~isscalar(x)
        eid = 'NotScalar';
        msg = 'Value must be scalar';
        throwAsCaller(MException(eid,msg))
    end
end