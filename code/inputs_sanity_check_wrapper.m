function [in] = inputs_sanity_check_wrapper(in)
    [in,errors,warnings] = inputs_sanity_check(in);
    if ~isempty(errors) > 0
        error(errors{1})
    end
    if ~isempty(warnings) > 0
        warnings{1}
    end
end