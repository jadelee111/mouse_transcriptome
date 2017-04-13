function out_str=json_get_property(in_str, property_str)
    id_a=strfind(in_str,['"' property_str '":']);
    if isempty(id_a)
        out_str='';
        return;
    else
        str=in_str(id_a(1):end);
        id_b=strfind(str,',');
        if isempty(id_b)
            id_b=strfind(str,'}');
            if isempty(id_b)
                out_str='';
                return;
            end
        end
    end
    out_str=str((length(property_str)+4):(id_b(1)-1));
end