function out_str=xml_get_property(in_str, property_str)
    id_a=strfind(in_str,['<' property_str '>']);
    id_b=strfind(in_str,['</' property_str '>']);
    if isempty(id_a) || isempty(id_b)
        out_str='';
    else
        out_str=in_str((id_a(1)+length(property_str)+2):(id_b(1)-1));
    end
end
