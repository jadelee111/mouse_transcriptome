function imageCoord=aibs_portal_referenceToImage(reference_space_id,x,y,z,section_data_set_id)
    imageCoord=[];
    apiPath = 'http://api.brain-map.org/api/v2/';
    urlQuerry=[apiPath 'reference_to_image/' num2str(reference_space_id) '.xml?x=' ...
        num2str(x) '&y=' num2str(y) '&z=' num2str(z) '&section_data_set_ids=' num2str(section_data_set_id)];
    str=urlread(urlQuerry);
    tmp=xml_get_property(str,'section-data-set-id');
    if isempty(tmp)
        return;
    end
    dataId=str2double(tmp);
    if dataId~=section_data_set_id
        return;
    end
    tmp=xml_get_property(str,'section-image-id');
    if isempty(tmp)
        return;
    end
    imageId=str2double(tmp);
    tmp=xml_get_property(str,'section-number');
    if isempty(tmp)
        return;
    end
    sectionNum=str2double(tmp);
    tmp=xml_get_property(str,'x');
    if isempty(tmp)
        return;
    end
    x=str2double(tmp);
    tmp=xml_get_property(str,'y');
    if isempty(tmp)
        return;
    end
    y=str2double(tmp);
    
    imageCoord.dataId=dataId;
    imageCoord.imageId=imageId;
    imageCoord.sectionNum=sectionNum;
    imageCoord.x=x;
    imageCoord.y=y;
end