
function atlasCoord=aibs_portal_imageToAtlas(section_image_id, x, y, atlas_id)
    atlasCoord=[];
    apiPath = 'http://api.brain-map.org/api/v2/';
    urlQuerry=[apiPath 'image_to_atlas/' num2str(section_image_id) '.xml?x=' ...
        num2str(x) '&y=' num2str(y) '&atlas_id=' num2str(atlas_id)];
    str=urlread(urlQuerry);
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
    
    atlasCoord.atlasId=atlas_id;
    atlasCoord.imageId=imageId;
    atlasCoord.sectionNum=sectionNum;
    atlasCoord.x=x;
    atlasCoord.y=y;
end