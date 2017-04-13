function imageInfo=aibs_portal_loadAtlasImageInfo(section_image_id, atlasName)
	if nargin<2
        atlasName='Atlas';
    end
    imageInfo=[];
    apiPath = 'http://api.brain-map.org/api/v2/';    
    urlQuerry=[apiPath 'data/AtlasImage/' num2str(section_image_id) '.xml?include=alternate_images'];
    str=urlread(urlQuerry);
    
    tmp=xml_get_property(str,'id');
    if isempty(tmp)
        return;
    end
    imageId=str2double(tmp);
    if imageId~=section_image_id
        return;
    end
    tmp=xml_get_property(str,'path');
    if isempty(tmp)
        return;
    end
    path=tmp;
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
    stridx=strfind(str,['<image-type>' atlasName]);
    if isempty(stridx)
        return;
    end
    str=str(stridx(1):end);
    tmp=xml_get_property(str,'path');
    if isempty(tmp)
        return;
    end
    annoPath=tmp;
    
    imageInfo.imageId=imageId;
    imageInfo.path=path;
    imageInfo.annoPath=annoPath;
    imageInfo.x=x;
    imageInfo.y=y;
end