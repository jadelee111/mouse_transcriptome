function imageInfo=aibs_portal_loadImageInfo(section_data_set_id, section_image_id)
	imageInfo=[];
    apiPath = 'http://api.brain-map.org/api/v2/';
    urlQuerry=[apiPath 'data/SectionDataSet/' num2str(section_data_set_id) '.json?include=section_images,genes'];
    str=urlread(urlQuerry);
    image_idx=strfind(str,['"id":' num2str(section_image_id) ',']);
    if isempty(image_idx)
        return;
    end
    start_idx=strfind(str,'{');
    end_idx=strfind(str,'}');
    id1=find(start_idx<image_idx(1));
    id2=find(end_idx>image_idx(1));
    if isempty(id1) || isempty(id2) || start_idx(id1(end))>end_idx(id2(1))
        return;
    end
    str=str(start_idx(id1(end)):end_idx(id2(1)));

    tmp=json_get_property(str,'id');
    if isempty(tmp)
        return;
    end
    imageId=str2double(tmp);
    if imageId ~= section_image_id
        return;
    end
    tmp=json_get_property(str,'path');
    if isempty(tmp)
        return;
    end
    path=tmp(2:end-1);
    tmp=json_get_property(str,'x');
    if isempty(tmp)
        return;
    end
    x=str2double(tmp);
    tmp=json_get_property(str,'y');
    if isempty(tmp)
        return;
    end
    y=str2double(tmp);
    tmp=json_get_property(str,'image_width');
    if isempty(tmp)
        return;
    end
    w=str2double(tmp);
    tmp=json_get_property(str,'image_height');
    if isempty(tmp)
        return;
    end
    h=str2double(tmp);
    
    imageInfo.imageId=imageId;
    imageInfo.h=h;
    imageInfo.w=w;
    imageInfo.path=path;
    imageInfo.x=x;
    imageInfo.y=y;
end