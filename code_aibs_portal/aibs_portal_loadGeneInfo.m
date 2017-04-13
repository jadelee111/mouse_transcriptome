function geneInfo=aibs_portal_loadGeneInfo(section_image_data_id)
	geneInfo=[];
    apiPath = 'http://api.brain-map.org/api/v2/';    
    urlQuerry=[apiPath 'data/SectionDataSet/query.xml?id=' num2str(section_image_data_id) '&include=genes'];
    str=urlread(urlQuerry);
    
    acronym=xml_get_property(str,'acronym');
    name=xml_get_property(str,'name');
    entrez_id=xml_get_property_entrez_id(str,'entrez-id');
    
    if isempty(acronym) || isempty(name)
        return;
    end
    
    geneInfo.dataId=section_image_data_id;
    geneInfo.acronym=acronym;
    geneInfo.name=name;
    geneInfo.entrez_id=entrez_id;
end