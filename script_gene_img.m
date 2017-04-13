%% Code for second half of webpage:annoation, regionally expressed genes 
% Calculate the slice with the highest expressions and query the AMBA
% database to obtain related images. Then find the relatively high and
% low expressed genes, genes that use the dictionary and btain the images 
% by querying the AMBA

basename='coronal';
load ./data/data_set_id.mat; %gene ids
section_data_id=foldername;
sizeGrid = [67 41 58];
res=200;
referenceSpaceId=10;
atlasId=1;
atlasName='Atlas - Adult Mouse';
dictSize=[100,200,400,600,800,1000];

%load mask
load coronal_master_mask.mat

%path of aibs api portal
addpath('../code_aibs_portal');

%load gene expression matrix
VAR='masked_data_sel_gene_coronal';
gmat=load('masked_data_sel_gene_coronal_final.mat',VAR);
gmat=gmat.(VAR);
gmat_mean=mean(gmat,2);
gmat_std=std(gmat,0,2);
x=gmat;
%normalize the gene expression matrix
for i=1:size(x,1)
    x(i,:)=(x(i,:)-gmat_mean)./gmat_std;
end
%find the l2 norm of normalized gene expression for each gene
xnorm=zeros(1,size(x,1));
for kk=1:size(x,1)
    xnorm(1,kk)=norm(x(kk,:));
end

%generate images and htmls
tmpimg=zeros(sizeGrid);
for dsize=dictSize
    clear amat;
    clear dmat;
    disp(dsize);
    name=[basename '_dsize_' num2str(dsize) '_lambda_150_iter_1000_Amat'];
    mkdir(['images\gene_' name]);
    %data
    amat=load(['results\' name '.txt']);
    dmat=load(['results\' basename '_dsize_' num2str(dsize) '_lambda_150_iter_1000_Dmat.txt']);

    tmp=amat;
    tmp(amat<2)=0;
    dexpress=tmp*gmat';
    dexpress=dexpress./repmat(sum(tmp,2),[1 size(dexpress,2)]);
    [dexpress_sort,dexpress_p]=sort(dexpress);
    
    for i=1:dsize
        disp([dsize i]);
        c = clock;
        disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
        %generate image of gene expressions
        [ignore,p1]=sort(dmat(:,i),'descend');
        [ignore,p2]=sort(amat(i,amat(i,:)>5));
        clear plotmat ignore
        plotmat=gmat(:,amat(i,:)>5)';
        plotmat=plotmat(p2,:);
        imwrite(plotmat,hot,['images\gene_' name '\sorted_' num2str(i) '.png']);
        
        %find the lists of genes for highly/low expressed
        generate=(dmat(:,i)-meand)./stdd;
        geneList=section_data_id(generate<-3);
        tmp=generate(generate<-3);
        [ignore p]=sort(tmp);
        geneList=geneList(p);
        geneList_r=section_data_id(generate>3);
        tmp=generate(generate>3);
        [ignore p]=sort(tmp,'descend');
        geneList_r=geneList_r(p);
        dspan=2;
        tmp=[];
        for rank=1:dspan
            tmp=[tmp,find(dexpress_p(rank,:)==i)];
        end
        geneList=section_data_id(tmp);
        tmp=[];
        for rank=dsize:-1:dsize-dspan+1
            tmp=[tmp,find(dexpress_p(rank,:)==i)];
        end
        geneList_r=section_data_id(tmp);
        
        %generate image for annotation
        %find the slice of major expression
        tmpimg(:)=0;
        tmpimg(master_mask)=amat(i,:);
        tmp=sum(sum(tmpimg,3),2);
        [ignore,x_atlas]=max(tmp);
        tmp=sum(tmpimg(x_atlas,:,:),3);
        y_1=min(find(tmp>2));
        y_2=max(find(tmp>2));
        y_atlas=(y_1+y_2)/2;
        height=y_2-y_1;

        tmp=sum(tmpimg(x_atlas,:,:),2);
        z_1=min(find(tmp>2));
        z_2=max(find(tmp>2));
        z_atlas=(z_1+z_2)/2;
        width=z_2-z_1;

        plotslice=tmpimg(x_atlas,y_1:y_2,z_2:-1:z_1);
        plotslice=reshape(plotslice,[height+1,width+1]);
        imwrite(plotslice,hot,['images\gene_' name '\slice_' num2str(i) '.png']);

        %micron resolution
        x_ref=x_atlas*res;
        y_ref=y_atlas*res;
        z_ref=z_atlas*res;
        width=width*res;
        height=height*res;
        
        imgwidth=width;
        maxwidth=350;
        downsample=0;
        while(imgwidth>maxwidth)
            imgwidth=imgwidth/2;
            downsample=downsample+1;
        end
        
        imgwidth=round(imgwidth);

        %write the first half of the webpage on structures
        fid=fopen(['images\html_' name '\' num2str(i) '.html'],'w');
        fprintf(fid,'<html><body>');
        fprintf(fid,'Spatial map visualization and covered brain regions.\n');
        fprintf(fid,'(Download: <a href="../%s/amat%d.zip">zip file</a>.)\n',name,i);
        fprintf(fid,'<div><iframe style="border:none" src="%d_anno.html" width="1500" height="500" ></iframe></div>\n',i);
        fprintf(fid,'<br/><h1>Anatomical atlas and the differentially expressed genes<h1/>\n');
        %turn a row
        fprintf(fid,'<h2>anatomy<h2/><br/>\n');
        fprintf(fid,'<table><tr>\n');
        %query the AMBA database and write the second half of webpage
        imageCoord=aibs_portal_referenceToImage(referenceSpaceId,x_ref,y_ref,z_ref,section_data_id(1));
        atlasCoord=aibs_portal_imageToAtlas(imageCoord.imageId, imageCoord.x, imageCoord.y, atlasId);
        atlasInfo=aibs_portal_loadAtlasImageInfo(atlasCoord.imageId, atlasName);
        nissleUrl=aibs_portal_makeCenteredImageUrl(downsample,atlasInfo.path,atlasCoord.x+atlasInfo.x,atlasCoord.y+atlasInfo.y,width,height);
        annoUrl=aibs_portal_makeCenteredImageUrl(downsample,atlasInfo.annoPath,atlasCoord.x+atlasInfo.x,atlasCoord.y+atlasInfo.y,width,height);
        info{1}.name=['Nissl stained image'];
        info{1}.url=['http://atlas.brain-map.org/atlas?atlas=' num2str(atlasId) '#atlas=' num2str(atlasId) '&plate=' num2str(atlasInfo.imageId)];
        info{2}.name=['ontology'];
        info{2}.url=['http://atlas.brain-map.org/atlas?atlas=' num2str(atlasId) '#atlas=' num2str(atlasId) '&plate=' num2str(atlasInfo.imageId)];
        fprintf(fid,'<td><a href="%s"><img src=%s/></a></td>\n',info{1}.url,nissleUrl);
        fprintf(fid,'<td><a href="%s"><img src=%s/></a></td>\n',info{2}.url,annoUrl);
        fprintf(fid,'<td><img src=%s width=%d height=%d /></td>',['..\gene_' name '\slice_' num2str(i) '.png'],imgwidth,round(height*imgwidth/width));
        info{3}.name=['spatial map'];
        info{3}.url=['../' name '/amat' num2str(i) '.zip'];
        infoid=4;
        fprintf(fid,'</tr>\n<tr>\n');
        for k=1:infoid-1
            fprintf(fid,'<td align=center><a href="%s">%s</a></td>\n',info{k}.url,info{k}.name);
        end
        fprintf(fid,'</tr>\n</table>\n');
        
        fprintf(fid,'<br/><h2>relatively highly expressed genes<h2/><br/>\n');
        fprintf(fid,'<table><tr>\n');
        infoid=1;
        accuwidth=0;
        for geneId=geneList_r'
            accuwidth=accuwidth+imgwidth+5;
            if(accuwidth>1500)
                %turn a row
                fprintf(fid,'</tr>\n<tr>\n');
                for k=1:infoid-1
                    fprintf(fid,'<td align=center><a href="%s">%s</a></td>\n',info{k}.url,info{k}.name);
                end
                fprintf(fid,'</tr>\n<tr>\n');
                
                infoid=1;
                clear info;
                accuwidth=imgwidth+5;
            end
            
            imageCoord=aibs_portal_referenceToImage(referenceSpaceId,x_ref,y_ref,z_ref,geneId);
            imageInfo=aibs_portal_loadImageInfo(imageCoord.dataId, imageCoord.imageId);
            imageUrl=aibs_portal_makeCenteredImageUrl(downsample,imageInfo.path,imageCoord.x+imageInfo.x,imageCoord.y+imageInfo.y,width,height);
            geneInfo=aibs_portal_loadGeneInfo(geneId);
            info{infoid}.name=geneInfo.acronym;
            info{infoid}.url=['http://mouse.brain-map.org/experiment/show/' num2str(geneId)];
            fprintf(fid,'<td><a href="%s"><img src=%s/></a></td>\n',info{infoid}.url,imageUrl);
            infoid=infoid+1;
        end
        if infoid>1
            fprintf(fid,'</tr>\n<tr>\n');
                for k=1:infoid-1
                    fprintf(fid,'<td align=center><a href="%s">%s</a></td>\n',info{k}.url,info{k}.name);
                end
                fprintf(fid,'</tr>\n');
                
                infoid=1;
                clear info;
        end
        fprintf(fid,'</table>\n');
        
        fprintf(fid,'<br/><h2>relatively none expressed genes<h2/><br/>\n');
        fprintf(fid,'<table><tr>\n');
        infoid=1;
        accuwidth=0;
        for geneId=geneList'
            accuwidth=accuwidth+imgwidth+5;
            if(accuwidth>1500)
                %turn a row
                fprintf(fid,'</tr>\n<tr>\n');
                for k=1:infoid-1
                    fprintf(fid,'<td align=center><a href="%s">%s</a></td>\n',info{k}.url,info{k}.name);
                end
                fprintf(fid,'</tr>\n<tr>\n');
                
                infoid=1;
                clear info;
                accuwidth=imgwidth+5;
            end
            
            imageCoord=aibs_portal_referenceToImage(referenceSpaceId,x_ref,y_ref,z_ref,geneId);
            imageInfo=aibs_portal_loadImageInfo(imageCoord.dataId, imageCoord.imageId);
            imageUrl=aibs_portal_makeCenteredImageUrl(downsample,imageInfo.path,imageCoord.x+imageInfo.x,imageCoord.y+imageInfo.y,width,height);
            geneInfo=aibs_portal_loadGeneInfo(geneId);
            info{infoid}.name=geneInfo.acronym;
            info{infoid}.url=['http://mouse.brain-map.org/experiment/show/' num2str(geneId)];
            fprintf(fid,'<td><a href="%s"><img src=%s/></a></td>\n',info{infoid}.url,imageUrl);
            infoid=infoid+1;
        end
        if infoid>1
            fprintf(fid,'</tr>\n<tr>\n');
                for k=1:infoid-1
                    fprintf(fid,'<td align=center><a href="%s">%s</a></td>\n',info{k}.url,info{k}.name);
                end
                fprintf(fid,'</tr>\n');
                
                infoid=1;
                clear info;
        end
        fprintf(fid,'</table>\n');
        % write the part for 'genes that use the dictionary'
        error=dmat(:,i).*amat(i,:);
        error_abs_sum=sum(error.^2,2);
        error_abs_sum=sqrt(error_abs_sum);
        error_abs_sum=error_abs_sum./f_norm';
        [val,idx]=sort(error_abs_sum,'descend');
        max_size=sum(val>0.1);
        if max_size<10
            max_size=sum(val>0.05);
        end
        if max_size>50
            max_size=50;
        end
        %find the list of genes
        geneList=section_data_id(idx(1:max_size));
        fprintf(fid,'<br/><h2>Genes that use this dictionary</h2>\n');
        fprintf(fid,'*value in the parentheses show the weight of the dictionary\n');
        fprintf(fid,'<table><tr>\n');
        infoid=1;
        accuwidth=0;
        count=0;
        for geneId=geneList'
            count=count+1;
            accuwidth=accuwidth+imgwidth+5;
            if(accuwidth>1500)
                %turn a row
                fprintf(fid,'</tr>\n<tr>\n');
                for k=1:infoid-1
                    fprintf(fid,'<td align=center><a href="%s">%s</a></td>\n',info{k}.url,info{k}.name);
                end
                fprintf(fid,'</tr>\n<tr>\n');
                
                infoid=1;
                clear info;
                accuwidth=imgwidth+5;
            end
            
            imageCoord=aibs_portal_referenceToImage(referenceSpaceId,x_ref,y_ref,z_ref,geneId);
            imageInfo=aibs_portal_loadImageInfo(imageCoord.dataId, imageCoord.imageId);
            imageUrl=aibs_portal_makeCenteredImageUrl(downsample,imageInfo.path,imageCoord.x+imageInfo.x,imageCoord.y+imageInfo.y,width,height);
            geneInfo=aibs_portal_loadGeneInfo(geneId);
            info{infoid}.name=[geneInfo.acronym,' (',num2str(val(count),'%5.3f'),')',];
            info{infoid}.url=['http://mouse.brain-map.org/experiment/show/' num2str(geneId)];
            fprintf(fid,'<td><a href="%s"><img src=%s/></a></td>\n',info{infoid}.url,imageUrl);
            infoid=infoid+1;
        end
        if infoid>1
            fprintf(fid,'</tr>\n<tr>\n');
                for k=1:infoid-1
                    fprintf(fid,'<td align=center><a href="%s">%s</a></td>\n',info{k}.url,info{k}.name);
                end
                fprintf(fid,'</tr>\n');
                
                infoid=1;
                clear info;
        end
        fprintf(fid,'</table>\n');
        fclose(fid);
    end
end