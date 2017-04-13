%% The code writes the indexing pages that organize the components into 13 main brain regions
%it generates comp_dsize_*.html and hier_region_1.html

sizeGrid = [67 41 58];
res=200;
offset=[0 0 0];
basename='coronal';
volnum=3;
%13 regions: Isocortex, OLF, HPF, CTXsp, STR, PAL, TH, HY, MB, HB, CB, fiber tracts, VS
regionList={315 698 1089 703 477 803 549 1097 313 1065 512 1009 73};
dictionaryList=[0 100 200 400 600 800 1000];

load ./data/coronal_master_mask.mat

%load annotation info on the heirachial structures
annoinfo=[];
annoids=[];
fp=fopen('./data/structures.csv','r');
fgetl(fp);
k=1;
while ~feof(fp)
    strline=fgetl(fp);
    strbreak=[];
    strmask=0;
    for i=1:length(strline)
        if(strmask==0)
            if(strline(i)==',')
                strbreak=[strbreak i];
            end
            if(strline(i)=='"')
                strmask=1;
            end
        else
            if(strline(i)=='"')
                strmask=0;
            end
        end
    end
    annoinfo{k}.id=str2double(strline(1:strbreak(1)-1));
    annoinfo{k}.name=strrep(strline(strbreak(2)+1:strbreak(3)-1),'"','');
    annoinfo{k}.acronym=strline(strbreak(3)+1:strbreak(4)-1);
    annoinfo{k}.tree=strline(strbreak(12)+1:strbreak(13)-1);
    annoinfo{k}.color=strline(strbreak(13)+1:strbreak(14)-1);
    annoids(k)=annoinfo{k}.id;
    parent=0;
    for i=1:length(regionList)
        tmp_flag=0;
        for j=1:length(regionList{i})
            if(~isempty(strfind(annoinfo{k}.tree, ['/' num2str(regionList{i}(j)) '/'])))
                parent=regionList{i}(j);
                tmp_flag=1;
                break;
            end
        end
        if(tmp_flag)
            break;
        end
    end
    annoinfo{k}.parent=parent;
    k=k+1;
end
fclose(fp);

%name of selected regions
regionName=[];
for i=1:length(regionList)
    aidx=find(annoids==regionList{i}(1));
    regionName{i}=annoinfo{aidx}.name;
end

%load and assign the child regions to the parent regions
load ./data/annotaion_100_downsample.mat;
mask=zeros(sizeGrid);
mask(1:66,1:40,1:57)=downsampled_data;
for i=1:length(mask(:))
    if(mask(i)~=0)
        aidx=find(annoids==mask(i));
        mask(i)=annoinfo{aidx}.parent;
    end
end

%save image masks for each of 13 regions
mkdir('images\anno');
for i=1:length(regionList)
    tmp=mask==regionList{i}(1);
    for j=2:length(regionList{i})
        tmp(mask==regionList{i}(j))=1;
    end
    fid = fopen(['images\anno\amat' num2str(i) '.raw'],'w');
    fwrite(fid,abs(tmp),'uint8');
    fclose(fid);
    fid = fopen(['images\anno\amat' num2str(i) '.mhd'],'w');
    fprintf(fid,'ObjectType = Image\nNDims = 3\nBinaryData = True\n');
    fprintf(fid,'BinaryDataByteOrderMSB = False\nCompressedData = False\n');
    fprintf(fid,'TransformMatrix = 1 0 0 0 1 0 0 0 1\nOffset = %d %d %d\nCenterOfRotation = 0 0 0\n',offset(1),offset(2),offset(3)); 
    fprintf(fid,'AnatomicalOrientation = RAI\nElementSpacing = %d %d %d\n',res,res,res);
    fprintf(fid,'DimSize = %d %d %d\nElementType = MET_UCHAR\nElementDataFile = ',sizeGrid(1),sizeGrid(2),sizeGrid(3));
    fprintf(fid,['amat' num2str(i) '.raw']);
    fclose(fid);
end

%arrange images
mats_all=[];
mats_all{1}.name='anno';
mats_all{1}.mat=zeros(length(regionList),sum(master_mask))>0;
mats_all{1}.hdis=zeros(length(regionList),sum(master_mask));
for i=1:length(regionList)
    tmp=mask==regionList{i}(1);
    for j=2:length(regionList{i})
        tmp(mask==regionList{i}(j))=1;
    end
    mats_all{1}.mat(i,:)=tmp(master_mask);
    %hausdorff distance
    tmpdis=bwdist(reshape(tmp,sizeGrid));
    mats_all{1}.hdis(i,:)=tmpdis(master_mask);
end
i=2;
for k=dictionaryList(2:end)
    disp([i k]);
    name=[basename '_dsize_' num2str(k) '_lambda_150_iter_1000_Amat'];
    mat=load(['results\' name '.txt']);
    mats_all{i}.name=name;
    mats_all{i}.mat=mat>5;
    mats_all{i}.parent=cell(size(mats_all{i}.mat,1),1);
    mats_all{i}.child=cell(size(mats_all{i}.mat,1),1);
    i=i+1;
end

%find corresponding component
for i=1:length(mats_all)
    mats_all{i}.parent=cell(size(mats_all{i}.mat,1),1);
    mats_all{i}.child=cell(size(mats_all{i}.mat,1),1);
end

%distance measure and assign the component to the 13 brain regions
for i=2:length(mats_all)
    for j=1
        disp([i j]);
        for jj=1:size(mats_all{j}.mat,1)
            for ii=1:size(mats_all{i}.mat,1)
                %tmp=max(mats_all{j}.hdis(jj,mats_all{i}.mat(ii,:)));
                tmp=sum(mats_all{j}.hdis(jj,mats_all{i}.mat(ii,:))>volnum)/sum(mats_all{i}.mat(ii,:));
                if(tmp<0.1)
                    mats_all{j}.child{jj}=[mats_all{j}.child{jj}; i, ii];
                    mats_all{i}.parent{ii}=[mats_all{i}.parent{ii}; j, jj];
                end
            end
        end
    end
end

%html
for k=2:length(dictionaryList)
    fid=fopen(['images\comp_dsize_' num2str(dictionaryList(k)) '.html'],'w' );
    fprintf(fid,'<html><body style="background-color:white;">dictionary size %d<table><tbody>\n',dictionaryList(k));
    for i=1:length(regionList)
        aid=find(annoids==regionList{i}(1));
        fprintf(fid,'<tr>\n<td style="background-color:#%s;">%s</td>\n</tr>\n',annoinfo{aid}.color,regionName{i});
        fprintf(fid,'<tr>\n');
        fprintf(fid,'<td>annotation</td>\n<td><img src="screenshots_%s\\%s.res100_60.%.4d.jpg" width="100" height="60"/></td>',...
                mats_all{1}.name, mats_all{1}.name, i-1);
        fprintf(fid,'</tr>\n<tr><td>components</td>\n');
        tmpcount=0;
        for j=1:size(mats_all{1}.child{i},1)
            if(mats_all{1}.child{i}(j,1)==k)
                fprintf(fid,'<td><a href="html_%s\\%d.html">',mats_all{mats_all{1}.child{i}(j,1)}.name, mats_all{1}.child{i}(j,2));
                fprintf(fid,'<img src="screenshots_%s\\%s.res100_60.%.4d.jpg" width="100" height="60" title="%d"/></a></td>',...
                    mats_all{mats_all{1}.child{i}(j,1)}.name, mats_all{mats_all{1}.child{i}(j,1)}.name, mats_all{1}.child{i}(j,2)-1, mats_all{1}.child{i}(j,2));
                tmpcount=tmpcount+1;
                if(tmpcount==10)
                    fprintf(fid,'</tr><tr><td></td>\n');
                    tmpcount=0;
                end
            end
        end
        fprintf(fid,'<tr>\n');
    end
    fprintf(fid,'</tbody></table></body></html>\n');
    fclose(fid);
end

for i=1:length(regionList)
    fid=fopen(['images\hier_region_' num2str(i) '.html'],'w' );
    fprintf(fid,'<html><body style="background-color:white;">%s<table><tbody>\n',regionName{i});

    fprintf(fid,'<tr>\n');
    fprintf(fid,'<td>annotation</td>\n<td><img src="screenshots_%s\\%s.res100_60.%.4d.jpg" width="100" height="60" /></td>',...
            mats_all{1}.name, mats_all{1}.name, i-1);
    fprintf(fid,'</tr>\n<tr>\n');
    for j=1:size(mats_all{1}.child{i},1)
        if((j>1 && mats_all{1}.child{i}(j,1) ~= mats_all{1}.child{i}(j-1,1)) || j==1)
            fprintf(fid,'</tr>\n<tr>\n<td>dictionary size %d</td>\n',dictionaryList(mats_all{1}.child{i}(j,1)));
            tmpcount=0;
        end
        if(tmpcount==10)
            fprintf(fid,'</tr><tr><td></td>\n');
            tmpcount=0;
        end
        fprintf(fid,'<td><a href="html_%s\\%d.html">',mats_all{mats_all{1}.child{i}(j,1)}.name, mats_all{1}.child{i}(j,2));
        fprintf(fid,'<img src="screenshots_%s\\%s.res100_60.%.4d.jpg" width="100" height="60" title="%d"/></a></td>',...
            mats_all{mats_all{1}.child{i}(j,1)}.name, mats_all{mats_all{1}.child{i}(j,1)}.name, mats_all{1}.child{i}(j,2)-1,mats_all{1}.child{i}(j,2));
        tmpcount=tmpcount+1;
    end
    fprintf(fid,'<tr>\n');

    fprintf(fid,'</tbody></table></body></html>\n');
    fclose(fid);
end

%write the index.html
fid=fopen('images\index.html','w');
fprintf(fid,'<html><body">\n');
fprintf(fid,'<h1>Transcriptome architecture of adult mouse brain revealed by sparse coding of genome-wide in situ hybridization images</h1>\n');
for k=2:length(dictionaryList)
    fprintf(fid,'dictionary size: %d<br/>\n',dictionaryList(k));
    fprintf(fid,'&nbsp;&gt;<a href="comp_dsize_%d.html">anatomical classified components</a><br/>\n',dictionaryList(k));
    fprintf(fid,'&nbsp;&gt;<a href="%s_dsize_%d_lambda_150_iter_1000_Amat.html">all components, anatomical labeled</a><br/>\n',basename,dictionaryList(k));
    fprintf(fid,'<br/>\n');
end
fprintf(fid,'all dictionary sizes, annatomical classified components:<br/>\n');
for i=1:length(regionList)
    fprintf(fid,'&nbsp;&gt;<a href="hier_region_%d.html">%s</a><br/>\n',i,regionName{i});
end
fprintf(fid,'</tbody></table></body></html>\n');
fclose(fid);