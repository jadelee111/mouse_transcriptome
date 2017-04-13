%% Code for first half of webpage: "Spatial map visualization and covered brain regions"
% Load the coefficient matrix and visualize in 3D brain space. the percentage 
% of overlapping volume between the component and ARA was calculated and the 
% top 20 regions were tabulated 
 
%load annotation info
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
    k=k+1;
end
fclose(fp);

% load annotation images
sizeGrid = [67 41 58];
load ./data/annotaion_100_downsample.mat;
mask=zeros(sizeGrid);
mask(1:66,1:40,1:57)=downsampled_data;

load coronal_master_mask.mat %the reference mask after preprocessing

for k=[1 2 4 6 8 10]
    name=['coronal_dsize_' num2str(k) '00_lambda_150_iter_1000_Amat'];
    folder=['images\html_' name];
    mkdir(folder);
    D=load(['results\' name '.txt']); %load coefficient matrix
    %compare with reference atlas
    masked_anno=mask(master_mask);
    anno_idx=unique(masked_anno);
    numOfAnno=length(anno_idx);

    D_mask=abs(D)>5;
    overlap=zeros(size(D,1),numOfAnno);
    a_size=sum(D_mask,2);
    m_size=zeros(numOfAnno,1);
    for j=1:numOfAnno
        this_anno_mask=masked_anno==anno_idx(j);
        m_size(j)=sum(this_anno_mask);
        for i=1:size(D,1)
            overlap(i,j)=sum(D_mask(i,this_anno_mask));
        end
    end
    ratio_a=overlap./repmat(a_size,[1 numOfAnno]);
    ratio_m=overlap./repmat(m_size',[size(D,1) 1]);

    for i=1:(k*100) %size(D,1)
        imagename=sprintf('..\\screenshots_%s\\%s.%.4d.jpg',name,name,i-1);
        [s,p]=sort(overlap(i,:),'descend');
        fid=fopen([folder '\' num2str(i) '_anno.html'],'w');
        fprintf(fid,'<html><body style="background-color:white;"><table><tbody>\n');
        fprintf(fid,'<tr style="background-color:silver;" ><td rowspan="23"> <img src="%s" width="800" height="480" /><td>voxel#</td><td>region%%</td><td width=500px>Name</td></tr>\n',imagename);
        for j=1:20
            if(s(j)<1)
                break;
            end
            idx=p(j);
            aid=anno_idx(idx);
            if(aid<=0)
                continue;
            end
            aidx=find(annoids==aid);
            fprintf(fid,'<tr style="background-color:#%s;"><td>%d</td><td>%.3f</td><td>%s</td></tr>\n',...
                annoinfo{aidx}.color,s(j),ratio_m(i,idx),annoinfo{aidx}.name);
        end
        fprintf(fid,'</tbody></table>\n');
        fprintf(fid,'</body></html>\n');
        fclose(fid);
    end
    %create the indexing page for each dictionary size
    fid=fopen(['images\' name '.html'],'w');
    for i=1:k*100
        fprintf(fid,'<div>%d, <a href="html_%s\\%d.html">link: gene analysis</a></br><iframe style="border:none" src="html_%s/%d_anno.html" width="1500" height="500" ></iframe></div>\n',i,name,i,name,i);
    end
    fclose(fid);
end