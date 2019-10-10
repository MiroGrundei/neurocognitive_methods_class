function sfp2mat(sfppath, file_ID)
svdir='mat';
mkdir(sfppath, svdir)
infiles=fullfile(sfppath, '*.sfp');
files=dir(infiles);
% subjects={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20'};

for n=1:length(files)
    
    %subj_ID = char(regexp(files(n,1).name, '\d+', 'match'));
    
    file=fullfile(sfppath, files(n,1).name);
    
    fid=fopen(file);
    
    tmp=textscan(fid, '%s');
    tmp=tmp{1,1};
    for m=2:4
        for el=1:length(tmp)/4
            lin=(el-1)*4+m;
            matpos(el,m-1)=str2num(cell2mat(tmp(lin)));
        end
    end
    
    fidus=matpos(1:3,:);
    elpos=matpos(4:end,:);
    
    svfid=fullfile(sfppath, svdir, [file_ID  '_fid.mat']);
    save(svfid, 'fidus');
    svpos=fullfile(sfppath, svdir, [file_ID  '_pos.mat']);
    save(svpos, 'elpos');

end