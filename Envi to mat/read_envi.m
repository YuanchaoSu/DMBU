function [M samples lines bands] = read_envi(filename)
file_img=[filename,'.img'];
file_hdr=[filename,'.hdr'];
file_hdr_id = fopen(file_hdr,'rt');
if file_hdr_id==-1
    errordlg('lack.hdr file£¡','Error');
    return
end
input_data_hdr = fscanf(file_hdr_id,'%s');
fclose(file_hdr_id); 
samples = str2num(input_data_hdr(strfind(input_data_hdr,'samples=')+length('samples='):strfind(input_data_hdr,'lines=')-1));
lines = str2num(input_data_hdr(strfind(input_data_hdr,'lines=')+length('lines='):strfind(input_data_hdr,'bands=')-1));
bands = str2num(input_data_hdr(strfind(input_data_hdr,'bands=')+length('bands='):strfind(input_data_hdr,'headeroffset')-1));
data_type_ENVI = str2num(input_data_hdr(strfind(input_data_hdr,'datatype=')+length('datatype='):strfind(input_data_hdr,'interleave')-1));
interleave = input_data_hdr(strfind(input_data_hdr,'interleave=')+length('interleave='):strfind(input_data_hdr,'interleave=')+length('interleave=')+2);
switch data_type_ENVI
    case 1
        data_class_Matlab = 'int8';
    case 2
        data_class_Matlab = 'int16';
    case 3
        data_class_Matlab = 'int32';
    case 4
        data_class_Matlab = 'single';
    case 5
        data_class_Matlab = 'double';
    case 12
        data_class_Matlab = 'uint16';
    case 13
        data_class_Matlab = 'uint32';
    case 14
        data_class_Matlab = 'int64';
    case 15
        data_class_Matlab = 'uint64';
    otherwise
        errordlg('Data type error£¡','Error');
        return
end 
file_img_id = fopen(file_img,'r');
M = zeros(samples,lines,bands);
switch interleave
    case 'bsq'%BSQ
        for i=1:bands
            M(:,:,i)=fread(file_img_id,[samples,lines],data_class_Matlab);
        end
    case 'bil'%BIL
        for i=1:lines
            M(:,i,:)=fread(file_img_id,[samples,bands],data_class_Matlab);
            %input_data(:,i,:)=M;
        end
    case 'bip'%BIP
        for i=1:samples
            M(i,:,:)=fread(file_img_id,[bands,lines],data_class_Matlab);
            %input_data(i,:,:)=M';
        end
    otherwise
end
fclose(file_img_id); 
end