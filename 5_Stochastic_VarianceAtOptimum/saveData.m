function [] = saveData(filename,X,labels,format,flag)

if nargin < 4 || isempty(format)
    format = '%6.4e';
end

if nargin < 5 || flag
    f = fopen(filename,'w');
    if nargin >= 3
        for i = 1:size(X,2)
            fprintf(f,[labels{i} '\t']);            
        end
        fprintf(f,'\n');
    end
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            fprintf(f,[format '\t'],X(i,j));
        end
        fprintf(f,'\n');
    end
    fclose(f);
end

end