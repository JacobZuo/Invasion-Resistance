function [ComRes] = CellCompare(StrCell,StrCom)

ComRes=zeros(size(StrCell));
for i=1:size(StrCell,1)
    ComRes(i)=strcmp(StrCell{i},StrCom);
end

end

