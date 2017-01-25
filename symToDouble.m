function dblArray = symToDouble(symArray)
    dblArray = arrayfun(@double,symArray,'UniformOutput',false);
    if numel(dblArray)==1
        dblArray={double(children(symArray))};
    end
end