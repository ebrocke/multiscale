function [get add] = prealloc(ladd, size)
get = @get_;
add = @add_;
ladd_ = ladd;
size_ = size;
arr = zeros(ladd_,size_);

l_ = ladd_;
i_ = 0;
    function add_(data)
        data_len = length(data(:,1));
        if i_+data_len > l_
            enlarge()
        end
        arr(i_+1:i_+data_len,:) = data;
        i_ = i_+data_len;
    end
    function out = get_(last)
        if last
            out = arr(i_,:);
        else
            out = arr(1:i_,:); 
        end
    end
    function enlarge()
        arr = [arr; zeros(ladd_,size_)];
        l_ = l_+ladd_;
    end
end