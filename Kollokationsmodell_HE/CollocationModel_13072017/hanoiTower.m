

N = 3;
M = 3;
T = 256;

w = N:-1:1;
p = M-1:-1:0;
towers = cell(T+1,M);

% init
t = 1;
towers{1,1} = w;

while t < T
    
    [towers, t] = moveright(towers, t, M);

    [towers, t] = moveleft(towers, t, M);
    
end


function [towers, t] = moveright(towers, t, M)

for m = 1:M-1
    source = towers{t,m};
    target = towers{t,m+1};
    
    if ~isempty(source)
        if isempty(target)
            target = source(end);
            source(end) = [];
            if isempty(source)
                clear source
                source = [];
            end
            t = t+1;
            towers(t,:) = towers(t-1,:);
            towers{t,m} = source;
            towers{t,m+1} = target;
        elseif source(end) < min(target)
            target = [target source(end)];
            source(end) = [];
            if isempty(source)
                clear source
                source = [];
            end
            t = t+1;
            towers(t,:) = towers(t-1,:);
            towers{t,m} = source;
            towers{t,m+1} = target;
        end
    end
    
    if t > 2 && isequal(towers(t,:), towers(t-2,:)) ...
            && sum(~ismember([1,3],~find(cellfun(@isempty,towers(t-1,:)))))>1
        towers(t,:) = towers(t-1,:);
        t = t-1;
    end
end
end

function [towers, t] = moveleft(towers, t, M)

for m = M:-1:2
    source = towers{t,m};
    target = towers{t,m-1};
    
    if ~isempty(source)
        if isempty(target)
            target = source(end);
            source(end) = [];
            if isempty(source)
                clear source
                source = [];
            end
            t = t+1;
            towers(t,:) = towers(t-1,:);
            towers{t,m} = source;
            towers{t,m-1} = target;
        elseif source(end) < min(target)
            target = [target source(end)];
            source(end) = [];
            if isempty(source)
                clear source
                source = [];
            end
            t = t+1;
            towers(t,:) = towers(t-1,:);
            towers{t,m} = source;
            towers{t,m-1} = target;
        end
        
    end
    
    if t > 2 && isequal(towers(t,:), towers(t-2,:)) && numel(find(cellfun(@isempty,towers(t,:))))<2
        towers(t,:) = towers(t-1,:);
        t = t-1;
    end
end
end