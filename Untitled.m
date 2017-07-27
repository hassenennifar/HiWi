A = [3, 2, 1];
B = [];
C = [];

% initiate call from source A to target C with auxiliary B
move(5, A, C, B)

function move(n, source, target, auxiliary)
    if n > 0
        % move n - 1 disks from source to auxiliary, so they are out of the way
        move(n - 1, source, auxiliary, target)

        % move the nth disk from source to target
        if ~isempty(source)
            target = [target source(end)];
            source(end) = [];
        end

        % Display our progress
        disp(A)
        disp(B)
        disp(C)
        disp('________')

        % move the n - 1 disks that we left on auxiliary onto target
        move(n - 1, auxiliary, target, source)

    end
end