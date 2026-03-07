function x = proj_and_repair(x, xmin, xmax, dmin)
x = sort(x(:).');
x = min(max(x, xmin), xmax);

% forward spacing
for i=2:numel(x)
    if x(i) < x(i-1) + dmin
        x(i) = x(i-1) + dmin;
    end
end

% shift back if overflow
if x(end) > xmax
    shift = x(end) - xmax;
    x = x - shift;
    x = max(x, xmin);
    x = sort(x);
    for i=2:numel(x)
        if x(i) < x(i-1) + dmin
            x(i) = x(i-1) + dmin;
        end
    end
end

x = min(max(x, xmin), xmax);
end


