%this function rounds off numbers to the nearest thousandth
function roundnum = rounders(n)
    n = n*1e4;
    rn = round(n);
    roundnum = rn*1e-4;
end