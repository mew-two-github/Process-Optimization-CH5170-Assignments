function f = cost2(x)
    L1 = 300; L2 = 500; L3 = 450;
    d11 = 300; d12=350;
    d21 = 200;d22 = 250;
    d31 = 150; d32 = 200;
    f = 1.2654*(x(1)*d11^1.327+(L1-x(1))*d12^1.327);
    f = f+ 1.2654*(x(2)*d21^1.327+(L2-x(2))*d22^1.327);
    f = f+ 1.2654*(x(3)*d31^1.327+(L3-x(3))*d32^1.327);
end