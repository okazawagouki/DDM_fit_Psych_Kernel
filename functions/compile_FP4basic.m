if ismac
    % Code to run on Mac platform (Intel chip)
    optim_flags = "COPTIMFLAGS='-O3 -fwrapv -DNDEBUG'";
    eval(sprintf("mex -R2018a -v -std=c99 %s FP4basic.c", optim_flags));
    
elseif isunix
    % Code to run on Linux platform
    compile_flags = join([...
    "CFLAGS='$CFLAGS",...
    "-std=c99 -march=native -fopenmp -ftree-vectorize -ffast-math -fopt-info'",...
    ]);
    optim_flags = "COPTIMFLAGS='-O3 -fwrapv -DNDEBUG'";
    eval(sprintf("mex -R2018a -v %s %s FP4basic.c", compile_flags, optim_flags));

elseif ispc
    % Code to run on Windows platform (Intel C compiler, on sorting machine)
    compile_flags = join([...
    "COMPFLAGS='$COMPFLAGS",...
    "/arch:SKYLAKE /Qopenmp /Qopt-report-stdout'",...
    ]);
    optim_flags = "OPTIMFLAGS='/O3 /DNDEBUG'";
    eval(sprintf("mex -R2018a -v %s %s FP4basic.c", compile_flags, optim_flags));

else
    disp('Platform not supported')
end