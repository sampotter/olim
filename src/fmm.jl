module fmm

export march

function march(boundary::Array{Bool, 2};
               h=1, x0=1, y0=1, speed=false, method="basic")
    m, n = size(boundary)
    if !speed
        speed = ones(m, n)
    else
        error("speed function not implemented yet")
    end
    values = zeros(m, n)
    if method == "basic"
        methodnum = 0
    elseif method == "olim4mp0"
        methodnum = 1
    elseif method == "olim4rhr"
        methodnum = 2
    elseif method == "olim4rhrlut"
        methodnum = 3
    elseif method == "olim8mp0"
        methodnum = 4
    elseif method == "olim8mp1"
        methodnum = 5
    elseif method == "olim8rhr"
        methodnum = 6
    elseif method == "solim4mp0"
        methodnum = 7
    else
        error("unknown marcher method")
    end
    ccall(
        (:fmm, "../build/Release/libfmm.so"),
        Void,
        (Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Cdouble, Ptr{Cdouble}, Cint),
        values, boundary, m, n, h, speed, methodnum)
    return values
end

end
