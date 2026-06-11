#Functions
units = Units(:meV, :angstrom)

function oned(sys, qq, energies, eres)
    # Plot 1d cut for constant q (qq), in the energy range (energies). The energy resolution is (eres).
    swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
    kernel = gaussian(σ=eres)
    res = intensities(swt, qq; energies, kernel)

    fig = Figure()
    ax = Axis(fig[1, 1]; 
        xlabel="E [meV]", 
        ylabel="Intensity [arb. units]", 
        title="q = $(qq[1])" # Convert to string
    )

    lines!(res.energies, res.data[:, 1]; label="B=0T")

    ##
    B=3.5
    set_field!(sys, B/sqrt(1^2+0.9107^2)*[1,-0.9107,0] * units.T)

    sys_minb = reshape_supercell(sys, [1 1 1; -1 1 0; 0 0 1])
    randomize_spins!(sys_minb)
    minimize_energy!(sys_minb, maxiters=100000);

    swtb = SpinWaveTheory(sys_minb; measure=ssf_perp(sys_minb))

    res = intensities(swtb, qq; energies, kernel) #kT=1.380649*10^(-23)*290
    lines!(res.energies, res.data[:, 1]; label="B=3.5T")
    axislegend()
    fig
end


function colormap_sw()
    # Define modified version of colormap with white at the middle (to ressemble SpinW)
    base_cmap = get(ColorSchemes.inferno, LinRange(0, 1, 256))

    N = 256
    zero_index = Int(N)

    cmap = copy(base_cmap)
    cmap[zero_index] = RGB(1,1,1)  # Set white at midpoint (assumes data range is symmetric)
    return Reverse(cmap)
end


function disp_band(sys; titl=" ", ylim=[0,100], qs=[[2.5,2.5,0.5],[3,3,0],[2.5,2.5,0],[2.5,2.5,1]], colormax=150, Eres=0.26,cryst=cryst)
    # Plot dipersion, all modes are drawn as lineplot.
    mswt = SpinWaveTheory(sys; measure=ssf_perp(sys))
    path = q_space_path(cryst, qs, 500)

    res = intensities_bands(swt, path) # Simple lineplot
    plot_intensities(res; units, colormap=colormap_sw(), colorrange=(0,colormax),fwhm=Eres, ylims=(ylim[1],ylim[2]), title="Line dispersion: $(titl), fwmh=$(Eres) eV")
end 


function disp_neutron(sys, T, energies, Eres, titl, qs=[[1,-2,-1],[1,-2,1]], colormax=150, cryst=cryst)
    # Plot neutron scattering intensities, for a given temperature (T -> uses kb and scales only intensities), energy resolution (Eres). 
    swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
    path = q_space_path(cryst, qs, 700)
    kernel = gaussian(fwhm=Eres) #to fwh!!!

    Temp=T #only scales intensities (by qunatum thermal fluctuations(nB))
    res = intensities(swt, path; energies, kernel,kT=8.617333262*10^(-2)*Temp) # Boltmzann constant in meV/K
  
    title = titl
    

    plot_intensities(res; units, colorrange=(0,colormax),colormap=colormap_sw(), title)
end

function angles_crystplane(sys)
    # Print angles of spin-dipoles out of the crystal plane.
    fpath = "/Users/fridabirkedalnielsnen/Desktop/Sunny_vec.txt"

    open(fpath, "w") do f
        println(f, join(["a", "b", "c"], ", "))
    end

    for i in 1:16
        a=sys.dipoles[1,1,1,i][1]
        b=sys.dipoles[1,1,1,i][2]
        c=sys.dipoles[1,1,1,i][3]
        open(fpath, "a") do f
            println(f, join([a, b, c],","))
        end
        println("$(i) vector: ($(a),$(b),$(c)), ang_OP: $(atand(c/sqrt(a^2+b^2))) deg")
    end
end

function magnetization(sys, b=[0.1,20], dir=[[1, -1, 0]], fpath="/Users/fridabirkedalnielsnen/Desktop/mag_curve.txt", line_dot=[false],label="")
    # Calculate and plot magnetization (in range b) along selected directions (dir), save b values and m values in .txt file (fpath)
    Bs = b[1]:b[1]:b[2]
    count=0

    open(fpath, "w") do f
        println(f, join(Bs, ", "))
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; 
        xlabel="B [ T ]", 
        ylabel="M [ μB ]", 
        title="Magnetization$(label)")

    for d in dir
        count+=1
        M = Float64[]
        for B in Bs
            field_dir = [d[1]*7.9711, d[2]*8.4391, d[3]*6.0003]
            norm_field = B*field_dir / norm(field_dir)
            set_field!(sys, norm_field * units.T)
            sys_min = reshape_supercell(sys, [1 1 1; -1 1 0; 0 0 1])
            randomize_spins!(sys_min)
            minimize_energy!(sys_min, maxiters=100000)

            mval = sum(-S·norm_field/B for S in sys_min.dipoles)*2 #·2: (lande g factor) (pr spin)

            println("B = $(B) T, M = $(mval)")
            push!(M, mval)
        end

        col = get(ColorSchemes.jet1, count / length(dir))  # `total_count` is the number of total lines you're plotting

        lines!(ax, collect(Bs), M;
            label = "H ‖ $(d)",
            linestyle = line_dot[count] ? :dash : :solid,
            color = col
        )

        open(fpath, "a") do f
            println(f, join(M, ", "))
        end

    end

    axislegend(ax, position=(0,1))
    return fig
end

function save_int(sys; fpath="/Users/fridabirkedalnielsnen/Desktop/spaghetti_sunny_apr_B0.txt",size=[500,1000], Elim=[0,100], Eres=0.26, qs=[[2.5,2.5,0.5],[3,3,0],[2.5,2.5,0],[2.5,2.5,1]])
    # Save .txt file with intensities of spinwave colorplot
    swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
    path = q_space_path(cryst, qs, size[1])

    res = intensities(swt, path; energies=range(Elim[1],Elim[2],size[2]), kernel=gaussian(fwhm=Eres), kT=0)

    open(fpath, "w") do f
        for row in eachrow(res.data)
            println(f, join(row, ", "))
        end
    end
end

function save_res(res; fpath="/Users/fridabirkedalnielsnen/Desktop/LLD_0K.txt")
    open(fpath, "w") do f
        for row in eachrow(res.data)
            println(f, join(row, ", "))
        end
    end
end

function qticks(N, qs, start, npts, dir)
    tick_indices = round.(Int, range(start, npts, length=N))

    # Interpolate the q-coordinates at those indices
    q_start = qs[1] 
    q_end   = qs[2]

    if dir == 'h'
        tick_labels = ["$(round(q_start[1] + t*(q_end[1]-q_start[1]), digits=2))" for t in range(0, 1, length=N)]

    elseif dir == 'k'
        tick_labels = ["$(round(q_start[2] + t*(q_end[2]-q_start[2]), digits=2))" for t in range(0, 1, length=N)]

    elseif dir == 'l'
        tick_labels = ["$(round(q_start[3] + t*(q_end[3]-q_start[3]), digits=2))" for t in range(0, 1, length=N)]
    else
        println("Invalid direction. Use 'h', 'k', or 'l'")
        tick_labels = []
    end

    return tick_indices, tick_labels
end

function qticks_lang(N, path, npts, dir)
    tick_indices = round.(Int, range(1, npts, length=N))
    q_start = path[1]
    q_end   = path[end]

    tick_labels = if dir == "h"
        ["$(round(q_start[1] + t*(q_end[1]-q_start[1]), digits=2))" for t in range(0, 1, length=N)]
    elseif dir == "k"
        ["$(round(q_start[2] + t*(q_end[2]-q_start[2]), digits=2))" for t in range(0, 1, length=N)]
    elseif dir == "l"
        ["$(round(q_start[3] + t*(q_end[3]-q_start[3]), digits=2))" for t in range(0, 1, length=N)]
    else
        println("Invalid direction. Use \"h\", \"k\", or \"l\"")
        String[]
    end

    return tick_indices, tick_labels
end
