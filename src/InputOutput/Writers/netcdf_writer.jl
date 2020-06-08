using NCDatasets
using OrderedCollections

struct NetCDFWriter <: AbstractWriter end

# we order the time dimension last because Julia is column-major (see:
# https://github.com/Alexander-Barth/NCDatasets.jl/issues/87#issuecomment-636098859)

function write_data(nc::NetCDFWriter, filename, dims, varvals, simtime)
    Dataset(full_name(nc, filename), "c") do ds
        # define spatial and time dimensions
        for (dn, (dv, da)) in dims
            defDim(ds, dn, length(dv))
        end
        defDim(ds, "time", Inf) # Inf sets UNLIMITED dimension

        # include dimensions as variables
        for (dn, (dv, da)) in dims
            defVar(ds, dn, dv, (dn,), attrib = da)
        end
        defVar(
            ds,
            "time",
            [simtime],
            ("time",),
            attrib = OrderedDict(
                "units" => "seconds since 1900-01-01 00:00:00",
                "long_name" => "time",
            ),
        )

        # save fields as variables
        for (vn, (vd, vv, va)) in varvals
            defVar(ds, vn, vv, vd, attrib = va)
        end
    end
    return nothing
end

function aggregate_files(
    nc::NetCDFWriter,
    out_dir::String,
    out_prefix::String,
    name::String,
    starttime::String,
)
    fprefix =
        joinpath(out_dir, @sprintf("%s_%s_%s_num", out_prefix, name, starttime))
    fnames = filter(
        f -> startswith(f, fprefix) && endswith(f, ".nc"),
        readdir(out_dir, join = true),
    )

    # aggregate data from a multi-file dataset
    NCDataset(fnames, "a"; aggdim = "time", deferopen = false) do mfds
        write(fprefix[1:(end - 4)] * ".nc", mfds)
    end
end
