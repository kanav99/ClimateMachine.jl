using NCDatasets
using OrderedCollections

struct NetCDFWriter <: AbstractWriter end

# NB for some reason, Dataset() saves nc files with reverse dimensions than how set
# check using ncdump of single files or aggregation (this is definitely a NCDataset thing)
# If time != 1st dim then paraview (and like programs) doesn't recognise time

function write_data(nc::NetCDFWriter, filename, dims, varvals, simtime)
    Dataset(full_name(nc, filename), "c") do ds
        # define spatial and time dimensions
        for (dn, (dv, da)) in dims
            defDim(ds, dn, length(dv) )
        end
        defDim(ds, "time", Inf ) # Inf sets UNLIMITED dimension

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
