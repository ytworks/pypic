"""XDMF sidecar for 3-D viewing of field dumps in ParaView / VisIt.

ParaView reads the HDF5 data through the XDMF file. Dumps written with the default zstd codec
need the HDF5 filter plugins: ``export HDF5_PLUGIN_PATH=$(python -c "import hdf5plugin;
print(hdf5plugin.PLUGIN_PATH)")`` before starting ParaView (or write with output.codec: deflate).
XDMF names the slowest array axis Z: pypic's x axis appears as ParaView's Z axis.
"""
from __future__ import annotations

import os

import h5py

from ..output.names import expand_fields, record_of


def write_xdmf(run, diag: str, path: str | None = None) -> str:
    fs = run.fields(diag)
    comps = expand_fields(fs.components)
    out = path or os.path.join(run.path, f"{diag}.xdmf")
    lines = ['<?xml version="1.0" ?>', '<Xdmf Version="3.0">', "<Domain>",
             f'<Grid Name="{diag}" GridType="Collection" CollectionType="Temporal">']
    for it in fs.iterations:
        fn = os.path.relpath(fs.files[it], os.path.dirname(os.path.abspath(out)))
        with h5py.File(fs.files[it], "r") as f:
            g = f[f"data/{it}"]
            t = float(g.attrs["time"])
            rec0, c0 = record_of(comps[0])
            m = g["meshes"][rec0]
            ds = m[c0] if c0 else m
            shape = ds.shape
            sp, off = m.attrs["gridSpacing"], m.attrs["gridGlobalOffset"]
        dims = " ".join(str(n) for n in shape)
        lines += [f'<Grid Name="it{it}" GridType="Uniform"><Time Value="{t}"/>',
                  f'<Topology TopologyType="3DCoRectMesh" Dimensions="{dims}"/>',
                  '<Geometry GeometryType="ORIGIN_DXDYDZ">',
                  f'<DataItem Dimensions="3" Format="XML">{off[0]} {off[1]} {off[2]}</DataItem>',
                  f'<DataItem Dimensions="3" Format="XML">{sp[0]} {sp[1]} {sp[2]}</DataItem>', "</Geometry>"]
        for c in comps:
            rec, cc = record_of(c)
            h5path = f"/data/{it}/meshes/{rec}" + (f"/{cc}" if cc else "")
            lines += [f'<Attribute Name="{c.replace(":", "_")}" AttributeType="Scalar" Center="Node">',
                      f'<DataItem Dimensions="{dims}" NumberType="Float" Precision="4" Format="HDF">'
                      f"{fn}:{h5path}</DataItem>", "</Attribute>"]
        lines.append("</Grid>")
    lines += ["</Grid>", "</Domain>", "</Xdmf>"]
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    return out
