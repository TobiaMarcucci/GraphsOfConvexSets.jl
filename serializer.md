# The GCS-Serialization Format

A directory is structured as so:

```
gcs_name
|-- graph.lgz
|-- vertex_programs
| |-- <vertex_id>
| | |-- program.cbf
|-- edges_programs
| |-- <edge_name>
| | |-- map_vertices.json
| | |-- program.cbf
```

## Parsing the Graph

The graph is specified in the [graphml](http://graphml.graphdrawing.org/) file format. The file **must** be named graph.graphml.

Modification: The graph is specified in the `lgz` format, as it seems to be the default one with `Graphs.jl`. But it seems not difficult to take another format (e.g. `graphml` as you proposed), there exists a `GraphIO.jl` library that could be interseting to use. So it depends on what is the easiest for `Python` and/or `C++`.

## Parsing A Vertex Program.

Every vertex is defined in a folder called `vertices_programs/<vertex_id>`. In this folder there **must** be a file called `info.json`. This file has the following structure.

```
{
    "program": "" # specifies the relative path of the file in the folder containing the program. This includes the program's file extension
    "name": "" # this is optional if you want to give the vertex a different name than it's id.
}
```

Modification: Currently I do not have an `info.json` file, but again it should not be too difficult to do it. The programs are in files called `program.cbf`.

## Parsing An Edge Program.

Every edges is defined in a folder called `edges_programs/<edge_id>`. The `<edge_id>` is always given by the id of the source and target vertices from the graphml file. For example if the edge is specified as `<edge source="n0" target="n2"/>`, then the `<edge_id>` is given by `n0_n2`. In this folder there **must** be a file called `info.json`. This file has the following structure.

```
{
    "program": "" # specifies the relative path of the file in the folder containing the program. This includes the program's file extension
    "fmt": "" # one of csfsf, cbf, or ...
    "source_variable_mapping": {
        {
            <edge_variable>: <vertex_variable>
        }
    }
    "name": "" # this is optional if you want to give the edge a different name than it's id.
}
```

Modification: Same remark as for the vertices. There is also a file called `map_vertices.json` that for each variable of the edge program gives three values: a boolean indicating if the variable belongs to a vertex or not, the index of the vertex if it was true, and the index of the variable in the vertex program (again if it was true). If the variable belongs to no vertex, we have (false, null, null).
