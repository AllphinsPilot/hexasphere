# hexasphere
A module to create (almost) hexagonal grids on a sphere

## description of the grid

### base polyhedron

The grid is built by subdividing a base polyhedron into hexagon tiles.
The base polyhedron is an icosahedron, which comprises:
- 20 faces (equilateral triangles)
- 12 vertices
- 30 edges

At the vertices of the icosahedron, instead of an hexagon, a pentagon is constructed.

For a given grid resolution `n`, each edge of the icosahedron goes through exactly `n` hexagon centers (not counting the polygons found at the ends). The inside of faces is then filled with tile-centers by following a triangular pattern.
This allows to cover the icosahedron. Its surface is then mapped to the sphere using a projection.

![icogrid](https://user-images.githubusercontent.com/70936497/186472661-e6255c76-46ae-4ce8-9b13-32c0683ee48b.jpeg)
*Here, `n = 7`*

### projection

There are two available projections in `projection` module:
- `GnomonicProj`: a simple projection, which produces hexagonal tiles about 60% larger (in area) at the corners of a face than at its center.

![Gnomonic](https://user-images.githubusercontent.com/70936497/186617448-b5ba845c-f2d5-42b2-b319-a594d2041905.png)
*Principle of Gnomonic projection*

- `SnyderEAProj`: a more complex projection, slower to compute (roughly 3x slower than Gnomonic projection), but which preserves areas. The implementation is based on [Brenton R S Recht's blog](https://brsr.github.io/2021/08/31/snyder-equal-area.html). See there for more details.

### hexagon tile identifier

A tile identifier has the following pattern: `?XXXXX-YYYYY-ZZZZZ`
- `?` is one of the 20 letters `A ... T`, each letter corresponding to one face of the icosahedron
- `XXXXX`, `YYYYY`, `ZZZZZ` are the integer coordinates of the tile in the triangular mesh covering face `?`. An useful property holds:

`XXXXX + YYYYY + ZZZZZ = 2 * (n + 1)`

---

## getting started

1. Install the package with pip

`$ pip install hexasphere`

2. Import the library in python

`from hexasphere import hexgrid, projection`

## usage

### construction of a grid

1. Create a `HexGrid` object:

`my_grid = hexgrid.HexGrid()`

2. Instantiate a projection system `Projection` associated with this grid:

`my_projection = projection.MyProjection(my_grid)`

3. Provide the projection system to the grid:

`my_grid.projection = my_projection`

### playing with grid resolutions

- Compute closest grid resolution `n` for any desired hex dimension **(in kilometers)**:

```
n = my_grid.compute_n_for_radius(0.25)
n = my_grid.compute_n_for_height(0.25)
n = my_grid.compute_n_for_side(0.25)
```

- Retrieve average hex dimension **(in kilometers)** for any given resolution `n`:

```
my_grid.compute_radius_for_n(n)
my_grid.compute_height_for_n(n)
my_grid.compute_side_for_n(n)
```

### encoding and decoding

- To find the string identifier of the hexagon to which a geographic point `(lat, lon)` belongs, call:

```
hex_identifier = my_grid.latlon_to_hex(lat, lon, n, out_str=True)[0]
```

- To find the `(lat, lon)` coordinates of the center of an hex, call:

```
my_grid.hex_to_latlon(hex_identifier, in_str=True)
my_grid.hex_to_latlon(hex_identifier, n, in_str=True) # n is here not required
```

### retrieving shape data

One can also deal with an `Hexagon` object instead of an hexagon string identifier:

```
hex_object = my_grid.latlon_to_hex(lat, lon, n)[0]
hex_object = hexgrid.Hexagon(my_grid, str_id=hexagon_identifier)
```

The coordinates of the vertices of the corresponding shape can then be retrieved:

```
shape_coordinates = hex_object.retrieve_polygon(out_latlon=True)
```