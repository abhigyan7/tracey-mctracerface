# Tracey McTracerFace

Tracey McTracerFace is a raytracer written in C, following the excellent Raytracing in a Weekend book. 

![Final render: Lots of small spheres, each of them either matte, glossy, or glass. Three large spheres, each of them made of one of the three materials. Blue ground](final-render.png)

To create the final render as in the book, run `make prep && make release` followed by `bin/release/trmtrf > final-render.ppm`.

Known issue: spheres with negative radii don't quite work as they need to. 

