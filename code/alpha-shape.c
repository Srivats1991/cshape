#include "alpha-shape.h"
#include "shape.h"

/*-------------------------------------------------------------------*/
//headers from qhull
#include "libqhull/libqhull.h"
#include "libqhull/poly.h"
#include "libqhull/qset.h"
#include "libqhull/geom.h"


/*-------------------------------------------------------------------*/
//GLOBAL VARIABLES used to store the triangulation result
// (defined in shape.h/.c ).
extern tVertex vertices;  // Here vertices are stored.
extern tEdge edges;       // ...  edges are stored.
extern tFace faces;       // ...  triangles are stored.
extern tTetra tetras;     // ...  tetrahedra are stored.

// This function will estimate the tetrahedra circumsphere square
// radius.
// The logic is described here: 
// http://www.ambrsoft.com/TrigoCalc/Sphere/Spher3D_.htm
double EstimateR(tVertex tetra_points[])
{
    // Define the auxiliary variables.
    double x[4]; double y[4]; double z[4]; double t[4];
    double T, D, E, F, G;  
  
    // Fill x,y,z variables by x,y,z coordinates of the tetrahedron
    // vertices.
    x[0] = tetra_points[0]->v[0];
    x[1] = tetra_points[1]->v[0];
    x[2] = tetra_points[2]->v[0];
    x[3] = tetra_points[3]->v[0];
 
    y[0] = tetra_points[0]->v[1];
    y[1] = tetra_points[1]->v[1];
    y[2] = tetra_points[2]->v[1];
    y[3] = tetra_points[3]->v[1];
 
    z[0] = tetra_points[0]->v[2];
    z[1] = tetra_points[1]->v[2];
    z[2] = tetra_points[2]->v[2];
    z[3] = tetra_points[3]->v[2];

    t[0] = -(x[0]*x[0] + y[0]*y[0] + z[0]*z[0]);
    t[1] = -(x[1]*x[1] + y[1]*y[1] + z[1]*z[1]);
    t[2] = -(x[2]*x[2] + y[2]*y[2] + z[2]*z[2]);
    t[3] = -(x[3]*x[3] + y[3]*y[3] + z[3]*z[3]);
    
    T = x[0]*( y[1] * ( z[2] * 1    -   1  * z[3] ) - z[1] * ( y[2] * 1    -   1  * y[3]) +    1 * ( y[2]*z[3] - z[2]*y[3])) -
        y[0]*( x[1] * ( z[2] * 1    -   1  * z[3] ) - z[1] * ( x[2] * 1    -   1  * x[3]) +    1 * ( x[2]*z[3] - z[2]*x[3])) +
        z[0]*( x[1] * ( y[2] * 1    -   1  * y[3] ) - y[1] * ( x[2] * 1    -   1  * x[3]) +    1 * ( x[2]*y[3] - y[2]*x[3])) -
           1*( x[1] * ( y[2] * z[3] - z[2] * y[3] ) - y[1] * ( x[2] * z[3] - z[2] * x[3]) + z[1] * ( x[2]*y[3] - y[2]*x[3]));

    D = ( t[0]*( y[1]*(z[2]*1 - 1*z[3]) - z[1]*(y[2]*1 - 1*y[3]) + 1*(y[2]*z[3] - z[2]*y[3])) -
          y[0]*( t[1]*(z[2]*1 - 1*z[3]) - z[1]*(t[2]*1 - 1*t[3]) + 1*(t[2]*z[3] - z[2]*t[3])) +
          z[0]*( t[1]*(y[2]*1 - 1*y[3]) - y[1]*(t[2]*1 - 1*t[3]) + 1*(t[2]*y[3] - y[2]*t[3])) -
             1*( t[1]*(y[2]*z[3]-z[2]*y[3]) - y[1]*(t[2]*z[3]-z[2]*t[3] )+ z[1]*(t[2]*y[3]-y[2]*t[3])) ) / T;

    E = ( x[0]*( t[1]*(z[2]*1 - 1*z[3]) - z[1]*(t[2]*1 - 1*t[3]) + 1*(t[2]*z[3] - z[2]*t[3])) -
          t[0]*( x[1]*(z[2]*1 - 1*z[3]) - z[1]*(x[2]*1 - 1*x[3]) + 1*(x[2]*z[3] - z[2]*x[3])) +
          z[0]*( x[1]*(t[2]*1 - 1*t[3]) - t[1]*(x[2]*1 - 1*x[3]) + 1*(x[2]*t[3] - t[2]*x[3])) -
             1*( x[1]*(t[2]*z[3]-z[2]*t[3]) - t[1]*(x[2]*z[3]-z[2]*x[3] )+ z[1]*(x[2]*t[3]-t[2]*x[3])) ) / T;

    F = ( x[0]*( y[1]*(t[2]*1 - 1*t[3]) - t[1]*(y[2]*1 - 1*y[3]) + 1*(y[2]*t[3] - t[2]*y[3])) -
          y[0]*( x[1]*(t[2]*1 - 1*t[3]) - t[1]*(x[2]*1 - 1*x[3]) + 1*(x[2]*t[3] - t[2]*x[3])) +
          t[0]*( x[1]*(y[2]*1 - 1*y[3]) - y[1]*(x[2]*1 - 1*x[3]) + 1*(x[2]*y[3] - y[2]*x[3])) -
             1*( x[1]*(y[2]*t[3]-t[2]*y[3]) - y[1]*(x[2]*t[3]-t[2]*x[3] )+ t[1]*(x[2]*y[3]-y[2]*x[3])) ) / T;


    G = ( x[0]*( y[1]*(z[2]*t[3] - t[2]*z[3]) - z[1]*(y[2]*t[3] - t[2]*y[3]) + t[1]*(y[2]*z[3] - z[2]*y[3])) -
          y[0]*( x[1]*(z[2]*t[3] - t[2]*z[3]) - z[1]*(x[2]*t[3] - t[2]*x[3]) + t[1]*(x[2]*z[3] - z[2]*x[3])) +
          z[0]*( x[1]*(y[2]*t[3] - t[2]*y[3]) - y[1]*(x[2]*t[3] - t[2]*x[3]) + t[1]*(x[2]*y[3] - y[2]*x[3])) -
          t[0]*( x[1]*(y[2]*z[3]-z[2]*y[3]) - y[1]*(x[2]*z[3]-z[2]*x[3] )+ z[1]*(x[2]*y[3]-y[2]*x[3])) ) / T;

    // If T is equal to zero, then all 4-point lies on the same plane. 
    // Such "plane tetrahedron" may be circumsphered only
    // by infinitly large sphere - set the radius as maximal possible 
    // value.
    if(T == 0)
        return DBL_MAX;
    else
        return 0.25*(D*D + E*E + F*F - 4*G);
}

/*-------------------------------------------------------------------*/
void AlphaShape( unsigned int alpha )
{
	//implement 3d Alpha Shape alogorithm here
	//
	// put results (vertices, edges, faces, and tetrahedra)
	// in the global structures above
	// see main.c for detail on how qhull does it
	//
//(init)==============================================================
    // Define auxiliary variables
    tVertex  ptr_v;         // This would store one vertex (buffer).
    tVertex * all_v = NULL; // This would store all vertices (local 
                            // copy of all vertices set).
    int vsize = 0;          // This would store the size of vertices
                            // set.
    int id = 0;             // ...
    int vid = 0;            // ... these would be used as iterator.
    tTetra  face;           // This would store one tetrahedra 
                            // (buffer).
    double R2;              // This would store circumsphere square 
                            // radius

    // Define global varibles for qhull
    static char * options = (char *)"delaunay QJ Pp";  // The "qhull" 
                                              // configuration string.
    int curlong, totlong;
    coordT * pt = NULL;              //  ... (input for "qhull")
    facetT *facet = NULL;            //  ...
    vertexT *vertex = NULL;          //  ...
    vertexT **vertexp = NULL;        //  ...
    facetT *neighbor, **neighborp;   //  ... all theese variables 
                         // would be filled/used in communication
                         // with "qhull" package.



    // Count number of points
    // (just loop through all vartices).
    ptr_v = vertices;
    do
    {                                 
        vsize++;
        ptr_v = ptr_v->next;
    }
    while ( ptr_v != vertices );
    
    // Allocate memory to store local variables.
    pt = (coordT*)calloc(vsize * 4, sizeof(coordT));  // This would be 
                               // used as input for "qhull": 4D points.
    all_v = (tVertex*)calloc(vsize, sizeof(tVertex));
    // Use assert to interrupt the program in case of some problem 
    // with memory allocation.
    assert(pt && all_v);
//==(init*)===========================================================

    //
    // create points in 4D (x,y,z,x^2+y^2+z^2)
    //
//==(3d->4d)==========================================================
    // Fill the local copy of all vertices set 
    // + prepare the input for "qhull"
    ptr_v = vertices;
    do
    { 
        // Here "qhull" will obtain 4D points: 3 {x,y,z} 
        // coordinates + ...                                
        pt[id++] = ptr_v->v[0];
        pt[id++] = ptr_v->v[1];
        pt[id++] = ptr_v->v[2];
        //               ...4th for x^2+y^2+z^2.
        pt[id++] = ptr_v->v[0]*ptr_v->v[0] + 
                   ptr_v->v[1]*ptr_v->v[1] +
                   ptr_v->v[2]*ptr_v->v[2];

        // Fill local copy of all vertices set.
        all_v[ptr_v->vnum] = ptr_v;
        ptr_v = ptr_v->next;
    }
    while ( ptr_v != vertices );
//==(3d->4d*)=========================================================

    //
    // compute convex hull in 4D by calling qhull
    // use flags: static char * options=(char *)"delaunay QJ Pp";
    //
//==(qhull)===========================================================
    // Using qhull.

    // Specify the sources for "qhull <-> user" communication
    qh_init_A(stdin, stdout, stderr, 0, NULL);

    // Configure the "qhull" for performing specific task specified
    // by "options" string.
    qh_initflags(options);
    // Specify the "qhull" input: pt - array of 4D points (of size 
    // given by "vsize").
    qh_init_B (pt, vsize, 4, false);
    // Run the "qhull", wait for result.
    qh_qhull();
    // Request the result / check it.
    qh_check_output();
//==(qhull*)==========================================================

    //
    //loop through all faces and call MakeNullTetra() to make a tetra 
    //remember that this is in 4D, so each face is a tetrahedron
    //
    //to fill the teta: get vertices of facet and loop through 
    // each vertex
    //use FOREACHvertex_()
    //
//==(consume_result)==================================================

    // Loop through all facets obtained from "qhull": 
    // here each "facet" is a tetrahedron
    // (looping would be made by "facet" variable ), 
    // and parse each tetrahedon into triagnles and edges.
    FORALLfacets 
    {
            
        // Get the tetrahedron points ...
        tVertex tetra_points[4];
        vid = 0;
        FOREACHvertex_(facet->vertices)
        {
            //get vertex
            tetra_points[vid++] = all_v[qh_pointid(vertex->point)];
        }
        // ... and estimate it's circumsphere square radius.
        R2 = EstimateR(tetra_points);   

        // If estimated circumsphere square radius is greater
        // then alpha - skip this tetrahedron.
        if (R2 < alpha)
        {
            face = MakeNullTetra(); // This function "MakeNullTetra"
                                // creates new tetrahedron, 
                                // insesrts it into the GLOBAL VARIABLE
                                // and output it's pointer. So here 
                                // the tetrahedron is 
                                // actually sotored in the specific
                                // GLOBAL VARIABLE.

            // The next loop will just fill the tetrahedron added 
            // previously inserted into 
            // the G.V. by appropriate vertices obtained from qhull.
            vid = 0;
            FOREACHvertex_(facet->vertices)
            {
                // (Looping by all vertices of the "facet" and fill
                //  the tetrahedron).
                // (Looping is made by "vertex" variable).
                // ("qh_pointid(vertex->point)" is equal to the index 
                //   of the input point that is stored in the 
                // "vertex").

                face->vertex[vid++] = all_v[qh_pointid(vertex->point)];
            }

            // Loop by all neighbors of the current tetrahedron - 
            // look for tangent tetrahedras sides.
            // Tanget sides  - are triangles that should be stored in
            // respective G.V.
            FOREACHneighbor_(facet)
            {
                // (Looping is made by "neighbor" variable).

                // Just as before, estimate R2: here for neighbor
                tVertex tetra_points[4];
                vid = 0;
                FOREACHvertex_(neighbor->vertices)
                {
                    tetra_points[vid++] = all_v[qh_pointid(vertex->point)];
                }
                R2 = EstimateR(tetra_points);
                
                // (Possibly, trying to avoid of multiple counting of 
                // the same "tangent side").
                // Skip those tetrahedra that has too large R2.
                if(facet<neighbor && R2 < alpha)
                {
                    // Transform the neighbor tetrahedra given by
                    // "neighbor" into 4 vertices.
                    tVertex vertices[4];
                    vid = 0;
                    FOREACHvertex_(neighbor->vertices)
                    {
                        //get vertex
                       vertices[vid++] = all_v[qh_pointid(vertex->point)];
                    }
                    // (possibly, miss-step here: why not to just use
                    // "neighbor->vertex" instead of "vertices". If
                    // the idea, that "face->vertex" and  
                    // "neighbor->vertex" would be erased, 
                    // then it have no sense cause "all_v" also 
                    // would be "free"d and moreover, the   
                    // "QHullCreateFace" uses first argument's  
                    // vertices to create triangles, edges but not the 
                    // second).
                    QHullCreateFace(face->vertex, vertices);
                }
            }
        }
    }
//==(consume_result*)=================================================


//==(clean)===========================================================
    // Free allocated previously for "temp" purposes memory.
    free(pt);
    free(all_v);
    pt = NULL;
    all_v = NULL;
    
    // Request the "qhull" to free it's allocations.
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);
//==(clean*)==========================================================
}

