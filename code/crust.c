#include "crust.h"
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

// This function will estimate the circumsphere center for given by 4points
// tetrahedra
// (it uses the same idea as "EstimateR" in "alpha-shape.c").
coordT* EstimateCenter(tVertex tetra_points[])
{
    // Define the auxiliary variables.
    double x[4]; double y[4]; double z[4]; double t[4];
    double T, D, E, F, G;  
  
    coordT* center = (coordT*)calloc(3, sizeof(coordT));

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


    center[0] = -D/2; center[1] = -E/2; center[2] = -F/2;

        return center;
}

// This function will check if 4 given points lie on one plane
// (it uses the same idea as "EstimateR" in "alpha-shape.c").
int IsOnPlane(tVertex tetra_points[])
{
    // Define the auxiliary variables.
    double x[4]; double y[4]; double z[4];
    double T;

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

    T = x[0]*( y[1] * ( z[2] * 1    -   1  * z[3] ) - z[1] * ( y[2] * 1    -   1  * y[3]) +    1 * ( y[2]*z[3] - z[2]*y[3])) -
        y[0]*( x[1] * ( z[2] * 1    -   1  * z[3] ) - z[1] * ( x[2] * 1    -   1  * x[3]) +    1 * ( x[2]*z[3] - z[2]*x[3])) +
        z[0]*( x[1] * ( y[2] * 1    -   1  * y[3] ) - y[1] * ( x[2] * 1    -   1  * x[3]) +    1 * ( x[2]*y[3] - y[2]*x[3])) -
           1*( x[1] * ( y[2] * z[3] - z[2] * y[3] ) - y[1] * ( x[2] * z[3] - z[2] * x[3]) + z[1] * ( x[2]*y[3] - y[2]*x[3]));

    if(T==0)
        return 1;
    else
        return 0;
}

// This function will estimate the convex hull for vertices stored in
// "vertices" and find out which vertices belong to convex hull
// (or which vertices form conhex hull).
// The output is the array of averaged outer triangles normals in case
// if point belong to convex hull and NULL otherwise.
coordT ** ConvexHullBelonging()
{
//==(CHB_init)========================================================
    // Define auxiliary variables
    tVertex  ptr_v;         // This would store one vertex (buffer).
    tVertex * all_v = NULL; // This would store all vertices (local
                            // copy of all vertices set).
    int vsize = 0;          // This would store the size of vertices
                            // set.
    int id = 0;             // ...
    int b_id = 0;           // ...
    int vid = 0;            // ... these would be used as iterator.
    coordT ** ConvexBelongingNormal;  // Will store the result.

    // Define global varibles for qhull
    static char * options = (char *)"qhull QJ Pp";  // The "qhull" 
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
    pt = (coordT*)calloc(vsize * 3, sizeof(coordT));  // This would be 
                               // used as input for "qhull": 3D points.
    all_v = (tVertex*)calloc(vsize, sizeof(tVertex));
    ConvexBelongingNormal  = (coordT**)calloc(vsize, sizeof(coordT*)); 

    // Use assert to interrupt the program in case of some problem 
    // with memory allocation.
    assert(pt && all_v && ConvexBelongingNormal);
//==(CHB_init*)=======================================================

//==(CHB_construct points)============================================

    // Fill the local copy of all vertices set 
    // + prepare the input for "qhull"
    ptr_v = vertices;
    do
    { 
        // Here "qhull" will obtain 3D points: 3 {x,y,z} 
        // coordinates.                              
        pt[id++] = ptr_v->v[0];
        pt[id++] = ptr_v->v[1];
        pt[id++] = ptr_v->v[2];

        ConvexBelongingNormal[b_id++] = NULL;// By default all vertex                                  
                            //counted as not belonged to convex hull.
        

        // Fill local copy of all vertices set.
        all_v[ptr_v->vnum] = ptr_v;
        ptr_v = ptr_v->next;
    }
    while ( ptr_v != vertices );

//==(CHB_construct points*)===========================================

    //
    // compute convex hull in 3D by calling qhull
    //
//==(CHB_qhull)=======================================================
// Using qhull.

    // Specify the sources for "qhull <-> user" communication
    qh_init_A(stdin, stdout, stderr, 0, NULL);

    // Configure the "qhull" for performing specific task specified
    // by "options" string.
    qh_initflags(options);
    // Specify the "qhull" input: pt - array of 3D points (of size 
    // given by "vsize").
    qh_init_B (pt, vsize, 3, false);
    // Run the "qhull", wait for result.
    qh_qhull();
    // Request the result / check it.
    qh_check_output();
//==(CHB_qhull*)======================================================

//==(CHB_consume_result)==============================================
    
    // Loop through all facets obtained from "qhull": 
    // here each "facet" is a triangle.
    FORALLfacets 
    {    
        // Loop through each vertex of facet and mark it as such
        // that belongs to convex hull
        vid = 0;
        FOREACHvertex_(facet->vertices)
        {
            int neigh_counter = 0;
            ConvexBelongingNormal[b_id] = (coordT*)calloc(3, sizeof(coordT));
            
            qh_vertexneighbors();
            FOREACHneighbor_(vertex)
            {
                ConvexBelongingNormal[b_id][0] += neighbor->normal[0];
                ConvexBelongingNormal[b_id][1] += neighbor->normal[1];
                ConvexBelongingNormal[b_id][2] += neighbor->normal[2];
                neigh_counter++;
            }

            ConvexBelongingNormal[b_id][0] /= neigh_counter;
            ConvexBelongingNormal[b_id][1] /= neigh_counter;
            ConvexBelongingNormal[b_id][2] /= neigh_counter;
        }
        
    }
//==(CHB_consume_result*)=============================================


//==(CHB_clean)=======================================================
    // Free allocated previously for "temp" purposes memory.
    //free(pt);
    //free(all_v);
    pt = NULL;
    all_v = NULL;
   
    // Request the "qhull" to free it's allocations.
    //qh_freeqhull(!qh_ALL);
    //(&curlong, &totlong);
//==(CHB_clean*)======================================================

    return ConvexBelongingNormal;
}


// Structure to hold the pole's list.
typedef struct pole
{
    coordT * coord;
    struct pole* next;        
} pole_t;

// This function will find all poles: negative + positive.
pole_t* find_poles()
{
//==(pole_init)=======================================================
    // Define auxiliary variables
    tVertex  ptr_v;         // This would store one vertex (buffer).
    tVertex * all_v = NULL; // This would store all vertices (local
                            // copy of all vertices set).
    int vsize = 0;          // This would store the size of vertices
                            // set.
    int id = 0;             // ...
    int b_id = 0;           // ...
    int vid = 0;            // ... these would be used as iterator.
    tTetra  face;           // This would store one tetrahedra 
                            // (buffer).
    coordT** ConvexBelongingNormal;   // Will store the array of 
                            // averaged outer triangles normals
                            // in case if semi-indexed point lies
                            // on convex hull and NULL othervse.   

    coordT * center = NULL;

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

//==(pole_init*)======================================================

    // Find out which points form the convex hull 
    // (which lie on it and which not lie)
    // (It will be NULL if not lie, and averaged outer triangles 
    // normals if lie)
    ConvexBelongingNormal = ConvexHullBelonging();

    //
    // create points in 4D (x,y,z,x^2+y^2+z^2)
    //
//==(pole_3d->4d)====================================================
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
//==(pole_3d->4d*)====================================================

//==(pole_qhull)======================================================
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
//==(pole_qhull*)====================================================

//==(pole_consume_result)=============================================
    
    // Define variable that will store p+ poles.
    pole_t* p_plus = NULL;
    pole_t* p_plus_temp = NULL;

    // Define the variable that will store p- poles.
    pole_t* p_minus = NULL;
    pole_t* p_minus_temp = NULL;

    // Define the temp variable to store the n+ vector.
    coordT* n_vector = (coordT*)calloc(3, sizeof(coordT));

    // Loop through all points (through S)
    for(b_id = 0; b_id < vsize; b_id++)
    {
        coordT * center_max = NULL;
        double distance_max = 0;
        double distance = 0;
 
        // If the point dose not lie on the convex hull...
        if(ConvexBelongingNormal[b_id] == NULL)
        {
            // Loop through all Delaunay tetrahedr 
            // (Voronoi vertex == centers of tetrahedr circumsphere)
            FORALLfacets 
            {
                // Get the points (4 points) the form this tetrahedr
                tVertex tri_points[4];
                vid = 0;
                FOREACHvertex_(facet->vertices)
                {
                    tri_points[vid++] = all_v[qh_pointid(vertex->point)];
                }
                // Check if this tetrahedr has non-zero volume
                if(!IsOnPlane(tri_points))
                {
                    // If so, evaluate the Voronoi vertex (tetrahedr circumsphere center)                        
                    //center = qh_facetcenter(facet->vertices);
                    center = qh_getcentrum(facet);
                    //center = EstimateCenter(tri_points);                

                    // and, estimate the distance between this tetrahedr and that vertix:
                    distance = (center[0]-all_v[b_id]->v[0])*(center[0]-all_v[b_id]->v[0])+
                               (center[1]-all_v[b_id]->v[1])*(center[1]-all_v[b_id]->v[1])+
                               (center[2]-all_v[b_id]->v[2])*(center[2]-all_v[b_id]->v[2]);
                    // Upate the information about the farthest Voronoi point.
                    if (distance > distance_max)
                    {
                        center_max = center;
                        distance_max = distance;
                    }
                }
            }

            // If the farthest Voronoi point was found: estimate n vector, and store
            // the farthest point as p+ (positive pole).
            if(center_max != NULL)
            {
                // Store the n vector values (n+ = sp)
                n_vector[0] = center_max[0] - all_v[b_id]->v[0];
                n_vector[1] = center_max[1] - all_v[b_id]->v[1];
                n_vector[2] = center_max[2] - all_v[b_id]->v[2];

                // Store the pole value (p+)
                if(p_plus_temp==NULL)
                {
                    p_plus = (pole_t*)malloc(sizeof(pole_t));
                    p_plus->coord = (coordT*)calloc(3, sizeof(coordT));
                    (p_plus->coord)[0] = center_max[0];
                    (p_plus->coord)[1] = center_max[1];
                    (p_plus->coord)[2] = center_max[2];
                    p_plus->next = NULL;
                
                    p_plus_temp = p_plus;
                }
                else
                {
                    p_plus_temp->next = (pole_t*)malloc(sizeof(pole_t));
                    p_plus_temp = p_plus_temp->next;

                    p_plus_temp->coord = (coordT*)calloc(3, sizeof(coordT));
                    (p_plus_temp->coord)[0] = center_max[0];
                    (p_plus_temp->coord)[1] = center_max[1];
                    (p_plus_temp->coord)[2] = center_max[2];
                    p_plus_temp->next = NULL;
                }
            }

        }
        else
        {
            // If point lie on convex hull - n vector is counted as
            // averaged outer triangles normals.
            n_vector[0] = ConvexBelongingNormal[b_id][0];
            n_vector[1] = ConvexBelongingNormal[b_id][1];
            n_vector[2] = ConvexBelongingNormal[b_id][2];            
        }

        center_max = NULL;
        distance_max = 0;
        distance = 0;
        // Loop through all Delaunay tetrahedr (Voronoi vertex == centers of tetrahedr
        // circumsphere)
        FORALLfacets 
        {
            // Get the points (4 points) the form this tetrahedr
            tVertex tri_points[4];
            vid = 0;
            FOREACHvertex_(facet->vertices)
            {
                tri_points[vid++] = all_v[qh_pointid(vertex->point)];
            }
            // Check if this tetrahedr has non-zero volume
            if(!IsOnPlane(tri_points))
            {
                // If so, evaluate the Voronoi vertex (tetrahedr circumsphere center)                        
                //center = qh_facetcenter(facet->vertices);
                center = qh_getcentrum(facet);
                //center = EstimateCenter(tri_points);

                // and, estimate the distance between this tetrahedr and that vertix.
                distance = (center[0]-all_v[b_id]->v[0])*(center[0]-all_v[b_id]->v[0])+
                           (center[1]-all_v[b_id]->v[1])*(center[1]-all_v[b_id]->v[1])+
                           (center[2]-all_v[b_id]->v[2])*(center[2]-all_v[b_id]->v[2]);
                
                // And also estimate the cosine of the angle between n vector and vector
                // pointed from currently viewed Voronoi point to the point of interest
                // ( (n+)^(ps)).
                double scalar_product = (all_v[b_id]->v[0] - center[0]) * n_vector[0] + 
                                        (all_v[b_id]->v[1] - center[1]) * n_vector[1] +
                                        (all_v[b_id]->v[2] - center[2]) * n_vector[2];

                double ps_norm = (all_v[b_id]->v[0] - center[0]) * (all_v[b_id]->v[0] - center[0]) + 
                                 (all_v[b_id]->v[1] - center[1]) * (all_v[b_id]->v[1] - center[1]) +
                                 (all_v[b_id]->v[2] - center[2]) * (all_v[b_id]->v[2] - center[2]);
                ps_norm = sqrt(ps_norm);

                double n_norm = n_vector[0] * n_vector[0] + 
                                n_vector[1] * n_vector[1] +
                                n_vector[2] * n_vector[2];
                n_norm = sqrt(n_norm);                

                double cos_n_ps = scalar_product / (ps_norm * n_norm);
                
                // Upate the information about the farthest Voronoi point - count only those
                // points that lie very close to the back projection of n+ vector - that
                // is defined by the estimated cosine.
                if ((distance > distance_max) && (cos_n_ps >=-0.1) && (cos_n_ps <=0.1))
                {
                    center_max = center;
                    distance_max = distance;
                }
            }
        }

        // Store the negative pole if so was found.
        if(center_max != NULL)
        {
            if(p_minus_temp==NULL)
            {
                p_minus = (pole_t*)malloc(sizeof(pole_t));
                p_minus->coord = (coordT*)calloc(3, sizeof(coordT));
                (p_minus->coord)[0] = center_max[0];
                (p_minus->coord)[1] = center_max[1];
                (p_minus->coord)[2] = center_max[2];
                p_minus->next = NULL;
                
                p_minus_temp = p_minus;
            }
            else
            {
                p_minus_temp->next = (pole_t*)malloc(sizeof(pole_t));
                p_minus_temp = p_minus_temp->next;
                
                p_minus_temp->coord = (coordT*)calloc(3, sizeof(coordT));
                (p_minus_temp->coord)[0] = center_max[0];
                (p_minus_temp->coord)[1] = center_max[1];
                (p_minus_temp->coord)[2] = center_max[2];
                p_minus_temp->next = NULL;
            }
        }
    }


    printf("SIZE = %d\n",vsize);
    // Count the number of positive poles
    int p_plus_size = 0;
    p_plus_temp = p_plus;
    while ( p_plus_temp!=NULL )
    {                                 
        p_plus_size++;
        p_plus_temp = p_plus_temp->next;
    }
    printf("POSITIVE POLES: %d\n",p_plus_size);
    
    // Count the number of negative poles
    int p_minus_size = 0;
    p_minus_temp = p_minus;
    while ( p_minus_temp!=NULL )
    {                                 
        p_minus_size++;
        p_minus_temp = p_minus_temp->next;
    }
    printf("NEGATIVE POLES: %d\n",p_minus_size);

    pole_t* P = p_plus;
    p_plus_temp = P;
    while ( (p_plus_temp->next)!=NULL )
        p_plus_temp = p_plus_temp->next;
    p_plus_temp->next = p_minus;

 printf("%f %f %f \n",p_plus->coord[0],p_plus->coord[1],p_plus->coord[2]);
    

    // Count the number of negative poles
    int P_size = 0;
    p_plus_temp = P;
    while ( p_plus_temp!=NULL )
    {                                 
        P_size++;
        p_plus_temp = p_plus_temp->next;
    }
    //printf("POLES: %d\n",P_size);  

    P_size = 0;
    p_plus_temp = P;
    while ( p_plus_temp!=NULL )
    {                                 
        P_size++;
        p_plus_temp = p_plus_temp->next;
    }
    printf("POLES: %d\n",P_size);   

//==(pole_consume_result*)============================================


//==(pole_clean)======================================================
    // Free allocated previously for "temp" purposes memory.
    free(pt);
    free(all_v);
    pt = NULL;
    all_v = NULL;
    
    // Request the "qhull" to free it's allocations.
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);
//==(pole_clean*)=====================================================

    return P;
}

/*-------------------------------------------------------------------*/
void Crust( void )
{
    // Use qhull to compute Delaunay triangulation of the points.
    // put results (vertices, edges, faces, and tetrahedra)
    // in the global structures above
    // see main.c for detail on how qhull does it
//==(init)============================================================
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
    pole_t* P = find_poles();// This would store the array of poles
    
    // Count the number of poles
    int P_size = 0;
    pole_t* P_temp = P;
    while ( P_temp!=NULL )
    {                                 
        P_size++;
        P_temp = P_temp->next;
    }

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
    pt = (coordT*)calloc((vsize+P_size) * 4, sizeof(coordT));  // This would be 
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

    P_temp = P;
    while ( P_temp!=NULL )
    {
        // Here "qhull" will obtain 4D points: 3 {x,y,z} 
        // coordinates + ...                                
        pt[id++] = P_temp->coord[0];
        pt[id++] = P_temp->coord[1];
        pt[id++] = P_temp->coord[2];
        //               ...4th for x^2+y^2+z^2.
        pt[id++] = P_temp->coord[0]*P_temp->coord[0] + 
                   P_temp->coord[1]*P_temp->coord[1] +
                   P_temp->coord[2]*P_temp->coord[2]; 
   
        P_temp = P_temp->next;
    }


//printf("%d %d\n",vsize, P_size);
//    for(id = 0; id < 4*(vsize+3); id++)
//        if(id % 4 == 3)
//            printf("%f \n---------------\n",pt[id]);
//        else
//            printf("%f ",pt[id]);

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
    qh_init_B (pt, vsize+P_size, 4, false);
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
        tVertex tri_points[4];
        // The next loop will just fill the tetrahedron added 
        // previously inserted into 
        // the G.V. by appropriate vertices obtained from qhull.
        vid = 0;
        FOREACHvertex_(facet->vertices)
        {
 
            // (Looping by all vertices of the "facet" and fill the 
            // tetrahedron).
            // (Looping is made by "vertex" variable).
            // ("qh_pointid(vertex->point)" is equal to the index 
            //   of the input point that is stored in the "vertex").
            if(qh_pointid(vertex->point) < vsize)
                tri_points[vid++] = all_v[qh_pointid(vertex->point)];

        }

        if(vid==3)
        {
            MakeFace( tri_points[0], tri_points[1], tri_points[2], NULL );
        }
        if(vid==4)
        {
            MakeFace( tri_points[0], tri_points[1], tri_points[2], NULL );
            MakeFace( tri_points[1], tri_points[2], tri_points[3], NULL );
            MakeFace( tri_points[2], tri_points[3], tri_points[0], NULL );
            MakeFace( tri_points[3], tri_points[0], tri_points[1], NULL );
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
