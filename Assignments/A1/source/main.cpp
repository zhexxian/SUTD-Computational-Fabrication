//Computational Fabrication Assignment #1
// By David Levin 2014
#include <iostream>
#include <vector>
#include "../include/CompFab.h"
#include "../include/Mesh.h"

#define PI 3.14159265

//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle, 
     * 0 otherwise */
	
	// calculate normal of triangle
	CompFab::Vec3 v1v2 = triangle.m_v2 - triangle.m_v1;
	CompFab::Vec3 v1v3 = triangle.m_v3 - triangle.m_v1;
	CompFab::Vec3 normal = CompFab::operator%(v1v2, v1v3);
    

	// calculate determinant
	CompFab::Vec3 p = CompFab::operator%(ray.m_direction, v1v3);
	float det = CompFab::operator*(v1v2, p);

	// check if ray and plane are parallel
	if(det > -EPSILON && det < EPSILON) {
		return 0;
	}
	float inv_det = 1.0f / det;


	//calculate distance from V1 to ray origin
	CompFab::Vec3 T = ray.m_origin - triangle.m_v1;

	// compute u
	float u = CompFab::operator*(T,p) * inv_det;

	//The intersection lies outside of the triangle
	if(u < 0.f || u > 1.f) {
		return 0;
	}

	//Prepare to test v parameter
	CompFab::Vec3 Q = CompFab::operator%(T, v1v2);

	//Calculate V parameter and test bound
	float v = CompFab::operator*(ray.m_direction,Q) * inv_det;
 
	//The intersection lies outside of the triangle
	if(v < 0.f || u + v  > 1.f) {
		return 0;
	}

    // compute t (ray parameter)
	float t = CompFab::operator*(v1v3, Q) * inv_det;
	if(t > EPSILON) { //ray intersection
		 return 1;
	}

	return 0;
}

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    
    unsigned int numHits = 0;
    
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir, 
     * from voxel center voxelPos intersects the surface */

	CompFab::Ray voxelRay = CompFab::Ray(voxelPos, dir);

	for (int i=0; i<g_triangleList.size(); i++){
		numHits += rayTriangleIntersection(voxelRay, g_triangleList[i]);
	}
    return numHits;
}

bool loadMesh(char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    delete tempMesh;
    
    return true;
   
}

void saveVoxelsToObj(const char * outfile)
{
 
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}


int main(int argc, char **argv)
{

    //unsigned int dim = 32; //dimension of voxel grid (e.g. 32x32x32)
	unsigned int dim = 64;


    //Load OBJ
    if(argc < 3)
    {
        std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename \n";
        return 0;
    }
    
    std::cout<<"Load Mesh : "<<argv[1]<<"\n";
    loadMesh(argv[1], dim);
    
    

    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction(1.0,0.0,0.0);

	

	//EXTRA CREDIT: cast ray in different directions
	std::vector<CompFab::Vec3> newDirections;
	float theta_r, XTheta, YTheta, phi_r, XPhi, YPhi, X, Y, Z;
	int stepSize = 60;
	for(int theta=0; theta<360; theta=theta+stepSize){
		for(int phi=0; phi<360; phi=phi+stepSize){
			theta_r = theta * PI / 180.0;
			phi_r = phi * PI / 180.0;
			XTheta = cos(theta_r);
			YTheta = sin(theta_r);
			XPhi = cos(phi_r);
			YPhi = sin(phi_r);
			X = XTheta * XPhi;
			Y = YTheta * XPhi;
			Z = YPhi;
			newDirections.push_back(CompFab::Vec3(X, Y, Z));
		}
	}
				
	
    /********* ASSIGNMENT *********/
    /* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in g_triangleList */

    //Cast ray, check if voxel is inside or outside
    //even number of surface intersections = outside (OUT then IN then OUT)
    //odd number = inside (IN then OUT)

	for(unsigned int i=0; i<g_voxelGrid->m_dimX; i++){
		for(unsigned int j=0; j<g_voxelGrid->m_dimY; j++){
			for(unsigned int k=0; k<g_voxelGrid->m_dimZ; k++){
				voxelPos = (g_voxelGrid->m_lowerLeft)
					+CompFab::Vec3(
					((double)i)* (g_voxelGrid->m_spacing), 
					((double)j)* (g_voxelGrid->m_spacing), 
					((double)k)* (g_voxelGrid->m_spacing));
				// Normal, cast ray in a single direction (without extra credit)
				/*
				if(numSurfaceIntersections(voxelPos, direction) % 2 == 0){
					g_voxelGrid->isInside(i, j, k) = false;
				}
				else {
					g_voxelGrid->isInside(i, j, k) = true;
				}
				*/

				// Extra credit
				
				g_voxelGrid->isInside(i, j, k) = false;
				for(int index=0; index<newDirections.size(); index++){
					if (g_voxelGrid->isInside(i, j, k)){
						break;
					}
					else if(numSurfaceIntersections(voxelPos, newDirections[index]) % 2 !=0) {
						g_voxelGrid->isInside(i, j, k) = true;
					}
				}
				
				

			}
		}
	} 

    //Write out voxel data as obj
    saveVoxelsToObj(argv[2]);
    
    delete g_voxelGrid;
}