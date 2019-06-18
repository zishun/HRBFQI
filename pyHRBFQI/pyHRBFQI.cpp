// pyHRBFQI.cpp : pybind11 of HRBFQI.cpp
//

#ifdef WIN32
#include "stdafx.h"
#else
#include <cstring>
#include "string.h"
#endif
#include <iostream>
#include <string>
#include <time.h>
#include "../DataStructure/PointSet.h"
#include "../DataStructure/PolygonalMesh.h"
#include "../polygonizer/Polygonizer.h"
#include "../DataStructure/OctTree.h"
#include "../HRBFQI/HRBF.h"
#include "../HRBFQI/MeshCleaner.h"
#include "omp.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#define PI 3.14159265

using namespace std;
namespace py = pybind11;

std::tuple<py::array_t<double>,py::array_t<int>> hrbfqi(py::array_t<double> points, py::array_t<double> normals,
                           bool ifRescale, bool ifOutputDifferenceGradientNormal,
                           double supportSizeScale, double etaValue,
                           double gridEdgeLen, double confThreshold,
                           int componentSize, bool verbose)
{
    auto buf_pnts = points.request(), buf_nmls = normals.request();
    double *ptr_pnts = (double *) buf_pnts.ptr,
           *ptr_nmls = (double *) buf_nmls.ptr;

    if (buf_pnts.size != buf_nmls.size)
        throw std::runtime_error("Input shapes must match");

    PointSet *ps = new PointSet;
    float fScale = 1.0f;
    float fCenter[3] = {0.0f, 0.0f, 0.0f};

    int point_N = buf_pnts.shape[0];
    ps->setPointSize(point_N);

    int i;
    for(i=0; i<point_N; i++){
        ps->setPoint(i, ptr_pnts[i*3], ptr_pnts[i*3+1], ptr_pnts[i*3+2]);
    }
    for(i=0; i<point_N; i++){
        ps->setNormal(i, ptr_nmls[i*3], ptr_nmls[i*3+1], ptr_nmls[i*3+2]);
    }

    if(ifRescale)
    {
        fScale = ps->fitIntoBox(fCenter, 1.0f);
    }

    if (verbose)
    {
	cout<<ps->point_N<<" points input"<<endl;
    }

    //-------------------------------------------------------------------------
    //	Interpolation or Quasi-interpolation
    if (verbose)
    {
        cout<<"Fitting..."<<std::flush;
    }
    clock_t start = clock();

    HRBF *hrbf = new HRBF();
    hrbf->setPointSet(ps);

    //Estimation of recommended support size
    float T = 0.75*hrbf->getAveragedLeafSize();
//	cout<<"The initial support radius: "<< T <<endl;
	
    hrbf->support = supportSizeScale*T;
    T = hrbf->support;

    float normal_smooth;

    normal_smooth = etaValue;///(T*T);//100.0f/(T*T);
    hrbf->fit(hrbf->support, normal_smooth);
    if (verbose)
    {
        cout<<"Done!"<<endl;
    }
    if(!ifOutputDifferenceGradientNormal)
    {
        if (verbose)
        {
            cout<<"The time for fitting (including the octree construction): "<< (float)(clock()-start)/ CLOCKS_PER_SEC<<"s;"<<endl;
        }
        //---------------
        //	for debug
        //cout<<"--------------------------"<<endl;
        //cout<<"The final support radius: "<< hrbf->support<<endl;
        //int neighborNum = hrbf->getMaximalNeighborsInSupport(hrbf->support);	//	Compute a maximal neighbor number according to the given support size
        //cout<<"The maximal neighbor number: "<< neighborNum<< endl;
        //float rho_min = (5.0*neighborNum+sqrt(25.0*neighborNum*neighborNum+2240.0*(1+normal_smooth)))/(8.0*(1.0+normal_smooth));
        //cout<<"The parameter: "<< normal_smooth<< endl;//fit_dia->m_normal_smooth
        //cout<<"The minimal support radius: "<< rho_min<<endl;
    }
    //-------------------------------------------------------------------------
    //	Iso-surface extraction
    if (verbose)
    {
        cout<<"Surface extraction..."<<std::flush;
    }
    start = clock();
    float x_max, x_min, y_max, y_min, z_max, z_min;
    ps->getBound(x_min, x_max, y_min, y_max, z_min, z_max);
    float bound_X = x_max - x_min;
    float bound_Y = y_max - y_min;
    float bound_Z = z_max - z_min;
    float grid_size = gridEdgeLen;

    float space = grid_size;

    Polygonizer* poly = new Polygonizer;

    poly->spaceX = space;
    poly->spaceY = space;
    poly->spaceZ = space;

    poly->originX = x_min - 5*space;
    poly->originY = y_min - 5*space;
    poly->originZ = z_min - 5*space;

    poly->dimX = (int)((bound_X)/space) + 10;
    poly->dimY = (int)((bound_Y)/space) + 10;
    poly->dimZ = (int)((bound_Z)/space) + 10;

    poly->func = hrbf;

    PolygonalMesh *mesh;
    mesh = NULL;
    mesh = poly->dualContouring(0.001f, 0.01f);

    if (verbose)
    {
        cout<<"Done!"<<endl;
    }
    if(!ifOutputDifferenceGradientNormal)
    {
        if (verbose)
        {
            cout<<"The time for getting triangles: "<<((float)(clock()-start))/ CLOCKS_PER_SEC<<"s;"<<endl;
            cout<<"The extracted mesh has " << mesh->vertex_N<<" vertices and "<< mesh->face_N<<" triangles."<<endl;
        }
    }

    //---------------
    //	show the differences between gradients and normals
    if(ifOutputDifferenceGradientNormal)
    {
        double maxAngle, aveAngle;

        hrbf->differenceFuncGradientAndNormalPt(maxAngle, aveAngle);
        cout<<"The maximal angle between gradients and normals: "<<maxAngle<<endl;
        cout<<"The average angle between gradients and normals: "<< aveAngle << endl<< endl;

        // distances between points and the reconstructed mesh
        double maxDist, aveDist;
        mesh->distanceFromPts(maxDist, aveDist, ps);
        cout<<"The maximal distance between points and the reconstructed mesh: "<< maxDist<<endl;
        cout<<"The maximal distance between points and the reconstructed mesh: "<< aveDist << endl;
    }
    //-------------------------------------------------------------------------
    //	Mesh cleaning
    if(!ifOutputDifferenceGradientNormal)
    {
        if(confThreshold>0.0 || componentSize>0.0)
        {
            CMeshCleaner *meshCleaner = new CMeshCleaner;
            PolygonalMesh *newMesh;
            newMesh = NULL;

            if (verbose)
            {
                cout<<endl<<"********************************************************************"<<endl;
                cout<<"Mesh cleaning..."<<endl;
                cout<<"--------------------------"<<endl;
                cout<<"Removing low-confidence geometry (threshold "<< confThreshold << ")..." <<endl;
            }

            clock_t duration = 0;
            start = clock();
            newMesh = meshCleaner->removeLowConfidenceGeometry(mesh, hrbf, confThreshold);
            duration += clock()-start;

            if (verbose)
            {
                cout<<"Removed "<< mesh->vertex_N - newMesh->vertex_N << " vertices and " << mesh->face_N - newMesh->face_N<< " triangles."<<endl;
            }

            delete mesh;
            //delete ps;
            //delete hrbf;
            //delete poly;

            if (componentSize > 0)
            {
                if (verbose)
                {
                    cout<<"--------------------------"<<endl;
                    cout<<"Removing isolated components with < " << componentSize << " vertices..." <<endl;
                }

                mesh = newMesh;
                start = clock();
                newMesh = meshCleaner->removeSmallIsolatedComponent(mesh, componentSize);
                duration += clock()-start;

                if (verbose)
                {
                    cout<<"Removed "<< mesh->vertex_N - newMesh->vertex_N << " vertices and " << mesh->face_N - newMesh->face_N<< " triangles."<<endl;
                }

                delete mesh;
            }
            if (verbose)
            {
                cout<<"--------------------------"<<endl;
                cout<<"The time for mesh cleaning:"<<((float)(duration))/ CLOCKS_PER_SEC<<"s."<<endl;
            }
            mesh = newMesh;
            delete meshCleaner;
        }

        if(mesh == NULL || mesh->face_N == 0)
        {
            throw std::runtime_error("Fail to reconstruct mesh!");
        }

        float (*vertex)[3] = mesh->vertex;
        int (**face) = mesh->face;
        int vertex_N = mesh->vertex_N;
        int face_N = mesh->face_N;

        py::array_t<double> mesh_verts = py::array_t<double>(vertex_N*3);
        py::array_t<int>    mesh_faces = py::array_t<int>(face_N*3);
        mesh_verts.resize({vertex_N, 3});
        mesh_faces.resize({face_N,3});
        py::buffer_info buf_mesh_v = mesh_verts.request();
        py::buffer_info buf_mesh_f = mesh_faces.request();
        double* ptr_mesh_v = (double*)buf_mesh_v.ptr;
        int*    ptr_mesh_f =    (int*)buf_mesh_f.ptr;

        int i;
        for(i=0; i<vertex_N; i++){
            ptr_mesh_v[i*3] = vertex[i][0]/fScale+fCenter[0];
            ptr_mesh_v[i*3+1] = vertex[i][1]/fScale+fCenter[1];
            ptr_mesh_v[i*3+2] = vertex[i][2]/fScale+fCenter[2];
        }

        for(i=0; i<face_N; i++){
            ptr_mesh_f[i*3] = face[i][0];
            ptr_mesh_f[i*3+1] = face[i][1];
            ptr_mesh_f[i*3+2] = face[i][2];
        }

        if (verbose)
        {
            cout<<mesh->vertex_N <<" vertices and " << mesh->face_N << " triangles output."<<endl;
        }

        delete mesh;
        delete ps;
        delete hrbf;
        delete poly;

        auto mesh_tuple = std::make_tuple(mesh_verts, mesh_faces);
        return mesh_tuple;
    }
}

PYBIND11_MODULE(pyHRBFQI, m) {
    m.doc() = "HRBF-based surface reconstruction";

    m.def("hrbfqi", &hrbfqi, "HRBF-based surface reconstruction. <if_rescale> <if_output_difference_gradient_normal> <support_size_scale>"
              " <parameter_eta> <isosurface_extract_stepsize> <confidence_threshold> <component_size> <if_verbose_output>");
}
