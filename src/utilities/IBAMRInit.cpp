// Filename: IBFEModel.cpp
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

// Config files
#include "IBAMR_config.h"
#include "IBTK_config.h"
#include "SAMRAI_config.h"

// Headers for basic PETSc functions
#include "petscsys.h"

// Headers for basic SAMRAI objects
#include "BergerRigoutsos.h"
#include "CartesianGridGeometry.h"
#include "LoadBalancer.h"
#include "StandardTagAndInitialize.h"

// Headers for basic libMesh objects
#include "libmesh/boundary_info.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/gmv_io.h"

// Headers for application-specific algorithm/data structure objects
#include "boost/multi_array.hpp"
#include "ibamr/IBExplicitHierarchyIntegrator.h"
#include "ibamr/IBFEMethod.h"
#include "ibamr/INSCollocatedHierarchyIntegrator.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibtk/AppInitializer.h"
#include "ibtk/libmesh_utilities.h"
#include "ibtk/muParserCartGridFunction.h"
#include "ibtk/muParserRobinBcCoefs.h"

// include class declaration
#include "ibamr/IBAMRInit.h"

// Set up application namespace declarations
#include "ibamr/app_namespaces.h"

// define static class members
bool IBAMRInit::init_exists = false;
IBAMRInit * IBAMRInit::ibamr_init = NULL;
int IBAMRInit::argc = 0;
char** IBAMRInit::argv = NULL;

// Factory method
IBAMRInit IBAMRInit::getInstance(int num_args, char** args, Mesh * m)
{
    if (! init_exists)
    {
        IBAMRInit::argc = num_args;
        IBAMRInit::argv = args;
        ibamr_init = new IBAMRInit(m);
        init_exists = true;
        ibamr_init->parse_inputdb();
        return (*ibamr_init);
    }
    return (*ibamr_init);
}

// Constructor
IBAMRInit::IBAMRInit(Mesh * m){
    mesh = m;
} // IBAMRInit

IBAMRInit::~IBAMRInit()
{
    deallocateAppInitializer();
    init_exists = false;
}

void
IBAMRInit::deallocateAppInitializer()
{
    app_initializer.setNull();
}

//getters to access private and protected objects

Pointer<AppInitializer>
IBAMRInit::getAppInitializer(){
    return app_initializer;
}

Pointer<Database>
IBAMRInit::getInputDB(){
    return input_db;
}

Pointer<INSHierarchyIntegrator>
IBAMRInit::getIntegrator(){
    if (solver_type == "STAGGERED")
    {
        navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    }
    else if (solver_type == "COLLOCATED")
    {
        navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
            "INSCollocatedHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
    }
    else
    {
        TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                    << "Valid for options for SOLVER_TYPE are: COLLOCATED, STAGGERED");
    }
    return navier_stokes_integrator;
}

Pointer<IBFEMethod>
IBAMRInit::getIBFEMethod()
{
    return new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           mesh,
                           max_levels,
                           restart_enabled,
                           restart_read_dirname,
                           restart_restore_num);
}


Pointer<IBHierarchyIntegrator>
IBAMRInit::getExplicitTimeIntegrator( SAMRAI::tbox::Pointer< IBStrategy > ib_method_ops,
                              SAMRAI::tbox::Pointer< INSHierarchyIntegrator > ins_hier_integrator )
{
    return new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                    app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                    ib_method_ops,
                    ins_hier_integrator, restart_enabled);
}

Pointer<CartesianGridGeometry<NDIM> >
IBAMRInit::getCartesianGridGeometry()
{
        return new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
}


Pointer<PatchHierarchy<NDIM> >
IBAMRInit::getPatchHierarchy(Pointer<CartesianGridGeometry<NDIM> > grid_geometry){
    return new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
}

void
IBAMRInit::build_square(double xmin, double xmax, double ymin, double ymax ){
   //build_square (UnstructuredMesh &mesh, const unsigned int nx, const unsigned int ny,
   //const Real xmin=0., const Real xmax=1., const Real ymin=0., const Real ymax=1.,
   //        const ElemType type=INVALID_ELEM, const bool gauss_lobatto_grid=false)
    MeshTools::Generation::build_square((*mesh) , static_cast<int>(ceil(0.1) / ds),
                                        static_cast<int>(ceil(1.0 / ds)),
                                        xmin, xmax, ymin, ymax, Utility::string_to_enum<ElemType>(elem_type));

}

void
IBAMRInit::translate_mesh(double xdisplacement, double ydisplacement )
{
    // 2D translation
        for (MeshBase::node_iterator it = mesh->nodes_begin();
             it != mesh->nodes_end(); ++it)
        {
            Node* n = *it;
            libMesh::Point& X = *n;
            X(0) += xdisplacement;
            X(1) += ydisplacement;
        }
}

void
IBAMRInit::translate_mesh(double xdisplacement, double ydisplacement, double zdisplacement)
{
    //3D translation
        for (MeshBase::node_iterator it = mesh->nodes_begin();
             it != mesh->nodes_end(); ++it)
        {
            Node* n = *it;
            libMesh::Point& X = *n;
            X(0) += xdisplacement;
            X(1) += ydisplacement;
            X(2) += zdisplacement;
        }
}




// private methods
void
IBAMRInit::parse_inputdb()
{
    app_initializer       = new AppInitializer(argc, argv, "IB.log");
    input_db              = app_initializer->getInputDatabase();
    dump_viz_data         = app_initializer->dumpVizData();
    viz_dump_interval     = app_initializer->getVizDumpInterval();

    uses_visit            = dump_viz_data && app_initializer->getVisItDataWriter();
    uses_exodus           = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
    uses_gmv              = dump_viz_data && !app_initializer->getGMVFilename().empty();
    uses_forcing_function = input_db->keyExists("ForcingFunction");

    gmv_filename          = app_initializer->getGMVFilename();
    exodus_filename       = app_initializer->getExodusIIFilename();
    dump_postproc_data    = app_initializer->dumpPostProcessingData();
    dx                    = input_db->getDouble("DX");
    ds                    = input_db->getDouble("MFAC") * dx;
    elem_type             = input_db->getString("ELEM_TYPE");
    max_levels            = app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels");

    if (input_db->keyExists("enable_restart")){
        if ((input_db->getString("enable_restart")).compare("TRUE") == 0){
            restart_enabled = true;
        }
    }
    if (restart_enabled){
        dump_restart_data     = app_initializer->dumpRestartData();
        restart_dump_interval = app_initializer->getRestartDumpInterval();
        restart_dump_dirname  = app_initializer->getRestartDumpDirectory();
        restart_read_dirname  = app_initializer->getRestartReadDirectory();
        restart_restore_num   = app_initializer->getRestartRestoreNumber();
    }else {
        restart_dump_dirname  = "";
        restart_read_dirname  = "";
    }
    if (input_db->keyExists("mesh_option")){
        mesh_option       = input_db->getString("MESH_OPTION");
    } else {
        mesh_option = "";
    }

    if (input_db->keyExists("input_mesh_file")){
        input_mesh_filename = input_db->getString("INPUT_MESH_FILE");
    }

    postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
    postproc_data_dump_dirname  = app_initializer->getPostProcessingDataDumpDirectory();
    if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
    {
        Utilities::recursiveMkdir(postproc_data_dump_dirname);
    }

    dump_timer_data     = app_initializer->dumpTimerData();
    timer_dump_interval = app_initializer->getTimerDumpInterval();

    // Create a simple FE mesh with Dirichlet boundary conditions.
    // Note that boundary condition data must be registered with each FE
    // system before calling IBFEMethod::initializeFEData().
    solver_type         = app_initializer->getComponentDatabase("Main")->getString("solver_type");
}

void
IBAMRInit::registerVelocityInitialConditions(Pointer<CartesianGridGeometry<NDIM> > grid_geometry)
{
    if (input_db->keyExists("VelocityInitialConditions"))
    {
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
    }
}


void
IBAMRInit::registerPressureInitialConditions(Pointer<CartesianGridGeometry<NDIM> > grid_geometry)
{
    if (input_db->keyExists("PressureInitialConditions"))
    {
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);
    }
}
