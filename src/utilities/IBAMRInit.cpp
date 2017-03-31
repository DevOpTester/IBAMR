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

Pointer<INSHierarchyIntegrator>          IBAMRInit::navier_stokes_integrator = NULL;
Pointer<IBFEMethod>                      IBAMRInit::ib_method_ops = NULL;
SAMRAI::tbox::Pointer< INSHierarchyIntegrator > IBAMRInit::ins_hier_integrator = NULL;
Pointer<CartesianGridGeometry<NDIM> >    IBAMRInit::grid_geometry = NULL;
Pointer<IBHierarchyIntegrator>           IBAMRInit::time_integrator = NULL;
Pointer<PatchHierarchy<NDIM> >           IBAMRInit::patch_hierarchy = NULL;
Pointer<StandardTagAndInitialize<NDIM> > IBAMRInit::error_detector = NULL;
Pointer<BergerRigoutsos<NDIM> >          IBAMRInit::box_generator = NULL;
Pointer<LoadBalancer<NDIM> >             IBAMRInit::load_balancer = NULL;
Pointer<GriddingAlgorithm<NDIM> >        IBAMRInit::gridding_algorithm = NULL;

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
IBAMRInit::IBAMRInit(Mesh * m)
{
    mesh = m;
} // IBAMRInit

IBAMRInit::~IBAMRInit()
{
    deallocateAppInitializer();
    init_exists = false;
    argc = 0;
    argv = NULL;

    navier_stokes_integrator = NULL;
    ib_method_ops = NULL;
    ins_hier_integrator = NULL;
    grid_geometry = NULL;
    time_integrator = NULL;
    patch_hierarchy = NULL;
    error_detector = NULL;
    box_generator = NULL;
    load_balancer = NULL;
    gridding_algorithm = NULL;
}

void
IBAMRInit::deallocateAppInitializer()
{
    app_initializer.setNull();
}

//getters to access private and protected objects

Pointer<AppInitializer>
IBAMRInit::getAppInitializer()
{
    if (!app_initializer){
        TBOX_ERROR("App initializer not located!" <<"\n");
    }
    return app_initializer;
}

Pointer<Database>
IBAMRInit::getInputDB()
{
    if (!input_db){
        TBOX_ERROR("No input database located!" <<"\n");
    }
    return input_db;
}

Pointer<INSHierarchyIntegrator>
IBAMRInit::getIntegrator()
{
    if (!navier_stokes_integrator){
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
    }
    return navier_stokes_integrator;
}

Pointer<IBFEMethod>
IBAMRInit::getIBFEMethod()
{
    if (!ib_method_ops){
        if (input_db->keyExists("IBFEMethod")){
        ib_method_ops = new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           mesh,
                           max_levels,
                           restart_enabled,
                           restart_read_dirname,
                           restart_restore_num);
        }
        else
        {
            TBOX_ERROR("You made a call to getIBFEMethod"
                <<" but no IBFEMethod database was found in the input file"<< "\n");
        }
    }
    return ib_method_ops;

}



Pointer<IBHierarchyIntegrator>
IBAMRInit::getExplicitTimeIntegrator( SAMRAI::tbox::Pointer< IBStrategy > ib_method_ops,
                              SAMRAI::tbox::Pointer< INSHierarchyIntegrator > ins_hier_integrator )
{
    if (!ib_method_ops){
        ib_method_ops = getIBFEMethod();
    }
    if (!ins_hier_integrator){
        ins_hier_integrator = getIntegrator();
    }
    if (!time_integrator){
        if (input_db->keyExists("IBHierarchyIntegrator")){
            time_integrator = new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                    app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                    ib_method_ops,
                    ins_hier_integrator, restart_enabled);
        }
        else
        {
            TBOX_ERROR("You made a call to getExplicitTimeIntegrator"
                <<" but no IBHierarchyIntegrator database was found in the input file"<< "\n");
        }
    }
    return time_integrator;
}

Pointer<CartesianGridGeometry<NDIM> >
IBAMRInit::getCartesianGridGeometry()
{
    if (!grid_geometry){
        if (input_db->keyExists("CartesianGeometry")){
            grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        }
        else
        {
            TBOX_ERROR("You made a call to getCartesianGridGeometry"
                <<" but no CartesianGeometry database was found in the input file"<< "\n");
        }
    }
    return grid_geometry;
}


Pointer<PatchHierarchy<NDIM> >
IBAMRInit::getPatchHierarchy(Pointer<CartesianGridGeometry<NDIM> > grid_geometry)
{
    if (!patch_hierarchy){
        if (!grid_geometry){
            grid_geometry = getCartesianGridGeometry();
        }
        patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    }
    return patch_hierarchy;
}

Pointer<StandardTagAndInitialize<NDIM> >
IBAMRInit::getErrorDetector( Pointer<IBHierarchyIntegrator> time_integrator )
{
    if (!error_detector){
        if (!time_integrator){
            time_integrator = getExplicitTimeIntegrator();
        }
        if (input_db->keyExists("StandardTagAndInitialize")){
            error_detector = new StandardTagAndInitialize<NDIM>("StandardTagAndInitialzie",
                                                time_integrator,
                                                getAppInitializer()->getComponentDatabase("StandardTagAndInitialize"));
        }
        else
        {
            TBOX_ERROR("You made a call to getErrorDetector"
                <<" but no StandardTagAndInitialize database was found in the input file"<< "\n");
        }
    }
    return error_detector;
}

Pointer<BergerRigoutsos<NDIM> >
IBAMRInit::getBoxGenerator()
{
    if (!box_generator){
        box_generator = new BergerRigoutsos<NDIM>();
    }
    return box_generator;
}

Pointer<LoadBalancer<NDIM> >
IBAMRInit::getLoadBalancer()
{
    if (!load_balancer){
        if (input_db->keyExists("LoadBalancer")){

            load_balancer = new LoadBalancer<NDIM>("LoadBalancer", getAppInitializer()->getComponentDatabase("LoadBalancer"));
        }
        else
        {
            TBOX_ERROR("You made a call to getLoadBalancer"
                <<" but no LoadBalancer database was found in the input file"<< "\n");
        }
    }
    return load_balancer;
}


Pointer<GriddingAlgorithm<NDIM> >
IBAMRInit::getGriddingAlgorithm( Pointer<StandardTagAndInitialize<NDIM> > error_detector,
                                 Pointer<BergerRigoutsos<NDIM> > box_generator,
                                 Pointer<LoadBalancer<NDIM> > load_balancer)
{
    if (!gridding_algorithm){
        if (!error_detector){
            error_detector = getErrorDetector();
        }
        if (!box_generator){
            box_generator = getBoxGenerator();
        }
        if (!load_balancer){
            load_balancer = getLoadBalancer();
        }
        if (input_db->keyExists("GriddingAlgorithm")){
            gridding_algorithm = new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        getAppInitializer()->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        }
        else
        {
            TBOX_ERROR("You made a call to getGriddingAlgorithm"
                <<" but no GriddingAlgorithm database was found in the input file"<< "\n");
        }
    }
    return gridding_algorithm;
}

libMesh::Order
IBAMRInit::getPK1DevOrder(string DEFAULT)
{
    return Utility::string_to_enum<libMesh::Order>(getInputDB()->getStringWithDefault("PK1_DEV_QUAD_ORDER", DEFAULT));
}

libMesh::Order
IBAMRInit::getPK1DilOrder(string DEFAULT)
{

    return Utility::string_to_enum<libMesh::Order>(getInputDB()->getStringWithDefault("PK1_DIL_QUAD_ORDER", DEFAULT));
}

void
IBAMRInit::log_start(int iteration_num, double loop_time)
{
    pout << "\n";
    pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    pout << "At beginning of timestep # " << iteration_num << "\n";
    pout << "Simulation time is " << loop_time << "\n";
}

void
IBAMRInit::log_end(int iteration_num, double loop_time)
{
    pout << "\n";
    pout << "At end       of timestep # " << iteration_num << "\n";
    pout << "Simulation time is " << loop_time << "\n";
    pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    pout << "\n";
}

void
dump_data(int iteration_num, double loop_time)
{
/*
    IBAMRInit * self = ibamr_init;
    const bool last_step = !time_integrator->stepsRemaining();
    if (self->dump_viz_data && (iteration_num % self->viz_dump_interval == 0 || last_step))
    {
        pout << "\nWriting visualization files...\n\n";
        if (self->uses_visit)
        {
            time_integrator->setupPlotData();
            self->visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }
        if (uses_exodus)
        {
            exodus_io->write_timestep(exodus_filename,
                                        *equation_systems,
                                        iteration_num / viz_dump_interval + 1,
                                        loop_time);
        }
        if (uses_gmv)
        {
            std::ostringstream file_name;
            file_name << gmv_filename + "_" << std::setw(6) << std::setfill('0') << std::right << iteration_num;
            gmv_io->write_equation_systems(file_name.str() + ".gmv", *equation_systems);
        }
    }
    if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
    {
        pout << "\nWriting restart files...\n\n";
        RestartManager::getManager()->writeRestartFile(ibamr_init.restart_dump_dirname, iteration_num);
        ib_method_ops->writeFEDataToRestartFile(ibamr_init.restart_dump_dirname, iteration_num);
    }
    if (ibamr_init.dump_timer_data && (iteration_num % ibamr_init.timer_dump_interval == 0 || last_step))
    {
        pout << "\nWriting timer data...\n\n";
        TimerManager::getManager()->print(plog);
    }
    if (ibamr_init.dump_postproc_data && (iteration_num % ibamr_init.postproc_data_dump_interval == 0 || last_step))
    {
        pout << "\nWriting state data...\n\n";
        output_data(patch_hierarchy,
                    navier_stokes_integrator,
                    mesh,
                    equation_systems,
                    iteration_num,
                    loop_time,
                    ibamr_init.postproc_data_dump_dirname);
    }
*/
}

void
IBAMRInit::build_square(double xmin, double xmax, double ymin, double ymax )
{
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
