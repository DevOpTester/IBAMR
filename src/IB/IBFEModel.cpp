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
#include "ibamr/IBFEModel.h"

// Set up application namespace declarations
#include "ibamr/app_namespaces.h"

// Constructor
IBFEModel(int argc, const char * const * argv)
{
    if(!model_exists){
        LibMeshInit init(argc, argv);
        app_initializer = new AppInitializer(argc, argv, "IB.log");
        input_db = app_initializer->getInputDatabase();
        parse_inputdb();
        TheModel = this;
    } else {
       // throw exception?
       // or could go with no public constructor and a factory method instead
    }
}


~IBFEModel()
{
   //TODO: delete TheModel.
}

// public methods
void log_input_database()
{

}

void deallocate_init_objects()
{

}

void init_hierarchy_config()
{

}

void write_viz_files()
{

}

void write_restart_files()
{

}

void write_timer_data()
{

}

void write_postproc_data()
{

}

void output_data(const int iteration_num, const double loop_time)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

    // Write Lagrangian data.
    file_name = data_dump_dirname + "/" + "fe_mesh.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    file_name += ".xda";
    mesh.write(file_name);
    file_name = data_dump_dirname + "/" + "fe_equation_systems.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    equation_systems->write(file_name, (EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA));
    return;
} // output_data

void set_viz_options(Pointer<AppInitializer> app_initializer)
{

}

void build_geometry(Mesh& mesh, Pointer<Database> input_db)
{
#if (NDIM == 2)
        MeshTools::Generation::build_square(mesh,
                                            static_cast<int>(ceil(0.1 / ds)),
                                            static_cast<int>(ceil(1.0 / ds)),
                                            0.95,
                                            1.05,
                                            0.0,
                                            1,
                                            Utility::string_to_enum<ElemType>(elem_type));
#endif
#if (NDIM == 3)
        mesh.read(input_db->getString("MESH_FILENAME"));
/*        MeshTools::Generation::build_cube(mesh,
                                          static_cast<int>(ceil(0.1 / ds)),
                                          static_cast<int>(ceil(1.0 / ds)),
                                          static_cast<int>(ceil(1.0 / ds)),
                                          0.95,
                                          1.05,
                                          0.0,
                                          1,
                                          0.0,
                                          1,
                                          Utility::string_to_enum<ElemType>(elem_type));
*/

        for (MeshBase::node_iterator it = mesh.nodes_begin();
             it != mesh.nodes_end(); ++it)
        {
            Node* n = *it;
            libMesh::Point& X = *n;
            X(0) += 0.9;
            X(1) += 0.5;
            X(2) += 0.5;
        }

#endif
        const MeshBase::const_element_iterator end_el = mesh.elements_end();
        for (MeshBase::const_element_iterator el = mesh.elements_begin(); el != end_el; ++el)
        {
            Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor(side);
                if (at_mesh_bdry)
                {
                    BoundaryInfo* boundary_info = mesh.boundary_info.get();
#if (NDIM == 2)
                    if (boundary_info->has_boundary_id(elem, side, 0) || boundary_info->has_boundary_id(elem, side, 2))
                    {
                        boundary_info->add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID);
                    }
#endif
#if (NDIM == 3)
                    if (!(boundary_info->has_boundary_id(elem, side, 2) ||
                          boundary_info->has_boundary_id(elem, side, 4)))
                    {
                        boundary_info->add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID);
                    }
#endif

}

void parse_inputdb()
{

        dump_viz_data         = app_initializer->dumpVizData();
        viz_dump_interval     = app_initializer->getVizDumpInterval();
        uses_visit            = dump_viz_data && app_initializer->getVisItDataWriter();
        uses_exodus           = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        exodus_filename       = app_initializer->getExodusIIFilename();
        uses_gmv              = dump_viz_data && !app_initializer->getGMVFilename().empty();
        gmv_filename          = app_initializer->getGMVFilename();
        dump_restart_data     = app_initializer->dumpRestartData();
        restart_dump_interval = app_initializer->getRestartDumpInterval();
        restart_dump_dirname  = app_initializer->getRestartDumpDirectory();
        restart_read_dirname  = app_initializer->getRestartReadDirectory();
        restart_restore_num   = app_initializer->getRestartRestoreNumber();
        dump_postproc_data    = app_initializer->dumpPostProcessingData();
        dx                    = input_db->getDouble("DX");
        ds                    = input_db->getDouble("MFAC") * dx;
        elem_type             = input_db->getString("ELEM_TYPE");

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
        mesh(init.comm(), NDIM);
}
