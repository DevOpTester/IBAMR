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

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/periodic_boundary.h>
#include <libmesh/gmv_io.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>
class IBFEModel{
    private:
        void parse_inputdb();
    public:
        // constructor/deconstructor
        IBFEModel::IBFEModel(int argc, const char * const * argv);
        virtual ~IBFEModel();

        // standard options set in input file
        const bool   dump_viz_data;
        const int    viz_dump_interval;
        const bool   uses_visit;
        const bool   uses_exodus;
        const string exodus_filename;
        const bool   uses_gmv;
        const string gmv_filename;
        const bool   dump_restart_data;
        const int    restart_dump_interval;
        const string restart_dump_dirname;
        const string restart_read_dirname;
        const int    restart_restore_num;
        const bool   dump_postproc_data;
        const int    postproc_data_dump_interval;
        const string postproc_data_dump_dirname;
        const bool   dump_timer_data;
        const int    timer_dump_interval;
        const double dx;
        const double ds;
        const string elem_type;
        const string solver_type;
        const bool   uses_inital_velocity;
        const bool   uses_initial_pressure;
        const bool   uses_forcing_function;
        const string data_dump_dirname;

        // objects associated with this model
        Pointer<AppInitializer> app_initializer;
        Pointer<Database> input_db;
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        Pointer<IBFEMethod> ib_method_ops;
        Pointer<IBHierarchyIntegrator> time_integrator;
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry;
        Pointer<StandardTagAndInitialize<NDIM> > error_detector;
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy;
        Pointer<BergerRigoutsos<NDIM> > box_generator;
        Pointer<LoadBalancer<NDIM> > load_balancer;
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm;
        Pointer<CartGridFunction> u_init;
        Pointer<CartGridFunction> p_init;
        Pointer<CartGridFunction> forcing_function;
        Pointer<VisItDataWriter<NDIM> > viz_data_writer;
        AutoPtr<ExodusII_IO> exodus_io;
        AutoPtr<GMVIO> gmv_io;
        EquationSystems* equation_systems;
        Mesh mesh;

        // public methods
        void log_input_database();
        void deallocate_init_objects();
        void init_hierarchy_config();
        void write_viz_files();
        void write_restart_files();
        void write_timer_data();
        void write_postproc_data();
        void output_data(const int iteration_num, const double loop_time);
};

// Constructor
IBFEModel::IBFEModel(int argc, const char * const * argv){
            app_initializer = new AppInitializer(argc, argv, "IB.log");
            parse_inputdb();
        }
IBFEModel::~IBFEModel(){

}

// public methods
void log_input_database(){

}

void deallocate_init_objects(){

}

void init_hierarchy_config(){

}

void write_viz_files(){

}

void write_restart_files(){

}

void write_timer_data(){

}

void write_postproc_data(){

}

void output_data(const int iteration_num, const double loop_time){
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

void set_viz_options(Pointer<AppInitializer> app_initializer){

}

void build_geometry(Mesh& mesh, Pointer<Database> input_db){

}

