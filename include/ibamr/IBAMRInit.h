// Filename: IBFEModel.h
// Created on 5 Oct 2017 by Elijah DeLee
//
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

#ifndef included_IBAMR_IBAMRInit
#define included_IBAMR_IBAMRInit

/////////////////////////////// INCLUDES /////////////////////////////////////
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

// Set up application namespace declarations
#include "ibamr/app_namespaces.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

class IBAMRInit{
    private:
        static bool init_exists;
        static IBAMRInit * ibamr_init;  // singleton
         // private constructor, use factory method to get instance
        IBAMRInit();
        // internal methods
    protected:
        Pointer<IBFEMethod>             ibfe_method;
        Pointer<AppInitializer>         app_initializer;
        Pointer<Database>               input_db;
        static LibMeshInit *            libmesh_init;
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        Mesh *                          mesh;
        static int                      argc;
        static char**                   argv;


    public:
        static  IBAMRInit getInstance(int argc, char** argv, LibMeshInit * lm_init);
        virtual ~IBAMRInit();
        void parse_inputdb();
        void    deallocateAppInitializer();
        // standard options set in input file
        bool   dump_viz_data;
        int    viz_dump_interval;
        bool   uses_visit;
        bool   uses_exodus;
        string exodus_filename;
        bool   uses_gmv;
        string gmv_filename;
        bool   restart_enabled;
        bool   dump_restart_data;
        int    restart_dump_interval;
        string restart_dump_dirname;
        string restart_read_dirname;
        string bc_coefs_name;
        string bc_coefs_db_name;
        int    restart_restore_num;
        bool   dump_postproc_data;
        int    postproc_data_dump_interval;
        string postproc_data_dump_dirname;
        bool   dump_timer_data;
        int    timer_dump_interval;
        double dx;
        double ds;
        string elem_type;
        string solver_type;
        bool   uses_inital_velocity;
        bool   uses_initial_pressure;
        bool   uses_forcing_function;
        int    max_levels;
        string data_dump_dirname;
        string mesh_option;       // options include: cube, square, cylinder, sphere
        string input_mesh_filename;

        // methods to modify objects owned by init object
        void registerVelocityInitialConditions(Pointer<CartesianGridGeometry<NDIM> > grid_geometry);
        void registerPressureInitialConditions(Pointer<CartesianGridGeometry<NDIM> > grid_geometry);

        // methods to retrieve objects associated, throw useful
        // error messages when unable to retrieve
        Pointer<AppInitializer>                  getAppInitializer();
        Pointer<Database>                        getInputDB();
        Mesh                                     getMesh();
        Pointer<INSHierarchyIntegrator>          getIntegrator();
};
#endif //#ifndef included_IBAMR_IBAMRInit
