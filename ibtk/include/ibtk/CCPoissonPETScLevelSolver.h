// Filename: CCPoissonPETScLevelSolver.h
// Created on 30 Aug 2010 by Boyce Griffith
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

#ifndef included_IBTK_CCPoissonPETScLevelSolver
#define included_IBTK_CCPoissonPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "CellVariable.h"
#include "IntVector.h"
#include "RefineSchedule.h"
#include "VariableContext.h"
#include "ibtk/PETScLevelSolver.h"
#include "ibtk/PoissonSolver.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CCPoissonPETScLevelSolver is a concrete PETScLevelSolver for
 * solving elliptic equations of the form \f$ \mbox{$L u$} = \mbox{$(C I +
 * \nabla \cdot D \nabla) u$} = f \f$ on a \em single SAMRAI::hier::PatchLevel.
 *
 * This solver class uses the PETSc library to solve linear equations of the
 * form \f$ (C I + \nabla \cdot D \nabla ) u = f \f$, where \f$C\f$ and \f$D\f$
 * are scalars, and \f$u\f$ and \f$f\f$ are cell-centered arrays.  The
 * discretization is second-order accurate.
 *
 * Robin boundary conditions may be specified through the interface class
 * SAMRAI::solv::RobinBcCoefStrategy.
 *
 * The user must perform the following steps to use class
 * CCPoissonPETScLevelSolver:
 *
 * -# Create a CCPoissonPETScLevelSolver object.
 * -# Set the problem specification via setPoissonSpecifications(),
 *    setPhysicalBcCoef(), and setHomogeneousBc().
 * -# Initialize CCPoissonPETScLevelSolver object using the function
 *    initializeSolverState().
 * -# Solve the linear system using the member function solveSystem(), passing
 *    in SAMRAI::solv::SAMRAIVectorReal objects corresponding to \f$u\f$ and
 *    \f$f\f$.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 options_prefix = ""           // see setOptionsPrefix()
 ksp_type = "gmres"            // see setKSPType()
 initial_guess_nonzero = TRUE  // see setInitialGuessNonzero()
 rel_residual_tol = 1.0e-5     // see setRelativeTolerance()
 abs_residual_tol = 1.0e-50    // see setAbsoluteTolerance()
 max_iterations = 10000        // see setMaxIterations()
 enable_logging = FALSE        // see setLoggingEnabled()
 \endverbatim
 *
 * PETSc is developed at the Argonne National Laboratory Mathematics and
 * Computer Science Division.  For more information about \em PETSc, see <A
 * HREF="http://www.mcs.anl.gov/petsc/petsc-as">http://www.mcs.anl.gov/petsc/petsc-as</A>.
 */
class CCPoissonPETScLevelSolver : public PETScLevelSolver, public PoissonSolver
{
public:
    /*!
     * \brief Constructor.
     */
    CCPoissonPETScLevelSolver(const std::string& object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                              const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~CCPoissonPETScLevelSolver();

    /*!
     * \brief Static function to construct a CCPoissonPETScLevelSolver.
     */
    static SAMRAI::tbox::Pointer<PoissonSolver> allocate_solver(const std::string& object_name,
                                                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                                const std::string& default_options_prefix)
    {
        return new CCPoissonPETScLevelSolver(object_name, input_db, default_options_prefix);
    } // allocate_solver

protected:
    /*!
     * \brief Generate IS/subdomains for Schwartz type preconditioners.
     */
    void generateASMSubdomains(std::vector<std::set<int> >& overlap_is, std::vector<std::set<int> >& nonoverlap_is);
    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     */
    void initializeSolverStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                          const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverStateSpecialized().
     */
    void deallocateSolverStateSpecialized();

    /*!
     * \brief Copy a generic vector to the PETSc representation.
     */
    void copyToPETScVec(Vec& petsc_x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x);

    /*!
     * \brief Copy a generic vector from the PETSc representation.
     */
    void copyFromPETScVec(Vec& petsc_x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x);

    /*!
     * \brief Copy solution and right-hand-side data to the PETSc
     * representation, including any modifications to account for boundary
     * conditions.
     */
    void setupKSPVecs(Vec& petsc_x,
                      Vec& petsc_b,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CCPoissonPETScLevelSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CCPoissonPETScLevelSolver(const CCPoissonPETScLevelSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CCPoissonPETScLevelSolver& operator=(const CCPoissonPETScLevelSolver& that);

    /*!
     * \name PETSc objects.
     */
    //\{
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    std::vector<int> d_num_dofs_per_proc;
    int d_dof_index_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, int> > d_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > d_data_synch_sched, d_ghost_fill_sched;
    //\}
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CCPoissonPETScLevelSolver
