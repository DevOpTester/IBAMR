// Filename: CartGridFunctionSet.cpp
// Created on 05 Sep 2012 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ostream>
#include <string>
#include <vector>

#include "SAMRAI/hier/BasePatchLevel.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/math/HierarchyDataOpsManager.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/MultiblockDataTranslator.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/math/PatchEdgeDataBasicOps.h"
#include "SAMRAI/math/PatchFaceDataBasicOps.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/math/PatchNodeDataBasicOps.h"
#include "SAMRAI/math/PatchSideDataBasicOps.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartGridFunctionSet.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartGridFunctionSet::CartGridFunctionSet(const std::string& object_name) : CartGridFunction(object_name), d_fcns()
{
    // intentionally blank
    return;
} // CartGridFunctionSet

CartGridFunctionSet::~CartGridFunctionSet()
{
    // intentionally blank
    return;
} // ~CartGridFunctionSet

void CartGridFunctionSet::addFunction(Pointer<CartGridFunction> fcn)
{
    TBOX_ASSERT(fcn);
    d_fcns.push_back(fcn);
    return;
} // addFunction

bool CartGridFunctionSet::isTimeDependent() const
{
    for (unsigned int k = 0; k < d_fcns.size(); ++k)
    {
        if (d_fcns[k]->isTimeDependent()) return true;
    }
    return false;
} // isTimeDependent

void CartGridFunctionSet::setDataOnPatchHierarchy(const int data_idx,
                                                  Pointer<Variable > var,
                                                  Pointer<PatchHierarchy > hierarchy,
                                                  const double data_time,
                                                  const bool initial_time,
                                                  const int coarsest_ln_in,
                                                  const int finest_ln_in)
{
    TBOX_ASSERT(hierarchy);
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    const int cloned_data_idx = var_db->registerClonedPatchDataIndex(var, data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(cloned_data_idx);
    }
    Pointer<HierarchyDataOpsReal<double> > hier_data_ops =
        HierarchyDataOpsManager::getManager()->getOperationsDouble(var,
                                                                         hierarchy,
                                                                         /* get_unique */ true);
    if (!hier_data_ops)
    {
        TBOX_ERROR(d_object_name << "::setDataOnPatchHierarchy():\n"
                                 << "  unsupported data centering.\n");
    }
    hier_data_ops->resetLevels(coarsest_ln, finest_ln);
    TBOX_ASSERT(!d_fcns.empty());
    d_fcns[0]->setDataOnPatchHierarchy(data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);
    for (unsigned int k = 1; k < d_fcns.size(); ++k)
    {
        d_fcns[k]->setDataOnPatchHierarchy(
            cloned_data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);
        hier_data_ops->add(data_idx, data_idx, cloned_data_idx);
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(cloned_data_idx);
    }
    var_db->removePatchDataIndex(cloned_data_idx);
    return;
} // setDataOnPatchHierarchy

void CartGridFunctionSet::setDataOnPatchLevel(const int data_idx,
                                              Pointer<Variable > var,
                                              Pointer<PatchLevel > level,
                                              const double data_time,
                                              const bool initial_time)
{
    TBOX_ASSERT(level);
    Pointer<CellVariable<double> > cc_var = var;
    Pointer<EdgeVariable<double> > ec_var = var;
    Pointer<FaceVariable<double> > fc_var = var;
    Pointer<NodeVariable<double> > nc_var = var;
    Pointer<SideVariable<double> > sc_var = var;
    TBOX_ASSERT(cc_var || ec_var || fc_var || nc_var || sc_var);
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    const int cloned_data_idx = var_db->registerClonedPatchDataIndex(var, data_idx);
    level->allocatePatchData(cloned_data_idx);
    TBOX_ASSERT(!d_fcns.empty());
    d_fcns[0]->setDataOnPatchLevel(data_idx, var, level, data_time, initial_time);
    for (unsigned int k = 1; k < d_fcns.size(); ++k)
    {
        d_fcns[k]->setDataOnPatchLevel(cloned_data_idx, var, level, data_time, initial_time);
        for (PatchLevel::Iterator p(level); p; p++)
        {
            Pointer<Patch > patch = p();
            if (cc_var)
            {
                Pointer<CellData<double> > data = patch->getPatchData(data_idx);
                Pointer<CellData<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchCellDataBasicOps<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (ec_var)
            {
                Pointer<EdgeData<double> > data = patch->getPatchData(data_idx);
                Pointer<EdgeData<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchEdgeDataBasicOps<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (fc_var)
            {
                Pointer<FaceData<double> > data = patch->getPatchData(data_idx);
                Pointer<FaceData<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchFaceDataBasicOps<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (nc_var)
            {
                Pointer<NodeData<double> > data = patch->getPatchData(data_idx);
                Pointer<NodeData<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchNodeDataBasicOps<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (sc_var)
            {
                Pointer<SideData<double> > data = patch->getPatchData(data_idx);
                Pointer<SideData<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchSideDataBasicOps<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else
            {
                TBOX_ERROR(d_object_name << "::setDataOnPatchLevel():\n"
                                         << "  unsupported data centering.\n");
            }
        }
    }
    level->deallocatePatchData(cloned_data_idx);
    var_db->removePatchDataIndex(cloned_data_idx);
    return;
} // setDataOnPatchLevel

void CartGridFunctionSet::setDataOnPatch(int data_idx,
                                         Pointer<Variable > var,
                                         Pointer<Patch > patch,
                                         double data_time,
                                         bool initial_time,
                                         Pointer<PatchLevel > patch_level)
{
    TBOX_ASSERT(patch);
    Pointer<CellVariable<double> > cc_var = var;
    Pointer<EdgeVariable<double> > ec_var = var;
    Pointer<FaceVariable<double> > fc_var = var;
    Pointer<NodeVariable<double> > nc_var = var;
    Pointer<SideVariable<double> > sc_var = var;
    TBOX_ASSERT(cc_var || ec_var || fc_var || nc_var || sc_var);
    Pointer<PatchData > data = patch->getPatchData(data_idx);
    Pointer<PatchData > cloned_data;
    if (cc_var)
    {
        Pointer<CellData<double> > p_data = data;
        cloned_data = new CellData<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (ec_var)
    {
        Pointer<EdgeData<double> > p_data = data;
        cloned_data = new EdgeData<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (fc_var)
    {
        Pointer<FaceData<double> > p_data = data;
        cloned_data = new FaceData<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (nc_var)
    {
        Pointer<NodeData<double> > p_data = data;
        cloned_data = new NodeData<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (sc_var)
    {
        Pointer<SideData<double> > p_data = data;
        cloned_data = new SideData<double>(
            p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth(), p_data->getDirectionVector());
    }
    else
    {
        TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                                 << "  unsupported data centering.\n");
    }
    cloned_data->setTime(data->getTime());
    TBOX_ASSERT(!d_fcns.empty());
    d_fcns[0]->setDataOnPatch(data_idx, var, patch, data_time, initial_time, patch_level);
    cloned_data->copy(*data);
    // NOTE: We operate on data_idx instead of cloned_data_idx here because it
    // is not straightforward to add a cloned data index to a single patch.
    for (unsigned int k = 1; k < d_fcns.size(); ++k)
    {
        d_fcns[k]->setDataOnPatch(data_idx, var, patch, data_time, initial_time, patch_level);
        if (cc_var)
        {
            Pointer<CellData<double> > p_data = data;
            Pointer<CellData<double> > p_cloned_data = cloned_data;
            PatchCellDataBasicOps<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (ec_var)
        {
            Pointer<EdgeData<double> > p_data = data;
            Pointer<EdgeData<double> > p_cloned_data = cloned_data;
            PatchEdgeDataBasicOps<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (fc_var)
        {
            Pointer<FaceData<double> > p_data = data;
            Pointer<FaceData<double> > p_cloned_data = cloned_data;
            PatchFaceDataBasicOps<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (nc_var)
        {
            Pointer<NodeData<double> > p_data = data;
            Pointer<NodeData<double> > p_cloned_data = cloned_data;
            PatchNodeDataBasicOps<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (sc_var)
        {
            Pointer<SideData<double> > p_data = data;
            Pointer<SideData<double> > p_cloned_data = cloned_data;
            PatchSideDataBasicOps<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else
        {
            TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                                     << "  unsupported data centering.\n");
        }
    }
    data->copy(*cloned_data);
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////