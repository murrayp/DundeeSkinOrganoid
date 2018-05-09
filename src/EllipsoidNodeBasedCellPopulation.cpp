
#include "EllipsoidNodeBasedCellPopulation.hpp"

template<unsigned DIM>
EllipsoidNodeBasedCellPopulation<DIM>::EllipsoidNodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh,
                                                                        std::vector<CellPtr>& rCells,
                                                                        const std::vector<unsigned> locationIndices,
                                                                        bool deleteMesh)
    : NodeBasedCellPopulation<DIM>(rMesh,
                                   rCells,
								   locationIndices,
								   deleteMesh)
{
}

template<unsigned DIM>
EllipsoidNodeBasedCellPopulation<DIM>::EllipsoidNodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh)
    : NodeBasedCellPopulation<DIM>(rMesh)
{
    // No Validate() because the cells are not associated with the cell population yet in archiving
}

template<unsigned DIM>
void EllipsoidNodeBasedCellPopulation<DIM>::UpdateNodeLocations(double dt)
{
    NodeBasedCellPopulation<DIM>::UpdateNodeLocations(dt);
}

template<unsigned DIM>
CellPtr EllipsoidNodeBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell,
		                                               CellPtr pParentCell)
{
    auto p_new_cell_temp = NodeBasedCellPopulation<DIM>::AddCell(pNewCell, pParentCell);
    return p_new_cell_temp;
}

template<unsigned DIM>
void EllipsoidNodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

// Explicit instantiation
template class EllipsoidNodeBasedCellPopulation<1>;
template class EllipsoidNodeBasedCellPopulation<2>;
template class EllipsoidNodeBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EllipsoidNodeBasedCellPopulation)
