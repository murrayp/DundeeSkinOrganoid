
#include "EllipsoidNodeBasedCellPopulation.hpp"
#include "EllipsoidNodeAttributes.hpp"
#include "SkinOrganoidProperty.hpp"
#include "Debug.hpp"




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
CellPtr EllipsoidNodeBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, CellPtr pParentCell)
{

    auto pNewCellTemp=NodeBasedCellPopulation<DIM>::AddCell(pNewCell, pParentCell);

    // Get new node
    Node<DIM>* p_new_node = this->GetNodeCorrespondingToCell(pNewCellTemp);// new Node<DIM>(this->GetNumNodes(), daughter_position, false); // never on boundary

    p_new_node->AddNodeAttribute(0.0);

    double semi_major = (this->GetNodeCorrespondingToCell(pParentCell))->rGetNodeAttributes()[NA_SEMIMAJORAXIS];
    double semi_minor = (this->GetNodeCorrespondingToCell(pParentCell))->rGetNodeAttributes()[NA_SEMIMINORAXIS];


    p_new_node->rGetNodeAttributes()[NA_SEMIMAJORAXIS] = semi_major;
    p_new_node->rGetNodeAttributes()[NA_SEMIMINORAXIS] = semi_minor;


    // Get this cell's type six machine property data
    CellPropertyCollection collection = pParentCell->rGetCellPropertyCollection().template GetProperties<SkinOrganoidProperty>();
    if (collection.GetSize() != 1)
    {
        EXCEPTION("TypeSixMachineCellKiller cannot be used unless each cell has a SkinOrganoidProperty");
    }
    boost::shared_ptr<SkinOrganoidProperty> p_parent_property = boost::static_pointer_cast<SkinOrganoidProperty>(collection.GetProperty());
    unsigned parent_data = p_parent_property->GetCellDifferentiatedType();



    // remove copied daughter cell property and create a new cell property
    pNewCellTemp->template RemoveCellProperty<SkinOrganoidProperty>();


    MAKE_PTR(SkinOrganoidProperty, p_property);
    p_property->SetCellDifferentiatedType(parent_data);

    pNewCellTemp->AddCellProperty(p_property);




    return pNewCellTemp;

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
