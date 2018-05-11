
#include "CellEllipsoidWriter.hpp"
#include "EllipsoidNodeBasedCellPopulation.hpp"
#include "UblasVectorInclude.hpp"
#include "EllipsoidNodeAttributes.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellEllipsoidWriter<ELEMENT_DIM, SPACE_DIM>::CellEllipsoidWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("ellipsoid.dat")
{
    this->mVtkVectorCellDataName = "Ellipsoid";
    this->mOutputScalarData = false;
    this->mOutputVectorData = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> CellEllipsoidWriter<ELEMENT_DIM, SPACE_DIM>::GetVectorCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    c_vector<double, SPACE_DIM> ellipsoid = scalar_vector<double>(SPACE_DIM, DOUBLE_UNSET);

    if (dynamic_cast<EllipsoidNodeBasedCellPopulation<SPACE_DIM>*>(pCellPopulation))
    {
        unsigned node_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
        Node<SPACE_DIM>* p_node = pCellPopulation->GetNode(node_index);
        p_node->AddNodeAttribute(0.0);
        double a = p_node->rGetNodeAttributes()[NA_SEMIMAJORAXIS];
        double b = p_node->rGetNodeAttributes()[NA_SEMIMINORAXIS];

        ellipsoid(0) = a;
        if (SPACE_DIM > 1)
        {
        	ellipsoid(1) = a;
        }
        if (SPACE_DIM > 2)
        {
        	ellipsoid(2) = b; // c = a is assumed by Sutterlin et al
        }
    }
    return ellipsoid;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellEllipsoidWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    c_vector<double, SPACE_DIM> ellipsoid = GetVectorCellDataForVtkOutput(pCell, pCellPopulation);

    *this->mpOutStream << location_index << " " << cell_id << " ";
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_location[i] << " ";
    }
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << ellipsoid[i] << " ";
    }
}

// Explicit instantiation
template class CellEllipsoidWriter<1,1>;
template class CellEllipsoidWriter<1,2>;
template class CellEllipsoidWriter<2,2>;
template class CellEllipsoidWriter<1,3>;
template class CellEllipsoidWriter<2,3>;
template class CellEllipsoidWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellEllipsoidWriter)
