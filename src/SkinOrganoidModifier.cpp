
#include "SkinOrganoidModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "SkinOrganoidProperty.hpp"
#include "Exception.hpp"
#include "VtkMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "NodesOnlyMesh.hpp"
#include "UniformCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "EllipsoidNodeAttributes.hpp"

#include "EllipsoidNodeBasedCellPopulation.hpp"




template<unsigned DIM>
SkinOrganoidModifier<DIM>::SkinOrganoidModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}


template<unsigned DIM>
SkinOrganoidModifier<DIM>::~SkinOrganoidModifier()
{
}

template<unsigned DIM>
void SkinOrganoidModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void SkinOrganoidModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}



template<unsigned DIM>
void SkinOrganoidModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{


    ///\todo Make sure the cell population is updated?
    rCellPopulation.Update();
    //double dt = SimulationTime::Instance()->GetTimeStep();


    // Iterate over cell population and update machines
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    { 
        // Get this cell's type six machine property data
        CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<SkinOrganoidProperty>();
        if (collection.GetSize() != 1)
        {
            EXCEPTION("SkinOrganoidModifier cannot be used unless each cell has a SkinOrganoidProperty");
        }
        boost::shared_ptr<SkinOrganoidProperty> p_property = boost::static_pointer_cast<SkinOrganoidProperty>(collection.GetProperty());
        //unsigned cell_differentiated_type = p_property->GetCellDifferentiatedType();
        //unsigned& r_NumMachineFiresInThisTimeStep = p_property->rGetNumMachineFiresInThisTimeStep();

        // first  model for assigning cell types based on height from basal layer
        c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        double cell_height= cell_location(2);


        if (cell_height <1.0) //basal
        {
            p_property->SetCellDifferentiatedType(0u);

            boost::shared_ptr<AbstractCellProperty> p_stem_type =
                    cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<StemCellProliferativeType>();
            cell_iter->SetCellProliferativeType(p_stem_type);

        }
        else if (cell_height < 10.0) // spinosal
        {
            p_property->SetCellDifferentiatedType(1u);
            //if (cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
            {
                boost::shared_ptr<AbstractCellProperty> p_diff_type =
                                    cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<DifferentiatedCellProliferativeType>();
                            cell_iter->SetCellProliferativeType(p_diff_type);
            }


            // mess about with semi major axis

            //EllipsoidNodeBasedCellPopulation<DIM>& rellipsoid_pop=(dynamic_cast<EllipsoidNodeBasedCellPopulation<DIM>&>(rCellPopulation));

                EllipsoidNodeBasedCellPopulation<DIM>* p_ellipsoid_pop=(static_cast<EllipsoidNodeBasedCellPopulation<DIM>*>(&rCellPopulation));


                   Node<DIM>* p_node = p_ellipsoid_pop->GetNodeCorrespondingToCell(*cell_iter);


                   //double semi_major = p_node->rGetNodeAttributes()[NA_SEMIMAJORAXIS];
                   //double semi_minor = p_node->rGetNodeAttributes()[NA_SEMIMINORAXIS];

                   double semi_minor=4.0;
                   double semi_major=20.0;
                   p_node->AddNodeAttribute(0.0);
                   p_node->rGetNodeAttributes()[NA_SEMIMAJORAXIS] = semi_major;
                   p_node->rGetNodeAttributes()[NA_SEMIMINORAXIS] = semi_minor;

        }
        else if (cell_height < 14.0) // granular
        {
             p_property->SetCellDifferentiatedType(2u);
        }
        else  // Cornified
        {
            p_property->SetCellDifferentiatedType(3u);
        }
    }
}

    
template<unsigned DIM>
void SkinOrganoidModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class SkinOrganoidModifier<1>;
template class SkinOrganoidModifier<2>;
template class SkinOrganoidModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SkinOrganoidModifier)

