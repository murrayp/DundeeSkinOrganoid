
#include "SkinOrganoidModifier.hpp"
#include "RandomNumberGenerator.hpp"
#include "SkinOrganoidProperty.hpp"
#include "Exception.hpp"
#include "VtkMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "NodesOnlyMesh.hpp"
#include "UniformCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "EllipsoidNodeAttributes.hpp"
#include "EllipsoidNodeBasedCellPopulation.hpp"

template<unsigned DIM>
SkinOrganoidModifier<DIM>::SkinOrganoidModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
	mBasalCellSemiMajorAndMinorAxis(0)=0.5;
	mBasalCellSemiMajorAndMinorAxis(1)=0.5;

	mSpinosalCellSemiMajorAndMinorAxis(0)=0.6;
	mSpinosalCellSemiMajorAndMinorAxis(1)=0.4;

	mGranularCellSemiMajorAndMinorAxis(0)=0.7;
	mGranularCellSemiMajorAndMinorAxis(1)=0.2;



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
    double dt = SimulationTime::Instance()->GetTimeStep();

    // Iterate over cell population and compute time derivatives
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

        // Update intracellular calcium

        double intracellular_calcium = p_property->GetIntraCellularCalcium();


        double d_intraCellularCalcium_dt=0.0;


        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);

       // Compute this cell's average neighbouring Delta concentration and store in CellData
       if (!neighbour_indices.empty())
       {
           for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                iter != neighbour_indices.end();
                ++iter)
           {

               CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);

               CellPropertyCollection collection_neighbour = p_cell->rGetCellPropertyCollection().template GetProperties<SkinOrganoidProperty>();
               boost::shared_ptr<SkinOrganoidProperty> p_property_neighbour = boost::static_pointer_cast<SkinOrganoidProperty>(collection_neighbour.GetProperty());
               double intracellular_calcium_neighbour = p_property->GetIntraCellularCalcium();


               double cell_wise_contribution;
               double transport_constant=10.1;


               cell_wise_contribution = transport_constant*(intracellular_calcium_neighbour - intracellular_calcium);



               d_intraCellularCalcium_dt+=cell_wise_contribution ;

           }


       }

       double intracellular_calcium_deg_rate=0.05;
       double reaction_term = -intracellular_calcium*intracellular_calcium_deg_rate;
       d_intraCellularCalcium_dt+=reaction_term;
       p_property->SetD_IntraCellularCalciumD_t(d_intraCellularCalcium_dt);


    }



	// Update Cell properties
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


	        double intracellular_calcium = p_property->GetIntraCellularCalcium();
            double d_intraCellularCalcium_dt = p_property->GetD_IntraCellularCalciumD_t();
            double new_intracellular_calcium = intracellular_calcium + d_intraCellularCalcium_dt*dt;
            p_property->SetIntraCellularCalcium(new_intracellular_calcium);

	        // First  model for assigning cell types based on height from basal layer
	        c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

	        double cell_height = cell_location(2);

	        EllipsoidNodeBasedCellPopulation<DIM>* p_ellipsoid_pop=(static_cast<EllipsoidNodeBasedCellPopulation<DIM>*>(&rCellPopulation));


	        double basal_height_threshold= 1.0*mBasalCellSemiMajorAndMinorAxis(1);
	        double spinosal_height_threshold= basal_height_threshold+2.0*3.0*mSpinosalCellSemiMajorAndMinorAxis(1);
	        double granular_height_threshold= spinosal_height_threshold+ 2.0*mGranularCellSemiMajorAndMinorAxis(1) ;


	        if (cell_height <basal_height_threshold) //basal
	        {
	        	p_property->SetCellDifferentiatedType(0u);

	        	//boost::shared_ptr<AbstractCellProperty> p_stem_type =
	        	//cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<StemCellProliferativeType>();
	        	//cell_iter->SetCellProliferativeType(p_stem_type);
	        	Node<DIM>* p_node = p_ellipsoid_pop->GetNodeCorrespondingToCell(*cell_iter);
	        	p_node->AddNodeAttribute(0.0);
	        	p_node->rGetNodeAttributes()[NA_SEMIMAJORAXIS] = mBasalCellSemiMajorAndMinorAxis(0);
	        	p_node->rGetNodeAttributes()[NA_SEMIMINORAXIS] = mBasalCellSemiMajorAndMinorAxis(1);

	        	p_property->SetIntraCellularCalcium(2.0);

	        }
	        else if (cell_height >basal_height_threshold && cell_height < spinosal_height_threshold && p_property->GetCellDifferentiatedType()==0u ) // spinosal
	        {
	        	p_property->SetCellDifferentiatedType(1u);

                boost::shared_ptr<AbstractCellProperty> p_diff_type =
                cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<DifferentiatedCellProliferativeType>();
                 cell_iter->SetCellProliferativeType(p_diff_type);
                 cell_iter->SetBirthTime(1e4);


                 Node<DIM>* p_node = p_ellipsoid_pop->GetNodeCorrespondingToCell(*cell_iter);
                 p_node->AddNodeAttribute(0.0);

                 p_node->rGetNodeAttributes()[NA_SEMIMAJORAXIS] = mSpinosalCellSemiMajorAndMinorAxis(0);
                 p_node->rGetNodeAttributes()[NA_SEMIMINORAXIS] = mSpinosalCellSemiMajorAndMinorAxis(1);

        }
        else if (cell_height> spinosal_height_threshold && cell_height < granular_height_threshold && p_property->GetCellDifferentiatedType()==1u) // granular
        {
             	 p_property->SetCellDifferentiatedType(2u);

             boost::shared_ptr<AbstractCellProperty> p_diff_type =
                            cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->template Get<DifferentiatedCellProliferativeType>();
                             cell_iter->SetCellProliferativeType(p_diff_type);


             Node<DIM>* p_node = p_ellipsoid_pop->GetNodeCorrespondingToCell(*cell_iter);

             p_node->AddNodeAttribute(0.0);
             p_node->rGetNodeAttributes()[NA_SEMIMAJORAXIS] = mGranularCellSemiMajorAndMinorAxis(0);
             p_node->rGetNodeAttributes()[NA_SEMIMINORAXIS] = mGranularCellSemiMajorAndMinorAxis(1);
        }
        else if (cell_height > granular_height_threshold&&p_property->GetCellDifferentiatedType()==2u)  // Cornified
        {
            p_property->SetCellDifferentiatedType(3u);
        	p_property->SetIntraCellularCalcium(0.0);

        }

    }
}


template<unsigned DIM>
void SkinOrganoidModifier<DIM>::SetBasalCellSemiMajorAndMinorAxis(c_vector<double, 2> basalCellSemiMajorAndMinorAxis)
{
	mBasalCellSemiMajorAndMinorAxis=basalCellSemiMajorAndMinorAxis;
}

template<unsigned DIM>
void SkinOrganoidModifier<DIM>::SetSpinosalCellSemiMajorAndMinorAxis(c_vector<double, 2> spinosalCellSemiMajorAndMinorAxis)
{
	mSpinosalCellSemiMajorAndMinorAxis=spinosalCellSemiMajorAndMinorAxis;
}

template<unsigned DIM>
void SkinOrganoidModifier<DIM>::SetGranularCellSemiMajorAndMinorAxis(c_vector<double, 2> granularCellSemiMajorAndMinorAxis)
{
	mGranularCellSemiMajorAndMinorAxis=granularCellSemiMajorAndMinorAxis;
}

template<unsigned DIM>
c_vector<double, 2> SkinOrganoidModifier<DIM>::GetBasalCellSemiMajorAndMinorAxis()
{
	return mBasalCellSemiMajorAndMinorAxis;
}

template<unsigned DIM>
c_vector<double, 2> SkinOrganoidModifier<DIM>::GetSpinosalCellSemiMajorAndMinorAxis()
{
	return mSpinosalCellSemiMajorAndMinorAxis;
}

template<unsigned DIM>
c_vector<double, 2> SkinOrganoidModifier<DIM>::GetGranularCellSemiMajorAndMinorAxis()
{
	return mGranularCellSemiMajorAndMinorAxis;
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
