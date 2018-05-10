
#include "SutterlinBasementMembraneForce.hpp"
#include "EllipsoidNodeBasedCellPopulation.hpp"
#include "EllipsoidNodeAttributes.hpp"

template<unsigned DIM>
SutterlinBasementMembraneForce<DIM>::SutterlinBasementMembraneForce()
    : AbstractForce<DIM>(),
	  mDeltaAdh(1.3),
	  mKAdh(2.2e-3), // N m^{-1}
	  mKCBm(0.01)
{
}

template<unsigned DIM>
SutterlinBasementMembraneForce<DIM>::~SutterlinBasementMembraneForce()
{
}

template<unsigned DIM>
void SutterlinBasementMembraneForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // This force class is defined for EllipsoidNodeBasedCellPopulation only
    assert(dynamic_cast<EllipsoidNodeBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

		unsigned node_global_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

	    Node<DIM>* p_node = rCellPopulation.GetNode(node_global_index);
	    const double semi_major_axis = p_node->rGetNodeAttributes()[NA_SEMIMAJORAXIS];
	    const double semi_minor_axis = p_node->rGetNodeAttributes()[NA_SEMIMINORAXIS];
	    const c_vector<double, DIM>& r_c = p_node->rGetLocation();

    	///\todo only add this force to basal cells!
    	if (r_c(DIM-1) < 0.3)
    	{
			c_vector<double, DIM> r_bm = r_c;
			r_bm(DIM-1) = 0.0;

			c_vector<double, DIM> v_c_bm = r_bm - r_c;
			double d_c_bm = norm_2(v_c_bm);
			v_c_bm /= d_c_bm;

			double temp = v_c_bm(0)*v_c_bm(0)/(semi_major_axis*semi_major_axis);
			if (DIM > 1)
			{
				temp += v_c_bm(1)*v_c_bm(1)/(semi_minor_axis*semi_minor_axis);
			}
			if (DIM > 2)
			{
				temp += v_c_bm(2)*v_c_bm(2)/(semi_major_axis*semi_major_axis);
			}

			double d_opt = norm_2(v_c_bm/temp);
			double r_c = mDeltaAdh*d_opt;
			double d_adh = r_c; ///\todo check this...

			double A_con = M_PI*r_c*(r_c - d_c_bm);

//			double d_gap = d_c_bm - d_opt;
			double d_hat_gap = 1; ///\todo check how we should interpret equatoin (9) in Sutterlin et al in the basement membrane case

			double F_adh = 0.0;
			if (d_c_bm < d_adh)
			{
				F_adh = mKAdh*mKCBm*d_hat_gap*A_con; ///\todo Note: this is mistyped as A_adh by Sutterlin et al, we think
			}

			c_vector<double,DIM> force = F_adh*v_c_bm;

			rCellPopulation.GetNode(node_global_index)->AddAppliedForceContribution(force);
    	}
    }
}

template<unsigned DIM>
void SutterlinBasementMembraneForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DeltaAdh>" << mDeltaAdh << "</DeltaAdh>\n";
    *rParamsFile << "\t\t\t<KAdh>" << mKAdh << "</KAdh>\n";
    *rParamsFile << "\t\t\t<KCBm>" << mKCBm << "</KCBm>\n";

	// Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SutterlinBasementMembraneForce<1>;
template class SutterlinBasementMembraneForce<2>;
template class SutterlinBasementMembraneForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SutterlinBasementMembraneForce)
