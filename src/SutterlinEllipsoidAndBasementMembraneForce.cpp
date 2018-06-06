
#include "SutterlinEllipsoidAndBasementMembraneForce.hpp"
#include "EllipsoidNodeBasedCellPopulation.hpp"
#include "EllipsoidNodeAttributes.hpp"
#include "Debug.hpp"


template<unsigned DIM>
SutterlinEllipsoidAndBasementMembraneForce<DIM>::SutterlinEllipsoidAndBasementMembraneForce()
   : AbstractForce<DIM>(),
     mDeltaOl(1.0),
     mDeltaOlMax(0.2),
     mDOlMin(0.0), // micrometres
     mKPr(25.2e1), // N m^{-1}
     mDeltaAdh(0.2),
     mKAdh(25.0), // N m^{-1}
	 mKCBm(1.0)
{
	///\todo work out if any of these default parameter values need to be rescaled
}

template<unsigned DIM>
void SutterlinEllipsoidAndBasementMembraneForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // This force class is defined for EllipsoidNodeBasedCellPopulation only
    assert(dynamic_cast<EllipsoidNodeBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr);

    // Calculate forces between neighbouring cells
    AbstractCentreBasedCellPopulation<DIM,DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<DIM,DIM>*>(&rCellPopulation);

    std::vector< std::pair<Node<DIM>*, Node<DIM>* > >& r_node_pairs = p_static_cast_cell_population->rGetNodePairs();

    for (typename std::vector< std::pair<Node<DIM>*, Node<DIM>* > >::iterator iter = r_node_pairs.begin();
        iter != r_node_pairs.end();
        iter++)
    {
        std::pair<Node<DIM>*, Node<DIM>* > pair = *iter;

        unsigned node_a_index = pair.first->GetIndex();
        unsigned node_b_index = pair.second->GetIndex();

        // Calculate the force between nodes
        c_vector<double, DIM> force = CalculateForceBetweenNodes(node_a_index, node_b_index, rCellPopulation);
        for (unsigned j=0; j<DIM; j++)
        {
            assert(!std::isnan(force[j]));
        }

        // Add the force contribution to each node
        c_vector<double, DIM> negative_force = -1.0*force;
        pair.first->AddAppliedForceContribution(force);
        pair.second->AddAppliedForceContribution(negative_force);
    }

    // Calculate forces between cells and the basement membrane
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
				temp += v_c_bm(1)*v_c_bm(1)/(semi_major_axis*semi_major_axis);
			}
			if (DIM > 2)
			{
				temp += v_c_bm(2)*v_c_bm(2)/(semi_minor_axis*semi_minor_axis);
			}

			double d_opt = norm_2(v_c_bm/temp);
			double r_c = mDeltaAdh*d_opt;
			double d_adh = r_c; ///\todo check this...

			double A_con = M_PI*r_c*(r_c - d_c_bm);

//			double d_gap = d_c_bm - d_opt;
			double d_hat_gap = 1; ///\todo check how we should interpret equation (9) in Sutterlin et al in the basement membrane case

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
c_vector<double, DIM> SutterlinEllipsoidAndBasementMembraneForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                               unsigned nodeBGlobalIndex,
                                                                               AbstractCellPopulation<DIM>& rCellPopulation)
{
    // This force class is defined for EllipsoidNodeBasedCellPopulation only
    assert(dynamic_cast<EllipsoidNodeBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr);

    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    const double semi_major_axis_a = p_node_a->rGetNodeAttributes()[NA_SEMIMAJORAXIS];
    const double semi_minor_axis_a = p_node_a->rGetNodeAttributes()[NA_SEMIMINORAXIS];

    Node<DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);
    const double semi_major_axis_b = p_node_b->rGetNodeAttributes()[NA_SEMIMAJORAXIS];
    const double semi_minor_axis_b = p_node_b->rGetNodeAttributes()[NA_SEMIMINORAXIS];

    ///\todo Vary k_cn according to the types of these cells
    double k_cn = 1.0;

    // Get the node locations
    const c_vector<double, DIM>& r_node_a_location = p_node_a->rGetLocation();
    const c_vector<double, DIM>& r_node_b_location = p_node_b->rGetLocation();

    // Get the unit vector parallel to the line joining the two nodes (assuming no periodicities etc.)
    // Note: This is called v_cn_hat in Sutterlin et al.
    c_vector<double, DIM> unit_vector = r_node_b_location - r_node_a_location;

    // Calculate the distance between the two nodes
    // Note: This is called ||v_cn|| in Sutterlin et al.
    double distance_between_nodes = norm_2(unit_vector);

    unit_vector /= distance_between_nodes;





    double temp_a = unit_vector(0)*unit_vector(0)/(pow(semi_major_axis_a,2));
    double temp_b = unit_vector(0)*unit_vector(0)/(pow(semi_major_axis_b,2));
    if (DIM > 1)
    {
    	temp_a += unit_vector(1)*unit_vector(1)/(pow(semi_major_axis_a,2));
    	temp_b += unit_vector(1)*unit_vector(1)/(pow(semi_major_axis_b,2));
    }
    if (DIM > 2)
    {
    	temp_a += unit_vector(2)*unit_vector(2)/(pow(semi_minor_axis_a,2));
    	temp_b += unit_vector(2)*unit_vector(2)/(pow(semi_minor_axis_b,2));
    }


    double d_opt = norm_2(unit_vector/temp_a) + norm_2(-unit_vector/temp_b);
    double d_hat_opt = mDeltaOl*d_opt;

    double d_ol = d_hat_opt - distance_between_nodes;
    double d_ol_max = mDeltaOlMax*d_hat_opt;

    // Compute overlapping force magnitude
    double F_pr = 0.0;
    if ((d_ol >= mDOlMin) && (d_ol < d_ol_max))
    {
    	F_pr = mKPr*d_ol;
    }
    else if (d_ol >= d_ol_max)
    {
    	F_pr = mKPr*d_ol_max*exp(1.0*(d_ol/d_ol_max - 1));
    }
    //PRINT_VARIABLE(d_ol);
    //PRINT_VARIABLE(F_pr);

    // Compute adhesion force magnitude
    double d_gap = distance_between_nodes - d_opt;

    // Note: It is unclear from Sutterlin et al how d_seg_cn and d_seg_nc are actually defined,
    // but they look to be the semi-major axes of the neighouring cells
    double d_seg_cn = semi_major_axis_a;
    double d_seg_nc = semi_major_axis_b;

    double d_hat_gap = 1.0;
    if (((1 - mDeltaAdh)*d_seg_cn < d_gap) && (d_gap < (mDeltaAdh - 1)*d_seg_cn))
    {
    	d_hat_gap = fabs(sin(0.5*M_PI*d_gap/((mDeltaAdh - 1)*d_seg_cn)));
    }

    double d_adh = mDeltaAdh*d_opt;
    double r_c = mDeltaAdh*d_seg_cn;
    double r_n = mDeltaAdh*d_seg_nc;

    double v_nc_sq = pow(distance_between_nodes,2);
    double A_adh = (0.25*M_PI/v_nc_sq)*(2*v_nc_sq*(pow(r_c,2) + pow(r_n,2)) + 2*(pow(r_c,2))*(pow(r_n,2)) - pow(r_c,4) - pow(r_n,4) - pow(v_nc_sq,2));

    double F_adh = 0.0;
    if ((d_ol < 0) && (distance_between_nodes < d_adh))
    {
    	F_adh = mKAdh*k_cn*d_hat_gap*A_adh;
    }

    c_vector<double, DIM> F = (F_adh - F_pr)*unit_vector;

    return F;
}

template<unsigned DIM>
void SutterlinEllipsoidAndBasementMembraneForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DeltaOl>" << mDeltaOl << "</DeltaOl>\n";
    *rParamsFile << "\t\t\t<DeltaOlMax>" << mDeltaOlMax << "</DeltaOlMax>\n";
    *rParamsFile << "\t\t\t<DOlMin>" << mDOlMin << "</DOlMin>\n";
    *rParamsFile << "\t\t\t<KPr>" << mKPr << "</KPr>\n";
    *rParamsFile << "\t\t\t<DeltaAdh>" << mDeltaAdh << "</DeltaAdh>\n";
    *rParamsFile << "\t\t\t<KAdh>" << mKAdh << "</KAdh>\n";
    *rParamsFile << "\t\t\t<KCBm>" << mKCBm << "</KCBm>\n";

	// Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SutterlinEllipsoidAndBasementMembraneForce<1>;
template class SutterlinEllipsoidAndBasementMembraneForce<2>;
template class SutterlinEllipsoidAndBasementMembraneForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SutterlinEllipsoidAndBasementMembraneForce)
