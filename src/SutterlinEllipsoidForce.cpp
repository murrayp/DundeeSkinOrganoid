
#include "SutterlinEllipsoidForce.hpp"
#include "EllipsoidNodeBasedCellPopulation.hpp"
#include "EllipsoidNodeAttributes.hpp"

template<unsigned DIM>
SutterlinEllipsoidForce<DIM>::SutterlinEllipsoidForce()
   : AbstractTwoBodyInteractionForce<DIM>(),
     mDeltaOl(0.15),
     mDeltaOlMax(0.5),
     mDOlMin(0.1), // micrometres
     mKPr(2.2e-3), // N m^{-1}
     mDeltaAdh(1.3),
     mKAdh(2.2e-3) // N m^{-1}
{
	///\todo work out if any of these default parameter values need to be rescaled
}

template<unsigned DIM>
c_vector<double, DIM> SutterlinEllipsoidForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
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

    double temp_a = unit_vector(0)*unit_vector(0)/(semi_major_axis_a*semi_major_axis_a);
    double temp_b = unit_vector(0)*unit_vector(0)/(semi_major_axis_b*semi_major_axis_b);
    if (DIM > 1)
    {
    	temp_a += unit_vector(1)*unit_vector(1)/(semi_minor_axis_a*semi_minor_axis_a);
    	temp_b += unit_vector(1)*unit_vector(1)/(semi_minor_axis_b*semi_minor_axis_b);
    }
    if (DIM > 2)
    {
    	temp_a += unit_vector(2)*unit_vector(2)/(semi_major_axis_a*semi_major_axis_a);
    	temp_b += unit_vector(2)*unit_vector(2)/(semi_major_axis_b*semi_major_axis_b);
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
    	F_pr = mKPr*d_ol_max*exp(d_ol/d_ol_max - 1);
    }

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
void SutterlinEllipsoidForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DeltaOl>" << mDeltaOl << "</DeltaOl>\n";
    *rParamsFile << "\t\t\t<DeltaOlMax>" << mDeltaOlMax << "</DeltaOlMax>\n";
    *rParamsFile << "\t\t\t<DOlMin>" << mDOlMin << "</DOlMin>\n";
    *rParamsFile << "\t\t\t<KPr>" << mKPr << "</KPr>\n";
    *rParamsFile << "\t\t\t<DeltaAdh>" << mDeltaAdh << "</DeltaAdh>\n";
    *rParamsFile << "\t\t\t<KAdh>" << mKAdh << "</KAdh>\n";

	// Call method on direct parent class
    AbstractTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SutterlinEllipsoidForce<1>;
template class SutterlinEllipsoidForce<2>;
template class SutterlinEllipsoidForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SutterlinEllipsoidForce)
