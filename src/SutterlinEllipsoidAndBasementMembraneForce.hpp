
#ifndef SUTTERLINELLIPSOIDANDBASEMENTMEMBRANEFORCE
#define SUTTERLINELLIPSOIDANDBASEMENTMEMBRANEFORCE

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractTwoBodyInteractionForce.hpp"

/**
 * A two-body interaction force proposed for neighbouring ellipsoidal cells.
 * For further details, see:
 *
 * Sütterlin T, Tsingos E, Bensaci J, Stamatas GN, Grabe N. A 3D self-organizing
 * multicellular epidermis model of barrier formation and hydration with
 * realistic cell morphology based on EPISIM. Scientific Reports. 2017 Mar 6;
 * 7:43472. doi:10.1038/srep43472
 */
template<unsigned DIM>
class SutterlinEllipsoidAndBasementMembraneForce : public AbstractForce<DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mDeltaOl;
        archive & mDeltaOlMax;
        archive & mDOlMin;
        archive & mKPr;
        archive & mDeltaAdh;
        archive & mKAdh;
        archive & mKCBm;
    }

    /**
     * Parameter delta_ol in Sütterlin et al.
     * Defaults to the value 0.15.
     */
    double mDeltaOl;

    /**
     * Parameter delta_ol_max in Sütterlin et al.
     * Defaults to the value 0.5.
     */
    double mDeltaOlMax;

    /**
     * Parameter d_ol_min in Sütterlin et al.
     * Defaults to the value 0.1 micrometres.
     */
    double mDOlMin;

    /**
     * Parameter k_pr in Sütterlin et al.
     * Defaults to the value 2.2e-3 N m^{-1}.
     */
    double mKPr;

    /**
     * Parameter delta_adh in Sütterlin et al.
     * Defaults to the value 1.3.
     */
    double mDeltaAdh;

    /**
     * Parameter k_adh in Sütterlin et al.
     * Defaults to the value 2.2e-3 N m^{-1}.
     */
    double mKAdh;

    /**
     * Parameter k_c_bm in Sütterlin et al.
     * Defaults to the value 0.01.
     */
    double mKCBm;

public:

    ///\todo add get and set methods for these parameters

    /**
     * Constructor.
     */
    SutterlinEllipsoidAndBasementMembraneForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * @return the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode().
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     */
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SutterlinEllipsoidAndBasementMembraneForce)

#endif /*SUTTERLINELLIPSOIDANDBASEMENTMEMBRANEFORCE*/
