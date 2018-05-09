
#ifndef SUTTERLINBASEMENTMEMBRANEFORCE_HPP_
#define SUTTERLINBASEMENTMEMBRANEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

/**
 * A force proposed between each ellipsoidal cell and a basement membrane.
 * For further details, see:
 *
 * Sütterlin T, Tsingos E, Bensaci J, Stamatas GN, Grabe N. A 3D self-organizing
 * multicellular epidermis model of barrier formation and hydration with
 * realistic cell morphology based on EPISIM. Scientific Reports. 2017 Mar 6;
 * 7:43472. doi:10.1038/srep43472
 */
template<unsigned DIM>
class SutterlinBasementMembraneForce  : public AbstractForce<DIM>
{
private:

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mDeltaAdh;
        archive & mKAdh;
        archive & mKCBm;
    }

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

    /**
     * Constructor.
     */
    SutterlinBasementMembraneForce();

    /**
     * Destructor.
     */
    ~SutterlinBasementMembraneForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SutterlinBasementMembraneForce)

#endif /*SUTTERLINBASEMENTMEMBRANEFORCE_HPP_*/
