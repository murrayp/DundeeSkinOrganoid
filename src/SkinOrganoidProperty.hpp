
#ifndef SKINORGANOIDPROPERTY_HPP_
#define SKINORGANOIDPROPERTY_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"
#include "PetscTools.hpp"
#include <set>

/**
 * \todo Document class
 */
class SkinOrganoidProperty : public AbstractCellProperty
{
private:

    /** \todo Document member */
    unsigned mCellDifferentiatedType;

    double mIntraCellularCalcium;

    double md_IntraCellularCalciumd_t;


public:

    /**
     * Constructor.
     */
    SkinOrganoidProperty();

    /**
     * Destructor.
     */
    virtual ~SkinOrganoidProperty();

    /**
     * @return #mMachineData
     */
    unsigned GetCellDifferentiatedType();

    void SetCellDifferentiatedType(unsigned cell_differentiated_type);

    double GetIntraCellularCalcium();

    void SetIntraCellularCalcium(double intraCellularCalcium);

    double GetD_IntraCellularCalciumD_t();

    void SetD_IntraCellularCalciumD_t(double d_intraCellularCalcium_dt);







};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SkinOrganoidProperty)

#endif /* SKINORGANOIDPROPERTY_HPP_ */
