#include "SkinOrganoidProperty.hpp"
#include "Debug.hpp"


SkinOrganoidProperty::SkinOrganoidProperty()
    : AbstractCellProperty(),mCellDifferentiatedType(0u)
{
}

SkinOrganoidProperty::~SkinOrganoidProperty()
{
}

unsigned SkinOrganoidProperty::GetCellDifferentiatedType()
{
    return mCellDifferentiatedType;
}

void SkinOrganoidProperty::SetCellDifferentiatedType(unsigned cell_differentiated_type)
{
    mCellDifferentiatedType=cell_differentiated_type;
}

double SkinOrganoidProperty::GetIntraCellularCalcium()
{
    return mIntraCellularCalcium;
}

void SkinOrganoidProperty::SetIntraCellularCalcium(double intraCellularCalcium)
{
    mIntraCellularCalcium=intraCellularCalcium;
}



#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SkinOrganoidProperty)
