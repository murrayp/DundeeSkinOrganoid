#include "SkinOrganoidProperty.hpp"
#include "Debug.hpp"


SkinOrganoidProperty::SkinOrganoidProperty()
    : AbstractCellProperty(),mCellDifferentiatedType(0u),mIntraCellularCalcium(0.0),md_IntraCellularCalciumd_t(0.0)
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

double SkinOrganoidProperty::GetD_IntraCellularCalciumD_t()
{
    return md_IntraCellularCalciumd_t;
}

void SkinOrganoidProperty::SetD_IntraCellularCalciumD_t(double d_intraCellularCalciumd_t)
{
    md_IntraCellularCalciumd_t=d_intraCellularCalciumd_t;
}



#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SkinOrganoidProperty)
