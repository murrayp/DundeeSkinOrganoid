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



#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SkinOrganoidProperty)
