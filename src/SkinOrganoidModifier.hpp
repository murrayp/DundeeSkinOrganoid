
#ifndef SkinOrganoidModifier_HPP_
#define SkinOrganoidModifier_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * \todo Document class
 */
template<unsigned DIM>
class SkinOrganoidModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:

    ///\todo Implement archiving
	c_vector<double, 2> mBasalCellSemiMajorAndMinorAxis;
	c_vector<double, 2> mSpinosalCellSemiMajorAndMinorAxis;
	c_vector<double, 2> mGranularCellSemiMajorAndMinorAxis;








public:

    /**
     * Constructor.
     */
    SkinOrganoidModifier();

    /**
     * Destructor.
     */
    virtual ~SkinOrganoidModifier();

    /**
        * Overridden UpdateAtEndOfTimeStep() method.
        *
        * Specify what to do in the simulation at the end of each time step.
        *
        * @param rCellPopulation reference to the cell population
        */
       virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
        * Overridden SetupSolve() method.
        *
        * Specify what to do in the simulation before the start of the time loop.
        *
        * @param rCellPopulation reference to the cell population
        * @param outputDirectory the output directory, relative to where Chaste output is stored
        */
       virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);



    /**
     * Helper method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);


    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    void SetBasalCellSemiMajorAndMinorAxis(c_vector<double, 2>);
    void SetSpinosalCellSemiMajorAndMinorAxis(c_vector<double, 2>);
    void SetGranularCellSemiMajorAndMinorAxis(c_vector<double, 2>);

    c_vector<double, 2> GetBasalCellSemiMajorAndMinorAxis();
    c_vector<double, 2> GetSpinosalCellSemiMajorAndMinorAxis();
    c_vector<double, 2> GetGranularCellSemiMajorAndMinorAxis();


};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SkinOrganoidModifier)

#endif /*SkinOrganoidModifier_HPP_*/
