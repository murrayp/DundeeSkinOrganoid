
#ifndef ELLIPSOIDMODIFIER_HPP_
#define ELLIPSOIDMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

///\todo Document class
template<unsigned DIM>
class EllipsoidModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:

    ///\todo Implement archiving

    /** \todo Document member */
    std::string mOutputDirectory;

    /** Meta results file for VTK. */
    out_stream mpVtkMetaFile;

public:

    /**
     * Constructor.
     */
    EllipsoidModifier();

    /**
     * Destructor.
     */
    virtual ~EllipsoidModifier();

    void SetOutputDirectory(std::string directory);

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
    {
    }

    /**
     * Overridden UpdateAtEndOfOutputTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each output timestep,
     * after UpdateAtEndOfTimeStep() has been called.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /** \todo Document method */
    void WriteVtk(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Specify what to do in the simulation at the end of each time loop.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EllipsoidModifier)

#endif /*ELLIPSOIDMODIFIER_HPP_*/