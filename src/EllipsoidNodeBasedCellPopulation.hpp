
#ifndef ELLIPSOIDNODEBASEDCELLPOPULAITON_HPP_
#define ELLIPSOIDNODEBASEDCELLPOPULAITON_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "NodeBasedCellPopulation.hpp"

/**
 * The class EllipsoidNodeBasedCellPopulation extends NodeBasedCellPopulation to
 * account for ellipsoidal cells. For further details, see:
 *
 * Sütterlin T, Tsingos E, Bensaci J, Stamatas GN, Grabe N. A 3D self-organizing
 * multicellular epidermis model of barrier formation and hydration with
 * realistic cell morphology based on EPISIM. Scientific Reports. 2017 Mar 6;
 * 7:43472. doi:10.1038/srep43472
 */
template<unsigned DIM>
class EllipsoidNodeBasedCellPopulation : public NodeBasedCellPopulation<DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the nodes is handled by load/save_construct_data,
     * so we don't actually have to do anything here except delegate to the base class.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<NodeBasedCellPopulation<DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     *
     * Note that the cell population will take responsibility for freeing the memory used by the nodes.
     *
     * @param rMesh a mutable nodes-only mesh
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh whether to delete nodes-only mesh in destructor
     */
    EllipsoidNodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh,
                                     std::vector<CellPtr>& rCells,
                                     const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                                     bool deleteMesh=false);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable nodes-only mesh
     */
    EllipsoidNodeBasedCellPopulation(NodesOnlyMesh<DIM>& rMesh);

    /**
     * Overridden UpdateNodeLocations() method.
     *
     * @param dt time step
     */
    void UpdateNodeLocations(double dt);

    /**
     * Overridden AddCell() method.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pNewCell the cell to add
     * @param pParentCell pointer to a parent cell
     *
     * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
     */
    CellPtr AddCell(CellPtr pNewCell, CellPtr pParentCell=CellPtr());

    /**
     * Overridden OutputCellPopulationParameters() method.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EllipsoidNodeBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an EllipsoidNodeBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
	Archive & ar, const EllipsoidNodeBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
	// Save data required to construct instance
	const NodesOnlyMesh<DIM>* p_mesh = &(t->rGetMesh());
	ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise an EllipsoidNodeBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
	Archive & ar, EllipsoidNodeBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
	// Retrieve data from archive required to construct new instance
	NodesOnlyMesh<DIM>* p_mesh;
	ar >> p_mesh;

	// Invoke inplace constructor to initialise instance
	::new(t)EllipsoidNodeBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*ELLIPSOIDNODEBASEDCELLPOPULAITON_HPP_*/
