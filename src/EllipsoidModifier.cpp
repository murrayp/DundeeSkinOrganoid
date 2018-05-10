
#include "EllipsoidModifier.hpp"
#include "Exception.hpp"
#include "VtkMeshWriter.hpp"
#include "EllipsoidNodeAttributes.hpp"
#include "OutputFileHandler.hpp"
#include "EllipsoidNodeBasedCellPopulation.hpp"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"

template<unsigned DIM>
EllipsoidModifier<DIM>::EllipsoidModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mOutputDirectory("")
{
}

template<unsigned DIM>
EllipsoidModifier<DIM>::~EllipsoidModifier()
{
}

template<unsigned DIM>
void EllipsoidModifier<DIM>::SetOutputDirectory(std::string directory)
{
    mOutputDirectory = directory;
}

template<unsigned DIM>
void EllipsoidModifier<DIM>::WriteVtk(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // This force class is defined for EllipsoidNodeBasedCellPopulation only
    assert(dynamic_cast<EllipsoidNodeBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr);

#ifdef CHASTE_VTK
    // Store the present time as a string
    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    // Create mesh writer for VTK output
    VtkMeshWriter<DIM, DIM> mesh_writer(mOutputDirectory,
                                        "ellipsoid_results_"+time.str(),
                                        false);

    // Create vector to store VTK cell data
    std::vector<c_vector<double, DIM> > ellipsoid_data;
    ellipsoid_data.resize(rCellPopulation.GetNumNodes());

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
        double a = p_node->rGetNodeAttributes()[NA_SEMIMAJORAXIS];
        double b = p_node->rGetNodeAttributes()[NA_SEMIMINORAXIS];

        c_vector<double, DIM> ellipsoid;
        ellipsoid(0) = a;
        if (DIM > 1)
        {
        	ellipsoid(1) = b;
        }
        if (DIM > 2)
        {
        	ellipsoid(2) = a; // c = a is assumed by Sutterlin et al
        }

        // Populate data for VTK
        ellipsoid_data[node_index] = ellipsoid;
    }

    mesh_writer.AddPointData("ellipsoids", ellipsoid_data);

    /*
     * At present, the VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK, then visualized
     * as glyphs in Paraview.
     */
    NodesOnlyMesh<DIM>* p_temp_mesh = static_cast<NodesOnlyMesh<DIM>*>(&(rCellPopulation.rGetMesh()));
//    ///\todo Consider removing hardcoding of "1.0" below
//    temp_mesh.ConstructNodesWithoutMesh(machine_nodes, 1.0);
    mesh_writer.WriteFilesUsingMesh(*p_temp_mesh);

    *mpVtkMetaFile << "        <DataSet timestep=\"";
    *mpVtkMetaFile << num_timesteps;
    *mpVtkMetaFile << "\" group=\"\" part=\"0\" file=\"ellipsoid_results_";
    *mpVtkMetaFile << num_timesteps;
    *mpVtkMetaFile << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

template<unsigned DIM>
void EllipsoidModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (mOutputDirectory == "")
    {
       EXCEPTION("SetOutputDirectory() must be called on a EllipsoidModifier before it is passed to a simulation");
    }
    
    WriteVtk(rCellPopulation);
}

template<unsigned DIM>
void EllipsoidModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
#ifdef CHASTE_VTK
    // Create output files for the visualizer
    double time_now = SimulationTime::Instance()->GetTime();
    std::ostringstream time_string;
    time_string << time_now;

    if (mOutputDirectory == "")
    {
        EXCEPTION("SetOutputDirectory() must be called on EllipsoidModifier");
    }
    mOutputDirectory += "/ellipsoid_results_from_time_" + time_string.str();

    OutputFileHandler output_file_handler(mOutputDirectory, false);
    mpVtkMetaFile = output_file_handler.OpenOutputFile("ellipsoid_results.pvd");
    *mpVtkMetaFile << "<?xml version=\"1.0\"?>\n";
    *mpVtkMetaFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *mpVtkMetaFile << "    <Collection>\n";
#endif //CHASTE_VTK

    WriteVtk(rCellPopulation);
}

template<unsigned DIM>
void EllipsoidModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
#ifdef CHASTE_VTK
    *mpVtkMetaFile << "    </Collection>\n";
    *mpVtkMetaFile << "</VTKFile>\n";
    mpVtkMetaFile->close();
#endif //CHASTE_VTK
}
    
template<unsigned DIM>
void EllipsoidModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class EllipsoidModifier<1>;
template class EllipsoidModifier<2>;
template class EllipsoidModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EllipsoidModifier)

