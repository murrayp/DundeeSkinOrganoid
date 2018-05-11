
#ifndef TESTSUTTERLINELLIPSOIDCELLS_HPP_
#define TESTSUTTERLINELLIPSOIDCELLS_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "EllipsoidNodeBasedCellPopulation.hpp"
#include "SutterlinBasementMembraneForce.hpp"
#include "SutterlinEllipsoidForce.hpp"
#include "EllipsoidNodeAttributes.hpp"
#include "CellEllipsoidWriter.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestSutterlinEllipsoidCells : public AbstractCellBasedTestSuite
{
public:
    void TestSimulation()
    {
    	RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        // Create a mesh
        std::vector<Node<3>*> nodes;
        unsigned index = 0;
        unsigned cells_across = 5;
        double scaling = 1;
        for (unsigned i=0; i<cells_across; i++)
        {
            for (unsigned j=0; j<cells_across; j++)
            {
                for (unsigned k=0; k<cells_across; k++)
                {
                	Node<3>* p_node = new Node<3>(index, false,  (double) i * scaling + 0.05*p_gen->ranf() , (double) j * scaling + 0.05*p_gen->ranf(), (double) k * scaling + 0.05*p_gen->ranf());
                    nodes.push_back(p_node);
                    index++;
                }
            }
        }
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
        	mesh.GetNode(i)->AddNodeAttribute(0.0);
        	mesh.GetNode(i)->rGetNodeAttributes()[NA_SEMIMAJORAXIS] = 5; // micrometres
        	mesh.GetNode(i)->rGetNodeAttributes()[NA_SEMIMINORAXIS] = 2; // micrometres
        }

        // Create a vector of cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_type);

        // Create a cell population
        EllipsoidNodeBasedCellPopulation<3> cell_population(mesh, cells);

        cell_population.AddCellWriter<CellEllipsoidWriter>();

        // Create a simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSutterlinEllipsoidCells");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(30.0);

        // Pass force laws to the simulation
        MAKE_PTR(SutterlinEllipsoidForce<3>, p_ellipsoid_force);
        simulator.AddForce(p_ellipsoid_force);
        MAKE_PTR(SutterlinBasementMembraneForce<3>, p_bm_force);
        simulator.AddForce(p_bm_force);

        // Run the simulation
        simulator.Solve();

        // Avoid memory leaks
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /* TESTSUTTERLINELLIPSOIDCELLS_HPP_ */
