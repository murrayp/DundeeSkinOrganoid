
#ifndef TESTSKINORGANOID3D_HPP_
#define TESTSKINORGANOID3D_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "TrianglesMeshReader.hpp"
#include "OffLatticeSimulation.hpp"
#include "TrianglesMeshWriter.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "CellVolumesWriter.hpp"
#include "Debug.hpp"
#include "SkinOrganoidProperty.hpp"
#include "SkinOrganoidLabelWriter.hpp"

#include "CellLabelWriter.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "SkinOrganoidModifier.hpp"

#include "StemCellProliferativeType.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "EllipsoidNodeBasedCellPopulation.hpp"
#include "SutterlinEllipsoidAndBasementMembraneForce.hpp"
#include "EllipsoidNodeAttributes.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellEllipsoidWriter.hpp"
#include "UniformCellCycleModel.hpp"

// Cell population writers
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSkinOrganoid3d : public AbstractCellBasedWithTimingsTestSuite
{
private:
    double mLocationGhosts;
    double mLocationWithoutGhosts;

    MutableMesh<3,3>* Make3dMesh(unsigned width=3, unsigned height=3, unsigned depth=3)
    {
        MutableMesh<3,3>* p_mesh = new MutableMesh<3,3>;
        p_mesh->ConstructCuboid(width, height, depth);
        TrianglesMeshWriter<3,3> mesh_writer("", "3dSpringMesh");
        mesh_writer.WriteFilesUsingMesh(*p_mesh);

        return p_mesh;
    }

public:

    void TestSkinOrgnaoid3D()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        //TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh", "StartMesh");
        //mesh_writer.WriteFilesUsingMesh(mesh);

        PRINT_VARIABLE(mesh.GetNumNodes());
        // Create cells
        //std::vector<CellPtr> cells;
        //CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        //cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());


        // Create cells
      std::vector<CellPtr> cells;
      MAKE_PTR(WildTypeCellMutationState, p_state);
      MAKE_PTR(TransitCellProliferativeType, p_type);
      for (unsigned i=0; i<mesh.GetNumNodes(); i++)
      {
          UniformCellCycleModel* p_model = new UniformCellCycleModel();
          p_model->SetMinCellCycleDuration(1.0);
          p_model->SetMaxCellCycleDuration(1.6);
          CellPtr p_cell(new Cell(p_state, p_model));
          p_cell->SetCellProliferativeType(p_type);


          MAKE_PTR(SkinOrganoidProperty, p_property);
          p_property->SetCellDifferentiatedType(0u);
          p_cell->AddCellProperty(p_property);


          //double birth_time = -RandomNumberGenerator::Instance()->ranf();
          p_cell->SetBirthTime(-0.9);



          //p_label->SetCellDifferentiatedType(0u);
          //p_cell>AddCellProperty(p_label);

          cells.push_back(p_cell);
      }

      MeshBasedCellPopulation<3> cell_population(mesh, cells);
      cell_population.SetWriteVtkAsPoints(true);
      cell_population.AddCellWriter<SkinOrganoidLabelWriter>();
      //cell_population.AddCellWriter<CellLabelWriter>();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSkinOrganoid3D");

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);


        //c_vector<double,3> point = zero_vector<double>(2);
       //       c_vector<double,3> normal = zero_vector<double>(2);

              MAKE_PTR(SkinOrganoidModifier<3>, p_modifier);
              //p_modifier->SetOutputDirectory("TestS");
              simulator.AddSimulationModifier(p_modifier);

         //point(2)=-10.5;
         //normal(2) = -1.0;
         //MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc1, (&cell_population, point, normal)); // y>0
         //simulator.AddCellPopulationBoundaryCondition(p_bc1);
        // Test SetSamplingTimestepMultiple method
        //TS_ASSERT_EQUALS(simulator.mSamplingTimestepMultiple, 1u);
       //imulator.SetSamplingTimestepMultiple(2);
        //TS_ASSERT_EQUALS(simulator.mSamplingTimestepMultiple, 2u);

        // Uncommenting this line calls an error in accessing nodes in the vertex elements #
        //cell_population.AddPopulationWriter<VoronoiDataWriter>();

        simulator.SetEndTime(0.21);
        simulator.Solve();
    }

    void TestSutterlinWithCellDifferentiation()
   {
	   RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

	   // Create a mesh
	   std::vector<Node<3>*> nodes;
	   unsigned index = 0;
	   unsigned cells_across = 5;
	   double scaling = 1.0;
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
		   mesh.GetNode(0u)->rGetNodeAttributes().resize(2);
		   mesh.GetNode(i)->rGetNodeAttributes()[NA_SEMIMAJORAXIS] = 0.55; // micrometres
		   mesh.GetNode(i)->rGetNodeAttributes()[NA_SEMIMINORAXIS] = 0.5; // micrometres
	   }

	   // Create cells
	   std::vector<CellPtr> cells;
	   MAKE_PTR(WildTypeCellMutationState, p_state);
	   MAKE_PTR(StemCellProliferativeType, p_type);
	   for (unsigned i=0; i<mesh.GetNumNodes(); i++)
	   {
		   UniformCellCycleModel* p_model = new UniformCellCycleModel();
		   p_model->SetMinCellCycleDuration(8.0);
		   p_model->SetMaxCellCycleDuration(15.6);
		   CellPtr p_cell(new Cell(p_state, p_model));
		   p_cell->SetCellProliferativeType(p_type);


		   MAKE_PTR(SkinOrganoidProperty, p_property);
		   p_property->SetCellDifferentiatedType(0u);
		   p_property->SetIntraCellularCalcium(0.0);

		   p_cell->AddCellProperty(p_property);


		   double birth_time = -RandomNumberGenerator::Instance()->ranf();
		   p_cell->SetBirthTime(birth_time);



		   //p_label->SetCellDifferentiatedType(0u);
		   //p_cell>AddCellProperty(p_label);

		   cells.push_back(p_cell);
	   }

	   // Create a cell population
	   EllipsoidNodeBasedCellPopulation<3> cell_population(mesh, cells);
	   //cell_population.SetWriteVtkAsPoints(true);
	   cell_population.AddCellWriter<SkinOrganoidLabelWriter>();
	   cell_population.AddCellWriter<CellEllipsoidWriter>();

	   // Create a simulation
	   OffLatticeSimulation<3> simulator(cell_population);
	   simulator.SetOutputDirectory("TestSutterlinEllipsoidCellsahh");
	   simulator.SetSamplingTimestepMultiple(12);
	   simulator.SetEndTime(20.0);

	   // Pass force law to the simulation
	   MAKE_PTR(SutterlinEllipsoidAndBasementMembraneForce<3>, p_force);
	   simulator.AddForce(p_force);

	   MAKE_PTR(SkinOrganoidModifier<3>, p_modifier2);
	   simulator.AddSimulationModifier(p_modifier2);

	   c_vector<double,3> point1 = zero_vector<double>(3);
	   c_vector<double,3> normal1 = zero_vector<double>(3);
	   point1(2)=-0.5;
	   normal1(2) = -1.0;
	   MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc1, (&cell_population, point1, normal1)); // y>0
	   simulator.AddCellPopulationBoundaryCondition(p_bc1);

	   c_vector<double,3> point2 = zero_vector<double>(3);
	   c_vector<double,3> normal2 = zero_vector<double>(3);
	   point2(0)=-0.5;
	   normal2(0) = -1.0;

	   MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc2, (&cell_population, point2, normal2)); // y>0
	   simulator.AddCellPopulationBoundaryCondition(p_bc2);

	   c_vector<double,3> point3 = zero_vector<double>(3);
	   c_vector<double,3> normal3 = zero_vector<double>(3);
	   point3(0)=5.5;
	   normal3(0) = 1.0;

	   MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc3, (&cell_population, point3, normal3)); // y>0
	   simulator.AddCellPopulationBoundaryCondition(p_bc3);
	   /*
	   c_vector<double,3> point4 = zero_vector<double>(3);
	  c_vector<double,3> normal4 = zero_vector<double>(3);
	  point4(1)=-0.5;
	  normal4(1) = -1.0;

	  MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc4, (&cell_population, point4, normal4)); // y>0
	  simulator.AddCellPopulationBoundaryCondition(p_bc4);

	  c_vector<double,3> point5 = zero_vector<double>(3);
	  c_vector<double,3> normal5 = zero_vector<double>(3);
	  point5(1)=5.5;
	  normal5(1) = 1.0;

	  MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc5, (&cell_population, point5, normal5)); // y>0
	  simulator.AddCellPopulationBoundaryCondition(p_bc5);

*/
	   // Run the simulation
	   simulator.Solve();

//           for (EllipsoidNodeBasedCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
//                            cell_iter != simulator.rGetCellPopulation().End();
//                            ++cell_iter)
//            {
//                 CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<SkinOrganoidProperty>();
//                 boost::shared_ptr<SkinOrganoidProperty> p_property = boost::static_pointer_cast<SkinOrganoidProperty>(collection.GetProperty());
//                 //unsigned cell_differentiated_type = p_property->GetCellDifferentiatedType();
//                 //unsigned& r_NumMachineFiresInThisTimeStep = p_property->rGetNumMachineFiresInThisTimeStep();
//
//                 // first stupid model for assigning cell types based on height from basal layer
//                 c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
//
//                 double cell_height= cell_location(2);
//
//
//                 if (cell_height <1.0) //basal
//                 {    //p_property->SetCellDifferentiatedType(0u);
//                 }
//                 else if (cell_height < 10.0) // spinosal
//                 {
//                     TS_ASSERT_EQUALS(p_property->GetCellDifferentiatedType(1u),1u);
//                 }
//            }


	   // Avoid memory leaks
	   for (unsigned i=0; i<nodes.size(); i++)
	   {
		   delete nodes[i];
	   }
   }
};

#endif /*TESTSKINORGANOID3D_HPP_*/
