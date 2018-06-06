
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
#include "DifferentiatedCellProliferativeType.hpp"
#include "SkinOrganoidCentreBasedDivisionRule.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "DiffusionForce.hpp"





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

class TestSkinOrganoid3dForExhibition : public AbstractCellBasedWithTimingsTestSuite
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

    void TestSutterlinWithCellDifferentiationForExhibition()
       {
           RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

           c_vector<double, 2> basalCellSemiMajorAndMinorAxis;
           basalCellSemiMajorAndMinorAxis(0)=0.5;
           basalCellSemiMajorAndMinorAxis(1)=0.5;

           /*c_vector<double, 2> spinosalCellSemiMajorAndMinorAxis;
           spinosalCellSemiMajorAndMinorAxis(0)=0.5;
           spinosalCellSemiMajorAndMinorAxis(1)=0.5;

           c_vector<double, 2> granularCellSemiMajorAndMinorAxis;
           granularCellSemiMajorAndMinorAxis(0)=0.5;
           granularCellSemiMajorAndMinorAxis(1)=0.5;
    */

           // Create a mesh
           std::vector<Node<3>*> nodes;
           unsigned index = 0;
           unsigned cells_across = 25;
           unsigned cells_deep = 7;

           unsigned cells_up=1;
           double scaling = basalCellSemiMajorAndMinorAxis(0)*1.0;
           double packing_factor=0.5;
           for (unsigned i=0; i<cells_across; i++)
           {
               for (unsigned j=0; j<cells_deep; j++)
               {
                   for (unsigned k=0; k<cells_up; k++)
                   {
                       Node<3>* p_node = new Node<3>(index, false,  (double) i * scaling + 0.05*p_gen->ranf() , (double) j * scaling + 0.05*p_gen->ranf(), (double) k * scaling + 0.05*p_gen->ranf());
                       nodes.push_back(p_node);
                       index++;
                   }
               }
           }
           NodesOnlyMesh<3> mesh;
           mesh.ConstructNodesWithoutMesh(nodes, 1.0);

           for (unsigned i=0; i<mesh.GetNumNodes(); i++)
           {
               mesh.GetNode(i)->AddNodeAttribute(0.0);
               mesh.GetNode(i)->rGetNodeAttributes().resize(2);
               mesh.GetNode(i)->rGetNodeAttributes()[NA_SEMIMAJORAXIS] = basalCellSemiMajorAndMinorAxis(0); // micrometres
               mesh.GetNode(i)->rGetNodeAttributes()[NA_SEMIMINORAXIS] = basalCellSemiMajorAndMinorAxis(1); // micrometres
           }

           // Create cells
           std::vector<CellPtr> cells;
           MAKE_PTR(WildTypeCellMutationState, p_state);
           MAKE_PTR(StemCellProliferativeType, p_type);
           MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

           for (unsigned i=0; i<mesh.GetNumNodes(); i++)
           {
               UniformCellCycleModel* p_model = new UniformCellCycleModel();
               p_model->SetMinCellCycleDuration(8.0);
               p_model->SetMaxCellCycleDuration(8.5);
               CellPtr p_cell(new Cell(p_state, p_model));


               double height =mesh.GetNode(i)->rGetLocation()[2];

               if (height <1.0*basalCellSemiMajorAndMinorAxis(1)) //basal
               {
                   p_cell->SetCellProliferativeType(p_type);
               }
                   else
               {
                   p_cell->SetCellProliferativeType(p_diff_type);
               }

               MAKE_PTR(SkinOrganoidProperty, p_property);
               p_property->SetCellDifferentiatedType(0u);
               p_property->SetIntraCellularCalcium(0.0);

               p_cell->AddCellProperty(p_property);


               double birth_time = -8.25*RandomNumberGenerator::Instance()->ranf();
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

           boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new SkinOrganoidCentreBasedDivisionRule<3,3>());
           cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

           // Create a simulation
           OffLatticeSimulation<3> simulator(cell_population);
           simulator.SetOutputDirectory("TestSutterlinEllipsoidCellsForExhibition");
           simulator.SetSamplingTimestepMultiple(120);
           simulator.SetDt(1.0/120.0/10.0);

           simulator.SetEndTime(100.0);

           // Pass force law to the simulation
           MAKE_PTR(SutterlinEllipsoidAndBasementMembraneForce<3>, p_force);
           simulator.AddForce(p_force);

//           MAKE_PTR(DiffusionForce<3>, p_force_diff);
//           p_force_diff->SetAbsoluteTemperature(0.150);
//                      simulator.AddForce(p_force_diff);

           MAKE_PTR(SkinOrganoidModifier<3>, p_modifier2);
           simulator.AddSimulationModifier(p_modifier2);

           // Kill all cells moving past z=1;
           MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer,(&cell_population, 4.5*unit_vector<double>(3,2), unit_vector<double>(3,2)));
           simulator.AddCellKiller(p_cell_killer);
           c_vector<double,3> point1 = zero_vector<double>(3);
           c_vector<double,3> normal1 = zero_vector<double>(3);
           point1(2)=-0.05;
           normal1(2) = -1.0;
           MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc1, (&cell_population, point1, normal1)); // y>0
           p_bc1->SetUseJiggledNodesOnPlane(true);

           simulator.AddCellPopulationBoundaryCondition(p_bc1);

           c_vector<double,3> point2 = zero_vector<double>(3);
           c_vector<double,3> normal2 = zero_vector<double>(3);
           point2(0)=-0.05;
           normal2(0) = -1.0;
           MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc2, (&cell_population, point2, normal2)); // y>0
           p_bc2->SetUseJiggledNodesOnPlane(true);

           simulator.AddCellPopulationBoundaryCondition(p_bc2);

           c_vector<double,3> point3 = zero_vector<double>(3);
           c_vector<double,3> normal3 = zero_vector<double>(3);
           point3(0)=scaling*cells_across+packing_factor;
           normal3(0) = 1.0;
           MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc3, (&cell_population, point3, normal3)); // y>0
           p_bc2->SetUseJiggledNodesOnPlane(true);

           simulator.AddCellPopulationBoundaryCondition(p_bc3);

           c_vector<double,3> point4 = zero_vector<double>(3);
          c_vector<double,3> normal4 = zero_vector<double>(3);
          point4(1)=-0.05;
          normal4(1) = -1.0;

          MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc4, (&cell_population, point4, normal4)); // y>0
          p_bc4->SetUseJiggledNodesOnPlane(true);
          simulator.AddCellPopulationBoundaryCondition(p_bc4);

          c_vector<double,3> point5 = zero_vector<double>(3);
          c_vector<double,3> normal5 = zero_vector<double>(3);
          point5(1)=scaling*cells_deep+packing_factor;
          normal5(1) = 1.0;

          MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc5, (&cell_population, point5, normal5)); // y>0
          p_bc5->SetUseJiggledNodesOnPlane(true);

          simulator.AddCellPopulationBoundaryCondition(p_bc5);


           // Run the simulation
           simulator.Solve();

           //TS_ASSERT_EQUALS(simulator.GetNumBirths(),)

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

    void NoTestSutterlinWithCellDifferentiation()
   {
	   RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

	   c_vector<double, 2> basalCellSemiMajorAndMinorAxis;
	   basalCellSemiMajorAndMinorAxis(0)=0.5;
	   basalCellSemiMajorAndMinorAxis(1)=0.5;

	   /*c_vector<double, 2> spinosalCellSemiMajorAndMinorAxis;
	   spinosalCellSemiMajorAndMinorAxis(0)=0.5;
	   spinosalCellSemiMajorAndMinorAxis(1)=0.5;

	   c_vector<double, 2> granularCellSemiMajorAndMinorAxis;
	   granularCellSemiMajorAndMinorAxis(0)=0.5;
	   granularCellSemiMajorAndMinorAxis(1)=0.5;
*/

	   // Create a mesh
	   std::vector<Node<3>*> nodes;
	   unsigned index = 0;
	   unsigned cells_across = 5;
	   unsigned cells_up=1;
	   double scaling = basalCellSemiMajorAndMinorAxis(0)*1.0;
	   for (unsigned i=0; i<cells_across; i++)
	   {
		   for (unsigned j=0; j<cells_across; j++)
		   {
			   for (unsigned k=0; k<cells_up; k++)
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
		   mesh.GetNode(i)->rGetNodeAttributes().resize(2);
		   mesh.GetNode(i)->rGetNodeAttributes()[NA_SEMIMAJORAXIS] = basalCellSemiMajorAndMinorAxis(0); // micrometres
		   mesh.GetNode(i)->rGetNodeAttributes()[NA_SEMIMINORAXIS] = basalCellSemiMajorAndMinorAxis(1); // micrometres
	   }

	   // Create cells
	   std::vector<CellPtr> cells;
	   MAKE_PTR(WildTypeCellMutationState, p_state);
	   MAKE_PTR(StemCellProliferativeType, p_type);
	   MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

	   for (unsigned i=0; i<mesh.GetNumNodes(); i++)
	   {
		   UniformCellCycleModel* p_model = new UniformCellCycleModel();
		   p_model->SetMinCellCycleDuration(8.0);
		   p_model->SetMaxCellCycleDuration(8.5);
		   CellPtr p_cell(new Cell(p_state, p_model));


		   double height =mesh.GetNode(i)->rGetLocation()[2];

		   if (height <1.0*basalCellSemiMajorAndMinorAxis(1)) //basal
		   {
		       p_cell->SetCellProliferativeType(p_type);
		   }
		       else
		   {
			   p_cell->SetCellProliferativeType(p_diff_type);
		   }

		   MAKE_PTR(SkinOrganoidProperty, p_property);
		   p_property->SetCellDifferentiatedType(0u);
		   p_property->SetIntraCellularCalcium(0.0);

		   p_cell->AddCellProperty(p_property);


		   double birth_time = -8.25*RandomNumberGenerator::Instance()->ranf();
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
	   cell_population.SetMeinekeDivisionSeparation(0.1);





       boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new SkinOrganoidCentreBasedDivisionRule<3,3>());
       cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

	   // Create a simulation
	   OffLatticeSimulation<3> simulator(cell_population);
	   simulator.SetOutputDirectory("TestSutterlinEllipsoidCellsahh");
	   simulator.SetSamplingTimestepMultiple(100);
	   simulator.SetDt(1.0/120.0/5.0);

	   simulator.SetEndTime(130.0);

	   // Pass force law to the simulation
	   MAKE_PTR(SutterlinEllipsoidAndBasementMembraneForce<3>, p_force);
	   simulator.AddForce(p_force);

	  // MAKE_PTR(SkinOrganoidModifier<3>, p_modifier2);
	   //simulator.AddSimulationModifier(p_modifier2);



	   c_vector<double,3> point1 = zero_vector<double>(3);
	   c_vector<double,3> normal1 = zero_vector<double>(3);
	   point1(2)=-0.05;
	   normal1(2) = -1.0;
	   MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc1, (&cell_population, point1, normal1)); // y>0
	   simulator.AddCellPopulationBoundaryCondition(p_bc1);

	   c_vector<double,3> point2 = zero_vector<double>(3);
	   c_vector<double,3> normal2 = zero_vector<double>(3);
	   point2(0)=-0.05;
	   normal2(0) = -1.0;

	   MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc2, (&cell_population, point2, normal2)); // y>0
	   simulator.AddCellPopulationBoundaryCondition(p_bc2);

	   c_vector<double,3> point3 = zero_vector<double>(3);
	   c_vector<double,3> normal3 = zero_vector<double>(3);
	   point3(0)=2.5;
	   normal3(0) = 1.0;

	   MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc3, (&cell_population, point3, normal3)); // y>0
	      simulator.AddCellPopulationBoundaryCondition(p_bc3);

	   c_vector<double,3> point4 = zero_vector<double>(3);
	  c_vector<double,3> normal4 = zero_vector<double>(3);
	  point4(1)=-0.05;
	  normal4(1) = -1.0;

	  MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc4, (&cell_population, point4, normal4)); // y>0
	  simulator.AddCellPopulationBoundaryCondition(p_bc4);

	  c_vector<double,3> point5 = zero_vector<double>(3);
	  c_vector<double,3> normal5 = zero_vector<double>(3);
	  point5(1)=2.5;
	  normal5(1) = 1.0;

	  MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc5, (&cell_population, point5, normal5)); // y>0
	  simulator.AddCellPopulationBoundaryCondition(p_bc5);


	   // Run the simulation
	   simulator.Solve();

	   //TS_ASSERT_EQUALS(simulator.GetNumBirths(),)

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
