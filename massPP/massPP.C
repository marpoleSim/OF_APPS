#include "fvCFD.H"
#include "fvc.H"

int main(int argc, char *argv[])
{
  timeSelector::addOptions();
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  instantList timeDirs = timeSelector::select0(runTime, args);

  // calculate fluid domain
  double v=0;
  forAll(mesh.C(),cell){v += mesh.V()[cell];}

  // annular volume mass initialization
  double reactorMass[100]; 
  double reactorVol[100]; 
  double timeMass[100]; 
  memset(reactorMass, 0.0, sizeof(double)*100);  //set all elements to zero for initial values
  memset(reactorVol, 0.0, sizeof(double)*100);  //set all elements to zero for initial values
  memset(timeMass, 0.0, sizeof(double)*100);     //set all elements to zero for initial values

  int timeLast;
  cout << "Input the last time in ms: ";
  cin >> timeLast;

  float zReactor;
  cout << "Input the z (0.035m, above which is deemed as annular volume): ";
  cin >> zReactor;

  for (int timeI=1; timeI<timeLast; timeI = timeI +1)
  {
        // read density results
        runTime.setTime(timeDirs[timeI], timeI);
        mesh.readUpdate();

        volScalarField rho 
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

	// calculate mass
	int index = timeI-1;
        forAll(mesh.C(), cell){
            if(mesh.C()[cell].component(2) > zReactor) { 
                reactorMass[index] += rho[cell]*mesh.V()[cell];
                reactorVol[index] += mesh.V()[cell];
	    }
	}  

        timeMass[index] = runTime.value();

        Info << "Time = " << runTime.timeName() << endl;
  }

  // write mass change history to result file
  ofstream myfile;
  myfile.open("results.txt");
  for (int i=0; i<100; i++) {
      myfile <<  timeMass[i]  << "," << reactorMass[i] << "," << reactorVol[i] << '\n';
  }
  myfile.close();

  Info << "Volumen:" << v <<  "m3" << endl;
  Info << "finished" << endl;

  return 0;

}

