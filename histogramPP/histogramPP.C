#include "fvCFD.H"
#include "fvc.H"

int main(int argc, char *argv[])
{
  timeSelector::addOptions();
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  instantList timeDirs = timeSelector::select0(runTime, args);

  double v=0;
  forAll(mesh.C(),cell){v += mesh.V()[cell];}

  double binCount[30][56]; 
  memset(binCount, 0.0, sizeof(double)*30*56);  //set all elements to zero for initial values

  double binCountMethane[30][56]; 
  memset(binCountMethane, 0.0, sizeof(double)*30*56);  //set all elements to zero for initial values

  int timeLast;
  cout << "Type last time in ms: ";
  cin >> timeLast;

  int nbin = timeLast/5;

  for (int timeI=5; timeI<timeLast; timeI = timeI +5)
  {
        // read temperature results
        runTime.setTime(timeDirs[timeI], timeI);
        mesh.readUpdate();

        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

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

        volScalarField Ych4
        (
            IOobject
            (
                "CH4",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

	// do bin counting
	int index = timeI/5 - 1;
	int binSelection;
        forAll(mesh.C(), cell){
            binSelection = (T[cell]-950.0)/50;
	    binCount[index][binSelection] += mesh.V()[cell];
	    binCountMethane[index][binSelection] += mesh.V()[cell]*rho[cell]*Ych4[cell];
	}  
	
        Info << "Time = " << runTime.timeName() << endl;
  }

  ofstream myfile;
  myfile.open("results_vol.txt");
  for (int i=0; i<56; i++){
    myfile << "bin " << i << " = "; 
    for (int j=0; j<nbin; j++){
        myfile <<  binCount[j][i] << "," ;
    }
    myfile << '\n';
  }
  myfile.close();

  ofstream myfile1;
  myfile1.open("results_ch4mass.txt");
  for (int i=0; i<56; i++){
    myfile1 << "bin " << i << " = "; 
    for (int j=0; j<nbin; j++){
        myfile1 <<  binCountMethane[j][i] << "," ;
    }
    myfile1 << '\n';
  }
  myfile1.close();

  Info << "Volumen:" << v <<  "m3" << endl;
  Info << "finished" << endl;

  return 0;

}

