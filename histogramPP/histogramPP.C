#include "fvCFD.H"
#include "fvc.H"

int main(int argc, char *argv[])
{
  timeSelector::addOptions();
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  instantList timeDirs = timeSelector::select0(runTime, args);

  // calculate total volume
  double v=0;
  forAll(mesh.C(),cell){v += mesh.V()[cell];}

  // bin count, volume in maximum 56 temperature bins and 30 time stamps
  double binCount[30][56]; 
  memset(binCount, 0.0, sizeof(double)*30*56);  //set all elements to zero for initial values

  // bin count, methane mass in maximum 56 temperature bins and 30 time stamps 
  double binCountMethane[30][56]; 
  memset(binCountMethane, 0.0, sizeof(double)*30*56);  //set all elements to zero for initial values

  int timeLast;
  cout << "Type last time in ms: ";
  cin >> timeLast;

  // total number of time stamps
  //int nbin = timeLast/5;
  int nbin = timeLast;

  // number of timestamps can not exceed the limit
  if (nbin >=30)  {nbin=29;}

  // initial time
  //timeIni = 5;
  int timeIni = 1;

  // time step
  //timeStep = 5;
  int timeStep = 1;

  // last time
  timeLast = timeIni + timeStep*nbin;

  // marching time
  for (int timeI=timeIni; timeI<timeLast; timeI = timeI + timeStep)
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
	int index = timeI/timeStep - 1;

	// temperature bin
	int binSelected;
	float tempMin = 1300.0;
	float tempBinSize = 10.0;
	// bin are ,
	// 1. 1200 - 1205
	// 2. 1205 - 1210
	// 3. 1210 - 1215
	// ...

        forAll(mesh.C(), cell){

	    // example: if temperature is 1432K, the bin # is 46, which is from 1430K to 1435K. 
	
            binSelected = (T[cell]-tempMin)/tempBinSize;

	    if (binSelected < 0) { binSelected = 0;}    // there is no bin for temperature below tempMin
	    if (binSelected >= 56) { binSelected = 55;}  // there is no bin for temperature larger than tempLast
	    binCount[index][binSelected] += mesh.V()[cell];      // total volume in this temperature bin
	    binCountMethane[index][binSelected] += mesh.V()[cell]*rho[cell]*Ych4[cell];  // total methane mass in this temperature bin 
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

