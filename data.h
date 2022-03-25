//
//						 data.h
//

/************************************************************************
			DECLARATION OF THE CLASS DATA
************************************************************************/
// Class for Data
// A data has 5 fields: 
//       n (of type int), the number of data points 
//       X and Y (arrays of type reals), arrays of x- and y-coordinates
//       meanX and meanY (of type reals), coordinates of the centroid (x and y sample means)

class Data
{
public:

	int n;
	/* reals *X; */
	/* reals *Y; */
	/* reals *err; */
	reals meanX, meanY;
	
	vector<reals> X, Y, err;
	
	// constructors
	Data();
	Data(int N);
	Data(int N, reals X[], reals Y[], reals err[]);
	// Data(int N, reals X[], reals Y[]);
	
	// routines
	void means(void);
	void center(void);
	void scale(void);
	void print(void);
	
	// destructors
	~Data();
};


/************************************************************************
			BODY OF THE MEMBER ROUTINES
************************************************************************/
// Default constructor
Data::Data()
{
	n=0;
	/* X = new reals[n]; */
	/* Y = new reals[n]; */
	/* err = new reals[n]; */
	/* for (int i=0; i<n; i++) */
	/* { */
	/* 	X[i]=0.; */
	/* 	Y[i]=0.; */
	/* 	err[i]=1.; */
	/* } */

	X.clear();
	Y.clear();
	err.clear();
	for (int i=0; i<n; i++) {
	  X.push_back(0);
	  Y.push_back(0);
	  err.push_back(1);
	}
}

// Constructor with assignment of the field N
Data::Data(int N)
{
	n=N;
	/* X = new reals[n]; */
	/* Y = new reals[n]; */
	/* err = new reals[n]; */
	/* for (int i=0; i<n; i++) */
	/* { */
	/* 	X[i]=0.; */
	/* 	Y[i]=0.; */
	/* 	err[i]=1.; */
	/* } */
	
	X.clear();
	Y.clear();
	err.clear();
	for (int i=0; i<n; i++) {
	  X.push_back(0);
	  Y.push_back(0);
	  err.push_back(1);
	}
}

// Constructor with assignment of each field
Data::Data(int N, reals dataX[], reals dataY[], reals dataEY[])
{
	n=N;
	/* X = new reals[n]; */
	/* Y = new reals[n]; */
	/* err = new reals[n]; */
	/* for (int i=0; i<n; i++) */
	/* { */
	/* 	X[i]=dataX[i]; */
	/* 	Y[i]=dataY[i]; */
	/* 	err[i]=dataEY[i]; */
	/* } */
	
	X.clear();
	Y.clear();
	err.clear();
	for (int i=0; i<n; i++) {
	  X.push_back(dataX[i]);
	  Y.push_back(dataY[i]);
	  err.push_back(dataEY[i]);
	}
	
}

// Routine that computes the x- and y- sample means (the coordinates of the centeroid)

void Data::means(void)
{
	meanX=0.; meanY=0.;
	
	for (int i=0; i<n; i++)
	{
		meanX += X[i];
		meanY += Y[i];
	}
	meanX /= n;
	meanY /= n;
}

// Routine that centers the data set (shifts the coordinates to the centeroid)

void Data::center(void)
{
	reals sX=0.,sY=0.;  
	int i;
	
	for (i=0; i<n; i++)
	{
		sX += X[i];
		sY += Y[i];
	}
	sX /= n;
	sY /= n;
	
	for (i=0; i<n; i++)
	{
		X[i] -= sX;
		Y[i] -= sY;
	}
	meanX = 0.;
	meanY = 0.;
}

// Routine that scales the coordinates (makes them of order one)

void Data::scale(void)
{
	reals sXX=0.,sYY=0.,scaling;  
	int i;
	
	for (i=0; i<n; i++)
	{
		sXX += X[i]*X[i];
		sYY += Y[i]*Y[i];
	}
	scaling = sqrt((sXX+sYY)/n/Two);
	
	for (i=0; i<n; i++)
	{
		X[i] /= scaling;
		Y[i] /= scaling;
	}
}

// Printing routine

void Data::print(void)
{
	cout << endl << "The data set has " << n << " points with coordinates :"<< endl;
	
	for (int i=0; i<n-1; i++) cout << setprecision(7) << "(" << X[i] << ","<< Y[i] << "), ";
	
	cout << "(" << X[n-1] << ","<< Y[n-1] << ")\n";
}

// Destructor
Data::~Data()
{
  // cout << " dest called " << endl; 

  X.clear();
  Y.clear();
  err.clear();
      
  /* delete[] X; */
  /* delete[] Y; */
  /* delete[] err; */
  
  /* for (int i=0; i<n; i++) { */
  /*   if(X[i]) delete X[i]; */
  /*   if(Y[i]) delete Y[i]; */
  /*   if(err[i]) delete err[i]; */
  /* } */
  
  /* if(X) delete X; */
  /* if(Y) delete Y; */
  /* if(err) delete err; */
  
 
}


