// ==================================================================
// Liya Ding. 2016.04.
// This code is developed based on the frame work of MitoGraph. 
// It reuses a majority part of the code of MitoGraph.
// MitoGraph is written by Matheus P. Viana, 
// from Susanne Rafelski Lab, University of California Irvine
// ==================================================================

#include <list>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <vtkMath.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkLongArray.h>
#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkTransform.h>
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkVectorText.h>
#include <vtkInformation.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkInformationVector.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkKdTreePointLocator.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkTIFFWriter.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkDirectory.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtksys/SystemTools.hxx>
#include <vtkImageResample.h>
#include <string>

bool _DebugFlag = true;
double _rad = 0.150;
double _dxy, _dz = -1.0;
bool _export_graph_files = true;
bool _export_image_resampled = true;
bool _adaptive_threshold = false;
bool _scale_polydata_before_save = true;
bool _improve_skeleton_quality = false;
bool _export_nodes_label = true;
double _div_threshold = 0.1666667;
bool _checkonly = false;
bool _width = false;
bool _resampleonly = false;
bool _image_2_surface_flag = false;
bool _filename_specific = false;


    //               |------06------|
    //               |------------------------18------------------------|
    //               |---------------------------------------26----------------------------------|
int ssdx_sort[26] = { 0,-1, 0, 1, 0, 0,-1, 0, 1, 0,-1, 1, 1,-1,-1, 0, 1, 0, -1, 1, 1,-1,-1, 1, 1,-1};
int ssdy_sort[26] = { 0, 0,-1, 0, 1, 0, 0,-1, 0, 1,-1,-1, 1, 1, 0,-1, 0, 1, -1,-1, 1, 1,-1,-1, 1, 1};
int ssdz_sort[26] = {-1, 0, 0, 0, 0, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1, -1,-1,-1,-1, 1, 1, 1, 1};

// In order to acess the voxel (x,y,z) from ImageJ, I should use
// GetId(x,(Dim[1]-1)-y,z,Dim);
// Or change the volume orientation...

/**========================================================
 Auxiliar functions
 =========================================================*/


// This routine returns the x of the id-th point of a 3D volume
// of size Dim[0]xDim[1]xDim[2]
 int  GetX(vtkIdType id, int *Dim);

// This routine returns the y of the id-th point of a 3D volume
// of size Dim[0]xDim[1]xDim[2]
 int  GetY(vtkIdType id, int *Dim);

// This routine returns the z of the id-th point of a 3D volume
// of size Dim[0]xDim[1]xDim[2]
 int  GetZ(vtkIdType id, int *Dim);

// This routine returns the id of a point located at coordinate
// (x,y,z) of a 3D volume of size Dim[0]xDim[1]xDim[2]
 vtkIdType GetId(int x, int y, int z, int *Dim);

// Swap values
 void Swap(double *x, double *y);

// Simple sorting algorithm and the output is such that
// l3 >= l2 >= l3
 void Sort(double *l1, double *l2, double *l3);

// Calculate the Frobenius norm of a given 3x3 matrix
// http://mathworld.wolfram.com/FrobeniusNorm.html
 double FrobeniusNorm(double M[3][3]);

// This routine scales the polydata points to the correct dimension
// given by parameters _dxy and _dz.
 void ScalePolyData(vtkSmartPointer<vtkPolyData> PolyData);
 void PolyData2XYPixelScale(vtkSmartPointer<vtkPolyData> PolyData);

/* ================================================================
   IMAGE TRANSFORM
=================================================================*/

// This routine converts 16-bit volumes into 8-bit volumes by
// linearly scaling the original range of intensities [min,max]
// in [0,255] (http://rsbweb.nih.gov/ij/docs/guide/146-28.html)
   vtkSmartPointer<vtkImageData> Convert16To8bit(vtkSmartPointer<vtkImageData> Image);

// Apply a threshold to a ImageData and converts the result in
// a 8-bit ImageData.
   vtkSmartPointer<vtkImageData> BinarizeAndConvertToDouble(vtkSmartPointer<vtkImageData> Image, double threshold);

// Fill holes in the 3D image
   void FillHoles(vtkSmartPointer<vtkImageData> ImageData);

/* ================================================================
   I/O ROUTINES
=================================================================*/

// Different utilities of this tool, major function
   int utilities(const char FileName[]);

   void SaveImageData(vtkSmartPointer<vtkImageData> Image, const char FileName[], bool _resample = false);
   void SavePolyData(vtkSmartPointer<vtkPolyData> PolyData, const char FileName[], bool scale = _scale_polydata_before_save);

/* ================================================================*/

   void SaveImageData(vtkSmartPointer<vtkImageData> Image, const char FileName[], bool _resample) {

   	if(_DebugFlag){
   		printf("Saving ImageData File...\n");
   	}

   	vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();

   	if (_resample) {
   		if(_DebugFlag){
   			printf("\tResampling data...%f\t%f\n",_dxy,_dz);
   		}
   		vtkSmartPointer<vtkImageResample> Resample = vtkSmartPointer<vtkImageResample>::New();
   		Resample -> SetInterpolationModeToLinear();
   		Resample -> SetDimensionality(3);
   		Resample -> SetInputData(Image);
   		Resample -> SetAxisMagnificationFactor(0,1.0);
   		Resample -> SetAxisMagnificationFactor(1,1.0);
   		Resample -> SetAxisMagnificationFactor(2,_dz/_dxy);
   		Resample -> Update();

   		vtkSmartPointer<vtkImageData> ImageResampled = Resample -> GetOutput();
   		ImageResampled -> SetSpacing(1,1,1);

   		writer -> SetInputData(ImageResampled);
   	} else {

   		writer -> SetInputData(Image);

   	}

   	writer -> SetFileType(VTK_BINARY);
   	writer -> SetFileName(FileName);
   	writer -> Write();

   	if(_DebugFlag){
   		printf("\tFile Saved!\n");
   	}
   }

   void SavePolyData(vtkSmartPointer<vtkPolyData> PolyData, const char FileName[], bool scale) {

   	if(_DebugFlag){
   		printf("Saving PolyData from XYZ list...\n");
   	}

   	if(_DebugFlag){
   		printf("\t#Points in PolyData file: %llu.\n",(vtkIdType)PolyData->GetNumberOfPoints());
   	}

   	vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
   	Writer -> SetFileType(VTK_BINARY);
   	Writer -> SetFileName(FileName);
   	Writer -> SetInputData(PolyData);
   	Writer -> Write();

   	if(_DebugFlag){
   		printf("\tFile Saved!\n");
   	}
   }

/**========================================================
 Auxiliar functions
 =========================================================*/

 int GetX(vtkIdType id, int *Dim) {
 	return (int) ( id % (vtkIdType)Dim[0] );
 }

 int GetY(vtkIdType id, int *Dim) {
 	return (int) (  ( id % (vtkIdType)(Dim[0]*Dim[1]) ) / (vtkIdType)Dim[0]  );
 }

 int GetZ(vtkIdType id, int *Dim) {
 	return (int) ( id / (vtkIdType)(Dim[0]*Dim[1]) );
 }

 vtkIdType GetId(int x, int y, int z, int *Dim) {
 	return (vtkIdType)(x+y*Dim[0]+z*Dim[0]*Dim[1]);
 }

 vtkIdType GetReflectedId(int x, int y, int z, int *Dim) {
 	int rx = ceil(0.5*Dim[0]);
 	int ry = ceil(0.5*Dim[1]);
 	int rz = ceil(0.5*Dim[2]);
 	int sx = (x-(rx-0.5)<0) ? -rx : Dim[0]-rx;
 	int sy = (y-(ry-0.5)<0) ? -ry : Dim[1]-ry;
 	int sz = (z-(rz-0.5)<0) ? -rz : Dim[2]-rz;
 	return GetId(x-sx,y-sy,z-sz,Dim);
 }

 void Swap(double *x, double *y) {
 	double t = *y; *y = *x; *x = t;
 }

 void Sort(double *l1, double *l2, double *l3) {
 	if (fabs(*l1) > fabs(*l2)) Swap(l1,l2);
 	if (fabs(*l2) > fabs(*l3)) Swap(l2,l3);
 	if (fabs(*l1) > fabs(*l2)) Swap(l1,l2);
 }

 double FrobeniusNorm(double M[3][3]) {
 	double f = 0.0;
 	for (int i = 3;i--;)
 		for (int j = 3;j--;)
 			f += M[i][j]*M[i][j];
 		return sqrt(f);
 	}

 	void ScalePolyData(vtkSmartPointer<vtkPolyData> PolyData) {
 		double r[3];
 		vtkPoints *Points = PolyData -> GetPoints();
 		for (vtkIdType id = 0; id < Points -> GetNumberOfPoints(); id++) {
 			Points -> GetPoint(id,r);
 			Points -> SetPoint(id,_dxy*r[0],_dxy*r[1],_dz*r[2]);
 		}
 		Points -> Modified();
 	}


 	void PolyData2XYPixelScale(vtkSmartPointer<vtkPolyData> PolyData) {
 		double r[3];
 		vtkPoints *Points = PolyData -> GetPoints();
 		for (vtkIdType id = 0; id < Points -> GetNumberOfPoints(); id++) {
 			Points -> GetPoint(id,r);
 			Points -> SetPoint(id,r[0],r[1],r[2]*_dz/_dxy);
 		}
 		Points -> Modified();
 	}


/* ================================================================
   I/O ROUTINES
=================================================================*/

/* ================================================================
   IMAGE TRANSFORM
=================================================================*/

   vtkSmartPointer<vtkImageData> Convert16To8bit(vtkSmartPointer<vtkImageData> Image) {

    // 8-Bit images
   	if (Image -> GetScalarType() == VTK_UNSIGNED_CHAR) {

   		return Image;

    // 16-Bit images
   	} else if (Image -> GetScalarType() == VTK_UNSIGNED_SHORT) {

   		if(_DebugFlag){
   			printf("Converting from 16-bit to 8-bit...\n");
   		}

   		vtkSmartPointer<vtkImageData> Image8 = vtkImageData::New();
   		Image8 -> ShallowCopy(Image);

   		vtkDataArray *ScalarsShort = Image -> GetPointData() -> GetScalars();
   		unsigned long int N = ScalarsShort -> GetNumberOfTuples();
   		double range[2];
   		ScalarsShort -> GetRange(range);

   		if(_DebugFlag){
   			printf("\tOriginal intensities range: [%d-%d]\n",(int)range[0],(int)range[1]);
   		}

   		vtkSmartPointer<vtkUnsignedCharArray> ScalarsChar = vtkSmartPointer<vtkUnsignedCharArray>::New();
   		ScalarsChar -> SetNumberOfComponents(1);
   		ScalarsChar -> SetNumberOfTuples(N);

   		double x, y;
   		vtkIdType register id;
   		for ( id = N; id--; ) {
   			x = ScalarsShort -> GetTuple1(id);
   			y = 255.0 * (x-range[0]) / (range[1]-range[0]);
   			ScalarsChar -> SetTuple1(id,(unsigned char)y);
   		}
   		ScalarsChar -> Modified();

   		Image8 -> GetPointData() -> SetScalars(ScalarsChar);
   		return Image8;

    // Other depth
   	} else {
   		return NULL;
   	}
   }

   vtkSmartPointer<vtkImageData> BinarizeAndConvertToDouble(vtkSmartPointer<vtkImageData> Image, double threshold) {

   	vtkSmartPointer<vtkImageData> Image16 = vtkImageData::New();
   	Image16 -> ShallowCopy(Image);

   	vtkDataArray *ScalarsDouble = Image -> GetPointData() -> GetScalars();
   	unsigned long int N = ScalarsDouble -> GetNumberOfTuples();
   	double range[2];
   	ScalarsDouble -> GetRange(range);

   	vtkSmartPointer<vtkUnsignedCharArray> ScalarsChar = vtkSmartPointer<vtkUnsignedCharArray>::New();
   	ScalarsChar -> SetNumberOfComponents(1);
   	ScalarsChar -> SetNumberOfTuples(N);

   	double x;



   	if (threshold > 0) {
   		for ( vtkIdType id = N; id--; ) {
   			x = ScalarsDouble -> GetTuple1(id);
   			if (x< threshold-0.1) {
   				ScalarsChar -> SetTuple1(id,0);
   			} else {
   				if(x< threshold+0.1)
   				{
   					ScalarsChar -> SetTuple1(id,40000);
   				}
   				else
   				{
   					ScalarsChar -> SetTuple1(id,0);
   				}
   			}
   		}
   	} else {
   		for ( vtkIdType id = N; id--; ) {
   			x = ScalarsDouble -> GetTuple1(id);
   			ScalarsChar -> SetTuple1(id,(int)(40000*(x-range[0])/(range[1]-range[0])));
   		}        
   	}
   	ScalarsChar -> Modified();

   	Image16 -> GetPointData() -> SetScalars(ScalarsChar);
   	return Image16;

   }

   void FillHoles(vtkSmartPointer<vtkImageData> ImageData) {

   	if(_DebugFlag){
   		printf("\tSearching for holes in the image...\n");
   	}

   	int *Dim = ImageData -> GetDimensions();
   	vtkIdType N = ImageData -> GetNumberOfPoints();

   	vtkIdType i, s, ido, id;

   	int x, y, z;
   	double v, r[3];
   	bool find = true;
   	long long int ro[3];
   	long int scluster, label;
   	ro[0] = Dim[0] * Dim[1] * Dim[2];
   	ro[1] = Dim[0] * Dim[1];

   	vtkSmartPointer<vtkIdList> CurrA = vtkSmartPointer<vtkIdList>::New();
   	vtkSmartPointer<vtkIdList> NextA = vtkSmartPointer<vtkIdList>::New();
   	vtkSmartPointer<vtkLongArray> CSz = vtkSmartPointer<vtkLongArray>::New();
   	vtkSmartPointer<vtkLongArray> Volume = vtkSmartPointer<vtkLongArray>::New();
   	Volume -> SetNumberOfComponents(1);
   	Volume -> SetNumberOfTuples(N);
   	Volume -> FillComponent(0,0);

   	for (x = 1; x < Dim[0]-1; x++) {
   		for (y = 1; y < Dim[1]-1; y++) {
   			for (z = 1; z < Dim[2]-1; z++) {
   				id = ImageData -> FindPoint(x,y,z);
   				if ((unsigned short int)ImageData->GetScalarComponentAsDouble(x,y,z,0)) {
   					Volume -> SetTuple1(id,0);
   				} else {
   					Volume -> SetTuple1(id,1);
   				}
   			}
   		}
   	}

   	Volume -> Modified();

   	label = 0;
   	while (find) {
   		for (s = 0; s < CurrA->GetNumberOfIds(); s++) {
   			ido = CurrA -> GetId(s);
   			ImageData -> GetPoint(ido,r);
   			x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
   			for (i = 0; i < 6; i++) {
   				id = ImageData -> FindPoint(x+ssdx_sort[i],y+ssdy_sort[i],z+ssdz_sort[i]);
   				v = Volume -> GetTuple1(id);
   				if ((long int)v > 0) {
   					NextA -> InsertNextId(id);
   					Volume -> SetTuple1(id,-label);
   					scluster++;
   				}
   			}
   		}
   		if (!NextA->GetNumberOfIds()) {
   			find = false;
   			for (id=ro[0]; id--;) {
   				v = Volume -> GetTuple1(id);
   				if ((long int)v > 0) {
   					find = true;
   					ro[0] = id;
   					break;
   				}
   			}
   			if (label) {
   				CSz -> InsertNextTuple1(scluster);
   			}
   			if (find) {
   				label++;
   				scluster = 1;
   				Volume -> SetTuple1(id,-label);
   				CurrA -> InsertNextId(id);
   			}
   		} else {
   			CurrA -> Reset();
   			CurrA -> DeepCopy(NextA);
   			NextA -> Reset();
   		}
   	}

   	for (id = N; id--;) {
   		if ((long int)Volume->GetTuple1(id)<-1) {
   			ImageData -> GetPointData() -> GetScalars() -> SetTuple1(id,255);
   		}
   	}
   	ImageData -> GetPointData() -> GetScalars() -> Modified();

   	if(_DebugFlag){
   		printf("\tNumber of filled holes: %ld\n",(long int)CSz->GetNumberOfTuples()-1);
   	}

   }

/* ================================================================
=================================================================*/

   int utilities(const char FileName[]) {     

    // Loading multi-paged TIFF file (Supported by VTK 6.2 and higher)
   	char _fullpath[256];
   	sprintf(_fullpath,"%s",FileName);

   	vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
   	int errlog = TIFFReader -> CanReadFile(_fullpath);
    // File cannot be opened
   	if (!errlog) {
   		printf("File %s cannnot be opened.\n",_fullpath);
   		return -1;
   	}
   	TIFFReader -> SetFileName(_fullpath);
   	TIFFReader -> Update();

   	int *Dim = TIFFReader -> GetOutput() -> GetDimensions();

   	if(_DebugFlag){
   		printf("Segmentation results to paraview volume and surfaces V1.0 [DEBUG mode]\n");
   		printf("File name: %s\n",_fullpath);
   		printf("Volume dimensions: %dx%dx%d\n",Dim[0],Dim[1],Dim[2]);
   	}

    // Exporting resampled images
   	if (_export_image_resampled) {
   		sprintf(_fullpath,"%s_resampled.vtk",FileName);
   		SaveImageData(TIFFReader->GetOutput(),_fullpath,true);
   	}


   	if(_image_2_surface_flag){
   		vtkSmartPointer<vtkPolyData> AllNucContours= vtkSmartPointer<vtkPolyData>::New();


   		vtkSmartPointer<vtkAppendPolyData> appendFilter =
   		vtkSmartPointer<vtkAppendPolyData>::New();

   		for( int objectNumber = 1; objectNumber<= TIFFReader-> GetOutput() -> GetScalarRange()[1]; objectNumber++)
   		{

   			vtkSmartPointer<vtkImageData> Binary = BinarizeAndConvertToDouble(TIFFReader->GetOutput(),objectNumber);


	        //FillHoles(Binary);

   			vtkSmartPointer<vtkImageGaussianSmooth> SmoothImage = vtkSmartPointer<vtkImageGaussianSmooth>::New();
   			SmoothImage -> SetInputData(Binary);        
   			SmoothImage -> SetStandardDeviation(2.5);
   			SmoothImage -> Update();

	        // get smooth surface
   			vtkSmartPointer<vtkContourFilter> Filter = vtkSmartPointer<vtkContourFilter>::New();
   			Filter -> SetInputData(SmoothImage->GetOutput());

   			Filter -> SetValue(1,5);
   			Filter -> Update();

   			PolyData2XYPixelScale(Filter-> GetOutput() );


   			vtkSmartPointer<vtkIntArray> ID = vtkSmartPointer<vtkIntArray>::New();
   			ID->SetNumberOfValues(Filter -> GetOutput() -> GetNumberOfPoints() );

	        // Set the ID
   			for(int ii = 0; ii<Filter -> GetOutput() ->GetNumberOfPoints(); ii++ )
   			{
   				ID->SetValue(ii,objectNumber);
   			}

   			Filter -> GetOutput() -> GetPointData() -> SetScalars(ID);

	        // append to the group of polydata surfaces
   			vtkSmartPointer<vtkPolyData> input1 =
   			vtkSmartPointer<vtkPolyData>::New();


   			input1->ShallowCopy(Filter->GetOutput());


   			appendFilter->AddInputData(input1);
   			appendFilter->Update();
   		}
   		sprintf(_fullpath,"%s_smooth_surface_all.vtk",FileName);
   		SavePolyData(appendFilter->GetOutput(),_fullpath);
   	} 
   	return 0;
   }

/* ================================================================
   MAIN ROUTINE
=================================================================*/

   int main(int argc, char *argv[]) {     

   	int i;
   	char _impath[256];
   	sprintf(_impath,"./examples/");
   	char _filename[256];
   	sprintf(_filename,"0");

    // Collecting input parameters
   	for (i = 0; i < argc; i++) {
   		printf("This argument: %s \n",argv[i]);

   		if (!strcmp(argv[i],"-file")) {

   			sprintf(_filename,"%s",argv[i+1]);
   			_filename_specific = true;
   			printf("Specified one file %s",_filename);
   		}

   		if (!strcmp(argv[i],"-path")) {

   			sprintf(_impath,"%s",argv[i+1]);

   			std::string fileString = _impath;

            // if not ended with /, add /
   			if (!fileString.empty()){
   				char lastChar = *fileString.rbegin();
   				if( lastChar != '/'){
   					sprintf(_impath,"%s/",argv[i+1]);
   				}                  
   			}


   		}
   		if (!strcmp(argv[i],"-xy")) {
   			_dxy = atof(argv[i+1]);
   		}
   		if (!strcmp(argv[i],"-z")) {
   			_dz = atof(argv[i+1]);
   		}
   		if (!strcmp(argv[i],"-export_image_resampled")) {
   			_export_image_resampled = true;
   		}
   		if (!strcmp(argv[i],"-resampleonly")) {
   			_resampleonly = true;
   			_export_image_resampled = true;
   		}
   		if (!strcmp(argv[i],"-image_2_surface_flag")) {
   			_image_2_surface_flag = true;
   		}     
   		if (!strcmp(argv[i],"-DebugFlag")) {
   			_DebugFlag = true;
   		}      

   	}


   	if (_dz<0) {        
   		_dxy = 0.09;
   		_dz = 0.2;
   		_export_image_resampled = true;
   		_improve_skeleton_quality = true;
   		_DebugFlag = true;
   		_div_threshold = 0.4;
   		printf("Please, use -dxy and -dz to provide the pixel size.\n");
        //return -1;
   	}

   	vtkSmartPointer<vtkDirectory> directory = vtkSmartPointer<vtkDirectory>::New();  
   	int opened = directory->Open(_impath); 
   	if(!opened)
   	{
   		std::cout << "Invalid directory!" << std::endl;
   		return EXIT_FAILURE;
   	}

   	vtkSmartPointer<vtkLongArray>  TiffFilenameIndexList = vtkSmartPointer<vtkLongArray>::New();

   	if(_filename_specific)
   	{  
   		int numberOfFiles = 1;  
   		std::cout << "1 file " << numberOfFiles << std::endl;

   		TiffFilenameIndexList -> SetNumberOfComponents(1);
   		TiffFilenameIndexList -> SetNumberOfTuples(0);
   		std::string fileString = _impath;

   		if (!fileString.empty()){
   			char lastChar = *fileString.rbegin();
   			if( lastChar != '/'){
   				fileString += "/";
   			}
   		}
   		
   		fileString += _filename;    

   		TiffFilenameIndexList->InsertNextTuple1(1);  

   	}
   	else
   	{
   		int numberOfFiles = directory->GetNumberOfFiles();  
   		std::cout << "Number of all files: " << numberOfFiles << std::endl;

   		TiffFilenameIndexList -> SetNumberOfComponents(1);
   		TiffFilenameIndexList -> SetNumberOfTuples(0);

   		for (int i = 0; i < numberOfFiles; i++)
   		{
   			std::string fileString = _impath;

   			if (!fileString.empty()){
   				char lastChar = *fileString.rbegin();
   				if( lastChar != '/'){
   					fileString += "/";
   				}
   			}


   			fileString += directory->GetFile(i); 
   			std::string ext = vtksys::SystemTools::GetFilenameLastExtension(fileString);    


   			if( (ext == ".tif")||(ext == ".tiff")) {

   				std::cout << "Find tiff file(s): "  << std::endl;
   				std::cout << fileString << std::endl;

   				TiffFilenameIndexList->InsertNextTuple1(i);  
   			}
   		}  

   	}
    //

   	char _tifffilename[512];


   	for(int iTiffInd = 0 ; iTiffInd< TiffFilenameIndexList->GetNumberOfTuples();iTiffInd++) {
   		std::string stdstr_tifffilename = _impath;
   		if (!stdstr_tifffilename.empty()){
   			char lastChar = *stdstr_tifffilename.rbegin();
   			if( lastChar != '/'){
   				stdstr_tifffilename += "/";
   			}
   		}

   		if(_filename_specific)
   		{ 
   			stdstr_tifffilename += _filename;  
   			std::cout << "Work on tiff file: "  << std::endl;
   			std::cout << stdstr_tifffilename << std::endl;            
   			sprintf(_tifffilename,stdstr_tifffilename.c_str());

   		}   
   		else
   		{
   			stdstr_tifffilename += directory->GetFile(TiffFilenameIndexList->GetTuple1(iTiffInd));  
   			std::cout << "Work on tiff file: "  << std::endl;
   			std::cout << stdstr_tifffilename << std::endl;

   			sprintf(_tifffilename,stdstr_tifffilename.c_str());
   		}

    	// the core function for the tool
   		utilities(_tifffilename);

   		printf("%s\t[done]\n",_tifffilename);
   	}

   	return 0;
}
