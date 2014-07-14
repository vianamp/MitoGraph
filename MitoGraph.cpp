// ==============================================================
// MitoGraph: Quantifying Mitochondrial Content in Live Cells
// Developed by Matheus P. Viana - vianamp@gmail.com - 2014.05.28
// Susanne Rafelski Lab, University of California Irvine
// Please, check the documentation at ?
// ==============================================================

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
#include <vtkAppendPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>
#include <vtkPoints.h>

#include "MitoThinning.h"

    double _rad = 0.150;
    double _dxy, _dz = -1.0;
    bool _export_graph_files = true;
    bool _scale_polydata_before_save = true;
    bool _export_nodes_label = true;
    double _div_threshold = 0.1666667;

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

/* ================================================================
   IMAGE TRANSFORM
=================================================================*/

// This routine converts 16-bit volumes into 8-bit volumes by
// linearly scaling the original range of intensities [min,max]
// in [0,255] (http://rsbweb.nih.gov/ij/docs/guide/146-28.html)
vtkImageData *Convert16To8bit(vtkImageData *Image);

// Apply a threshold to a ImageData and converts the result in
// a 8-bit ImageData.
vtkImageData *BinarizeAndConvertDoubleToChar(vtkImageData *Image, double threshold);

/* ================================================================
   I/O ROUTINES
=================================================================*/

// Export maximum projection of a given ImageData as a
// PNG file.
void ExportMaxProjection(vtkImageData *Image, const char FileName[], bool binary);

/* ================================================================
   ROUTINES FOR VESSELNESS CALCUATION VIA DISCRETE APPROCH
=================================================================*/

// This routine uses a discrete differential operator to
// calculate the derivatives of a given 3D volume
void GetImageDerivativeDiscrete(vtkDataArray *Image, int *dim, char direction, vtkFloatArray *Derivative);

// This routine calculate the Hessian matrix for each point
// of a 3D volume and its eigenvalues (Discrete Approach)
void GetHessianEigenvaluesDiscrete(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3);

// Calculate the vesselness at each point of a 3D volume based
// based on the Hessian eigenvalues
void GetVesselness(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3);

// Calculate the vesselness over a range of different scales
int MultiscaleVesselness(const char FileName[], double _sigmai, double _dsigma, double _sigmaf, double *attibutes);

/* ================================================================
   DIVERGENCE FILTER
=================================================================*/

// This routine calculates the divergence filter of a 3D volume
// based on the orientation of the gradient vector field
void GetDivergenceFilter(int *Dim, vtkDoubleArray *Scalars);

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

/* ================================================================
   I/O ROUTINES
=================================================================*/

void ExportMaxProjection(vtkImageData *Image, const char FileName[]) {

    #ifdef DEBUG
        printf("Saving Max projection...\n");
    #endif

    int *Dim = Image -> GetDimensions();
    vtkSmartPointer<vtkImageData> MaxP = vtkSmartPointer<vtkImageData>::New();
    MaxP -> SetDimensions(Dim[0],Dim[1],1);
    vtkIdType N = Dim[0] * Dim[1];

    vtkSmartPointer<vtkUnsignedCharArray> MaxPArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
    MaxPArray -> SetNumberOfComponents(1);
    MaxPArray -> SetNumberOfTuples(N);

    int x, y, z;
    double v, vproj;
    for (x = Dim[0]; x--;) {
        for (y = Dim[1]; y--;) {
            vproj = 0;
            for (z = Dim[2]; z--;) {
                v = Image -> GetScalarComponentAsFloat(x,y,z,0);
                vproj = (v > vproj) ? v : vproj;
            }
            MaxPArray -> SetTuple1(MaxP->FindPoint(x,y,0),(unsigned char)vproj);
        }
    }
    MaxPArray -> Modified();

    MaxP -> GetPointData() -> SetScalars(MaxPArray);

    vtkSmartPointer<vtkPNGWriter> PNGWriter = vtkSmartPointer<vtkPNGWriter>::New();
    PNGWriter -> SetFileName(FileName);
    PNGWriter -> SetFileDimensionality(2);
    PNGWriter -> SetCompressionLevel(0);
    PNGWriter -> SetInputData(MaxP);
    PNGWriter -> Write();

    #ifdef DEBUG
        printf("File Saved!\n");
    #endif

}

/* ================================================================
   IMAGE TRANSFORM
=================================================================*/

vtkImageData *Convert16To8bit(vtkImageData *Image) {

    // 8-Bit images
    if (Image -> GetScalarType() == VTK_UNSIGNED_CHAR) {

        return Image;

    // 16-Bit images
    } else if (Image -> GetScalarType() == VTK_UNSIGNED_SHORT) {

        #ifdef DEBUG
            printf("Converting from 16-bit to 8-bit...\n");
        #endif

        vtkImageData *Image8 = vtkImageData::New();
        Image8 -> ShallowCopy(Image);

        vtkDataArray *ScalarsShort = Image -> GetPointData() -> GetScalars();
        unsigned long int N = ScalarsShort -> GetNumberOfTuples();
        double range[2];
        ScalarsShort -> GetRange(range);

        #ifdef DEBUG
            printf("\tOriginal intensities range: [%d-%d]\n",(int)range[0],(int)range[1]);
        #endif

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

vtkImageData *BinarizeAndConvertDoubleToChar(vtkImageData *Image, double threshold) {

    vtkImageData *Image8 = vtkImageData::New();
    Image8 -> ShallowCopy(Image);

    vtkDataArray *ScalarsDouble = Image -> GetPointData() -> GetScalars();
    unsigned long int N = ScalarsDouble -> GetNumberOfTuples();
    double range[2];
    ScalarsDouble -> GetRange(range);

    vtkSmartPointer<vtkUnsignedCharArray> ScalarsChar = vtkSmartPointer<vtkUnsignedCharArray>::New();
    ScalarsChar -> SetNumberOfComponents(1);
    ScalarsChar -> SetNumberOfTuples(N);
        
    double x;
    for ( vtkIdType id = N; id--; ) {
        x = ScalarsDouble -> GetTuple1(id);
        if (x<=threshold) {
            ScalarsChar -> SetTuple1(id,0);
        } else {
            ScalarsChar -> SetTuple1(id,255);
        }
    }
    ScalarsChar -> Modified();

    Image8 -> GetPointData() -> SetScalars(ScalarsChar);
    return Image8;

}

/* ================================================================
   ROUTINES FOR VESSELNESS CALCUATION VIA DISCRETE APPROCH
=================================================================*/

void GetImageDerivativeDiscrete(vtkDataArray *Image, int *dim, char direction, vtkFloatArray *Derivative) {
    #ifdef DEBUG
        printf("Calculating Image Derivatives (Discrete)...\n");
    #endif

    double d, f1, f2;
    vtkIdType i, j, k;
    if (direction=='x') {
        for (i = (vtkIdType)dim[0]; i--;) {
            for (j = (vtkIdType)dim[1]; j--;) {
                for (k = (vtkIdType)dim[2]; k--;) {
                    if (i==0) {
                        Image -> GetTuple(GetId(1,j,k,dim),&f1);
                        Image -> GetTuple(GetId(0,j,k,dim),&f2);
                    }
                    if (i==dim[0]-1) {
                        Image -> GetTuple(GetId(dim[0]-1,j,k,dim),&f1);
                        Image -> GetTuple(GetId(dim[0]-2,j,k,dim),&f2);
                    }
                    if (i>0&i<dim[0]-1) {
                        Image -> GetTuple(GetId(i+1,j,k,dim),&f1);
                        Image -> GetTuple(GetId(i-1,j,k,dim),&f2);
                        f1 /= 2.0; f2 /= 2.0;
                    }
                    d = f1 - f2;
                    Derivative -> SetTuple(GetId(i,j,k,dim),&d);
                }
            }
        }
    }

    if (direction=='y') {
        for (i = (vtkIdType)dim[0]; i--;) {
            for (j = (vtkIdType)dim[1]; j--;) {
                for (k = (vtkIdType)dim[2]; k--;) {
                    if (j==0) {
                        Image -> GetTuple(GetId(i,1,k,dim),&f1);
                        Image -> GetTuple(GetId(i,0,k,dim),&f2);
                    }
                    if (j==dim[1]-1) {
                        Image -> GetTuple(GetId(i,dim[1]-1,k,dim),&f1);
                        Image -> GetTuple(GetId(i,dim[1]-2,k,dim),&f2);
                    }
                    if (j>0&j<dim[1]-1) {
                        Image -> GetTuple(GetId(i,j+1,k,dim),&f1);
                        Image -> GetTuple(GetId(i,j-1,k,dim),&f2);
                        f1 /= 2.0; f2 /= 2.0;
                    }
                    d = f1 - f2;
                    Derivative -> SetTuple(GetId(i,j,k,dim),&d);
                }
            }
        }
    }

    if (direction=='z') {
        for (i = (vtkIdType)dim[0]; i--;) {
            for (j = (vtkIdType)dim[1]; j--;) {
                for (k = (vtkIdType)dim[2]; k--;) {
                    if (k==0) {
                        Image -> GetTuple(GetId(i,j,1,dim),&f1);
                        Image -> GetTuple(GetId(i,j,0,dim),&f2);
                    }
                    if (k==dim[2]-1) {
                        Image -> GetTuple(GetId(i,j,dim[2]-1,dim),&f1);
                        Image -> GetTuple(GetId(i,j,dim[2]-2,dim),&f2);
                    }
                    if (k>0&k<dim[2]-1) {
                        Image -> GetTuple(GetId(i,j,k+1,dim),&f1);
                        Image -> GetTuple(GetId(i,j,k-1,dim),&f2);
                        f1 /= 2.0; f2 /= 2.0;
                    }
                    d = f1 - f2;
                    Derivative -> SetTuple(GetId(i,j,k,dim),&d);
                }
            }
        }
    }

}

void GetHessianEigenvaluesDiscrete(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3) {
    #ifdef DEBUG
        printf("Calculating Hessian Eigeinvalues (Discrete)...\n");
    #endif

    int *Dim = Image -> GetDimensions();
    vtkIdType id, N = Image -> GetNumberOfPoints();
    double H[3][3], Eva[3], Eve[3][3], dxx, dyy, dzz, dxy, dxz, dyz, l1, l2, l3, frobnorm;

    #ifdef DEBUG
        printf("Calculating Gaussian Convolution...\n");
    #endif

    vtkSmartPointer<vtkImageGaussianSmooth> Gauss = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    Gauss -> SetInputData(Image);
    Gauss -> SetDimensionality(3);
    Gauss -> SetRadiusFactors(10,10,10);
    Gauss -> SetStandardDeviations(sigma,sigma,sigma);
    Gauss -> Update();

    vtkFloatArray *ImageG = (vtkFloatArray*) Gauss -> GetOutput() -> GetPointData() -> GetScalars();

    vtkSmartPointer<vtkFloatArray> Dx = vtkSmartPointer<vtkFloatArray>::New(); Dx -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dy = vtkSmartPointer<vtkFloatArray>::New(); Dy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dz = vtkSmartPointer<vtkFloatArray>::New(); Dz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxx = vtkSmartPointer<vtkFloatArray>::New(); Dxx -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dyy = vtkSmartPointer<vtkFloatArray>::New(); Dyy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dzz = vtkSmartPointer<vtkFloatArray>::New(); Dzz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxy = vtkSmartPointer<vtkFloatArray>::New(); Dxy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxz = vtkSmartPointer<vtkFloatArray>::New(); Dxz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dyz = vtkSmartPointer<vtkFloatArray>::New(); Dyz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Fro = vtkSmartPointer<vtkFloatArray>::New();
    Fro -> SetNumberOfComponents(1);
    Fro -> SetNumberOfTuples(N);

    GetImageDerivativeDiscrete(ImageG,Dim,'x',Dx);
    GetImageDerivativeDiscrete(ImageG,Dim,'y',Dy);
    GetImageDerivativeDiscrete(ImageG,Dim,'z',Dz);
    GetImageDerivativeDiscrete(Dx,Dim,'x',Dxx);
    GetImageDerivativeDiscrete(Dy,Dim,'y',Dyy);
    GetImageDerivativeDiscrete(Dz,Dim,'z',Dzz);
    GetImageDerivativeDiscrete(Dy,Dim,'x',Dxy);
    GetImageDerivativeDiscrete(Dz,Dim,'x',Dxz);
    GetImageDerivativeDiscrete(Dz,Dim,'y',Dyz);

    for ( id = N; id--; ) {
        l1 = l2 = l3 = 0.0;
        H[0][0]=Dxx->GetTuple1(id); H[0][1]=Dxy->GetTuple1(id); H[0][2]=Dxz->GetTuple1(id);
        H[1][0]=Dxy->GetTuple1(id); H[1][1]=Dyy->GetTuple1(id); H[1][2]=Dyz->GetTuple1(id);
        H[2][0]=Dxz->GetTuple1(id); H[2][1]=Dyz->GetTuple1(id); H[2][2]=Dzz->GetTuple1(id);
        frobnorm = FrobeniusNorm(H);
        if (H[0][0]+H[1][1]+H[2][2]<0.0) {
            vtkMath::Diagonalize3x3(H,Eva,Eve);
            l1 = Eva[0]; l2 = Eva[1]; l3 = Eva[2];
            Sort(&l1,&l2,&l3);
        }
        L1 -> SetTuple1(id,l1);
        L2 -> SetTuple1(id,l2);
        L3 -> SetTuple1(id,l3);
        Fro -> SetTuple1(id,frobnorm);
    }
    double ftresh,frobenius_norm_range[2];
    Fro -> GetRange(frobenius_norm_range);
    ftresh = sqrt(frobenius_norm_range[1]);

    for ( id = N; id--; ) {
        if ( Fro->GetTuple1(id) < ftresh) {
            L1 -> SetTuple1(id,0.0);
            L2 -> SetTuple1(id,0.0);
            L3 -> SetTuple1(id,0.0);
        }
    }
    L1 -> Modified();
    L2 -> Modified();
    L3 -> Modified();

}

/* ================================================================
   VESSELNESS ROUTINE
=================================================================*/

void GetVesselness(double sigma, vtkImageData *Image, vtkDoubleArray *L1, vtkDoubleArray *L2, vtkDoubleArray *L3) {

    #ifdef DEBUG
        printf("Calculating Vesselness...\n");
    #endif

    double c = 500.0;
    double beta = 0.5;
    double alpha = 0.5;
    double std = 2 * c * c;
    double rbd = 2 * beta * beta;
    double rad = 2 * alpha * alpha;
    double l1, l2, l3, ra, ran, rb, rbn, st, stn, ft_old, ft_new;
    vtkIdType id, N = Image -> GetNumberOfPoints();

    GetHessianEigenvaluesDiscrete(sigma,Image,L1,L2,L3);
    
    for ( id = N; id--; ) {
        l1 = L1 -> GetTuple1(id);
        l2 = L2 -> GetTuple1(id);
        l3 = L3 -> GetTuple1(id);
        if (l2<0&&l3<0) {

            ra = fabs(l2) / fabs(l3);
            ran = -ra * ra;

            rb = fabs(l1) / sqrt(l2*l3);
            rbn = -rb * rb;

            st = sqrt(l1*l1+l2*l2+l3*l3);
            stn = -st * st;

            ft_new = (1-exp(ran/rad)) * exp(rbn/rbd) * (1-exp(stn/std));

            //L1 is used to return vesselness values
            L1 -> SetTuple1(id,ft_new);

        } else L1 -> SetTuple1(id,0.0);
    }
    L1 -> Modified();
}

/* ================================================================
   DIVERGENCE FILTER
=================================================================*/

void GetDivergenceFilter(int *Dim, vtkDoubleArray *Scalars) {

    #ifdef DEBUG
        printf("Calculating Divergent Filter...\n");
    #endif

    vtkIdType id;
    int register j, i;
    int x, y, z, s = 2;
    double v, norm, V[6][3];
    int Dx[6] = {1,-1,0,0,0,0};
    int Dy[6] = {0,0,1,-1,0,0};
    int Dz[6] = {0,0,0,0,1,-1};
    int MI[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

    vtkSmartPointer<vtkDoubleArray> Div = vtkSmartPointer<vtkDoubleArray>::New();
    Div -> SetNumberOfTuples(Scalars->GetNumberOfTuples());

    for (z = s+1; z < Dim[2]-s-1; z++) {
        for (y = s+1; y < Dim[1]-s-1; y++) {
            for (x = s+1; x < Dim[0]-s-1; x++) {
                v = 0;
                id = GetId(x,y,z,Dim);
                if (Scalars->GetTuple1(id)) {
                    for (i = 0; i < 6; i++) {
                        for (j = 0; j < 3; j++) {
                            V[i][j]  = Scalars -> GetTuple1(GetId(x+s*Dx[i]+MI[j][0],y+s*Dy[i]+MI[j][1],z+s*Dz[i]+MI[j][2],Dim));
                            V[i][j] -= Scalars -> GetTuple1(GetId(x+s*Dx[i]-MI[j][0],y+s*Dy[i]-MI[j][1],z+s*Dz[i]-MI[j][2],Dim));
                        }
                        norm = sqrt(pow(V[i][0],2)+pow(V[i][1],2)+pow(V[i][2],2));
                        if (norm) {V[i][0]/=norm; V[i][1]/=norm; V[i][2]/=norm; }
                    }
                    v = (V[0][0]-V[1][0])+(V[2][1]-V[3][1])+(V[4][2]-V[5][2]);
                    v = (v<0) ? -v / 6.0 : 0.0;
                }
                Div -> InsertTuple1(id,v);
            }
        }
    }
    Div -> Modified();
    Scalars -> DeepCopy(Div);
    Scalars -> Modified();

}

/* ================================================================
   MULTISCALE VESSELNESS
=================================================================*/

int MultiscaleVesselness(const char FileName[], double _sigmai, double _dsigma, double _sigmaf, double *attributes) {     

    // Loading multi-paged TIFF file (Supported by VTK 6.2 and higher)
    char _fullpath[256];
    sprintf(_fullpath,"%s.tif",FileName);
    vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    int errlog = TIFFReader -> CanReadFile(_fullpath);
    // File cannot be opened
    if (!errlog) {
        printf("File %s connot be opened.\n",_fullpath);
        return -1;
    }
    TIFFReader -> SetFileName(_fullpath);
    TIFFReader -> Update();

    //DATA CONVERSION
    //---------------

    int *Dim = TIFFReader -> GetOutput() -> GetDimensions();

    #ifdef DEBUG
        printf("MitoGraph V2.0 [DEBUG mode]\n");
        printf("File name: %s\n",_fullpath);
        printf("Volume dimensions: %dx%dx%d\n",Dim[0],Dim[1],Dim[2]);
        printf("Scales to run: [%1.3f:%1.3f:%1.3f]\n",_sigmai,_dsigma,_sigmaf);
        printf("Threshold: %1.5f\n",_div_threshold);
    #endif

    // Conversion 16-bit to 8-bit
    vtkImageData *Image = Convert16To8bit(TIFFReader->GetOutput());

    if (!Image) printf("Format not supported.\n");

    //VESSELNESS
    //----------

    vtkIdType id, N = Image -> GetNumberOfPoints();

    vtkSmartPointer<vtkDoubleArray> AUX1 = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> AUX2 = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> AUX3 = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> VSSS = vtkSmartPointer<vtkDoubleArray>::New();

    AUX1 -> SetNumberOfTuples(N);
    AUX2 -> SetNumberOfTuples(N);
    AUX3 -> SetNumberOfTuples(N);
    VSSS -> SetNumberOfTuples(N);
    VSSS -> FillComponent(0,0);

    double sigma, vn, vo;

    for ( sigma = _sigmai; sigma <= _sigmaf+0.5*_dsigma; sigma += _dsigma ) {
        
        #ifdef DEBUG
            printf("Running sigma = %1.3f\n",sigma);
        #endif
        
        GetVesselness(sigma,Image,AUX1,AUX2,AUX3);
        
        for ( id = N; id--; ) {
            vn = AUX1 -> GetTuple1(id);
            vo = VSSS -> GetTuple1(id);
            if ( vn > vo ) {
                VSSS -> SetTuple1(id,vn);
            }
        }

    }
    VSSS -> Modified();

    //DIVERGENCE FILTER
    //-----------------

    GetDivergenceFilter(Dim,VSSS);

    vtkImageData *ImageEnhanced = vtkImageData::New();
    ImageEnhanced -> GetPointData() -> SetScalars(VSSS);
    ImageEnhanced -> SetDimensions(Dim);

    //CREATING SURFACE POLYDATA
    //-------------------------

    vtkSmartPointer<vtkContourFilter> Filter = vtkSmartPointer<vtkContourFilter>::New();
    Filter -> SetInputData(ImageEnhanced);
    Filter -> SetValue(1,_div_threshold);
    Filter -> Update();

    //SAVING SURFACE
    //--------------

    sprintf(_fullpath,"%s_surface.vtk",FileName);
    SavePolyData(Filter->GetOutput(),_fullpath,true);

    //SAVING BINARY IMAGEDATA
    //-----------------------

    vtkImageData *Binary = BinarizeAndConvertDoubleToChar(ImageEnhanced,_div_threshold);

    sprintf(_fullpath,"%s.png",FileName);
    ExportMaxProjection(Binary,_fullpath);

    Thinning3D(Binary,FileName,attributes);

    return 0;
}

/* ================================================================
   MAIN ROUTINE
=================================================================*/

int main(int argc, char *argv[]) {     

    int i;
    char _prefix[64];
    char _impath[128];
    sprintf(_impath,"");
    double _sigmai = 1.00;
    double _sigmaf = 1.50;
       int _nsigma = 6;

    // Collecting input parameters
    for (i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            sprintf(_impath,"%s//",argv[i+1]);
        }
        if (!strcmp(argv[i],"-xy")) {
            _dxy = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-z")) {
            _dz = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-rad")) {
            _rad = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-scales")) {
            _sigmai = atof(argv[i+1]);
            _sigmaf = atof(argv[i+2]);
            _nsigma = atoi(argv[i+3]);
        }
        if (!strcmp(argv[i],"-threshold")) {
            _div_threshold = (double)atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-scale_off")) {
            _scale_polydata_before_save = false;
        }        
        if (!strcmp(argv[i],"-graph_off")) {
            _export_graph_files = false;
        }
        if (!strcmp(argv[i],"-labels_off")) {
            _export_nodes_label = false;
        }
    }

    if (_dz<0) {
        printf("Please, use -dxy and -dz to provide the pixel size.\n");
        return -1;
    }

    // Generating list of files to run
    char _cmd[256];
    sprintf(_cmd,"ls %s*.tif | sed -e 's/.tif//' > %smitograph.files",_impath,_impath);
    system(_cmd);

    double _dsigma = (_sigmaf-_sigmai) / (_nsigma-1);

    // Generating summary file and writing the header
    char _summaryfilename[256];
    sprintf(_summaryfilename,"%ssummary.txt",_impath);
    FILE *fsummary = fopen(_summaryfilename,"w");
    fprintf(fsummary,"MitoGraph V2.0\n");
    fprintf(fsummary,"Folder: %s\n",_impath);
    fprintf(fsummary,"Pixel size: -xy %1.4fum, -z %1.4fum\n",_dxy,_dz);
    fprintf(fsummary,"Average tubule radius: -r %1.4fum\n",_rad);
    fprintf(fsummary,"Scales: -scales %1.2f",_sigmai);
    for ( double sigma = _sigmai+_dsigma; sigma < _sigmaf+0.5*_dsigma; sigma += _dsigma )
        fprintf(fsummary," %1.2f",sigma);
    fprintf(fsummary,"\nPost-divergence threshold: -threshold %1.5f\n",_div_threshold);
    time_t now = time(0);
    fprintf(fsummary,"%s\n",ctime(&now));
    fprintf(fsummary,"Image path\tsurface-volume (um3)\ttotal length (um)\tskeleton-volume (um3)\n");

    // Multiscale vesselness
    char _tifffilename[256];
    char _tifflistpath[128];
    double *attibutes = new double[3];
    sprintf(_tifflistpath,"%smitograph.files",_impath);
    FILE *f = fopen(_tifflistpath,"r");
    while (fgets(_tifffilename,256, f) != NULL) {
        _tifffilename[strcspn(_tifffilename, "\n" )] = '\0';
        MultiscaleVesselness(_tifffilename,_sigmai,_dsigma,_sigmaf,attibutes);
        
        // Saving network attributes in the group file
        fprintf(fsummary,"%s\t%1.5f\t%1.5f\t%1.5f\n",_tifffilename,attibutes[0],attibutes[1],attibutes[2]);

        // Also printing on the screen
        printf("%s\t%1.5f\t%1.5f\n",_tifffilename,attibutes[0],attibutes[2]);

    }
    fclose(f);
    fclose(fsummary);

}
